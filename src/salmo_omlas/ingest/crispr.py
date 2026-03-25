"""Fetch exon/CDS sequence and persist SpCas9 guide candidates."""

from __future__ import annotations

import re
import sqlite3
from dataclasses import dataclass

import requests

from salmo_omlas.config import ENSEMBL_REST
from salmo_omlas.ingest._util import now_iso


@dataclass
class GuideCandidate:
    guide_sequence: str
    pam: str
    strand: str
    start: int
    end: int
    on_target_score: float
    off_target_risk: float
    exon_rank: int | None


def fetch_cds_sequence(ensembl_gene_id: str) -> tuple[str, str]:
    """Return (cds_sequence, molecule_id used) from Ensembl REST."""
    headers = {"Content-Type": "application/json"}
    # Lookup gene -> canonical transcript
    url = f"{ENSEMBL_REST}/lookup/id/{ensembl_gene_id}?expand=1"
    r = requests.get(url, headers=headers, timeout=60)
    r.raise_for_status()
    data = r.json()
    tid = data.get("canonical_transcript")
    if not tid:
        trs = data.get("Transcript") or []
        for tr in trs:
            if tr.get("is_canonical") == 1:
                tid = tr.get("id")
                break
        if not tid and trs:
            tid = trs[0].get("id")
    if not tid:
        raise ValueError("No transcript for gene")
    surl = f"{ENSEMBL_REST}/sequence/id/{tid}?type=cds"
    sr = requests.get(surl, headers=headers, timeout=60)
    sr.raise_for_status()
    sdata = sr.json()
    seq = str(sdata.get("seq", "")).upper()
    return seq, tid


def _gc_score(seq: str) -> float:
    if not seq:
        return 0.0
    g = seq.count("G") + seq.count("C")
    return g / len(seq)


def _offtarget_heuristic(seq: str, guide: str) -> float:
    """Higher = more risky (rough duplicate k-mer count in CDS)."""
    if len(guide) < 12:
        return 0.0
    seed = guide[-12:]
    return float(max(0, seq.count(seed) - 1))


def design_spcas9_guides(
    cds: str,
    *,
    max_guides: int = 20,
) -> list[GuideCandidate]:
    """Scan CDS for 20nt + NGG on forward strand (demo simplification)."""
    cds = cds.upper()
    guides: list[GuideCandidate] = []
    for i in range(0, len(cds) - 23):
        g = cds[i : i + 20]
        pam = cds[i + 20 : i + 23]
        if len(pam) < 3 or pam[1:] != "GG":
            continue
        if not re.fullmatch(r"[ATGC]+", g):
            continue
        on_score = 0.4 * _gc_score(g) + 0.3 * (1.0 - _offtarget_heuristic(cds, g) * 0.05)
        on_score = max(0.0, min(1.0, on_score))
        risk = min(1.0, _offtarget_heuristic(cds, g) * 0.15)
        guides.append(
            GuideCandidate(
                guide_sequence=g,
                pam=pam,
                strand="+",
                start=i + 1,
                end=i + 23,
                on_target_score=round(on_score, 4),
                off_target_risk=round(risk, 4),
                exon_rank=None,
            )
        )
        if len(guides) >= max_guides:
            break
    guides.sort(key=lambda x: (-x.on_target_score, x.off_target_risk))
    return guides[:max_guides]


def save_guides_for_gene(
    conn: sqlite3.Connection,
    ensembl_gene_id: str,
    *,
    max_guides: int = 15,
) -> int:
    cur = conn.cursor()
    row = cur.execute("SELECT id FROM genes WHERE ensembl_gene_id=?", (ensembl_gene_id,)).fetchone()
    if not row:
        raise ValueError(f"Unknown gene {ensembl_gene_id}")
    gid = row[0]
    cds, _tid = fetch_cds_sequence(ensembl_gene_id)
    guides = design_spcas9_guides(cds, max_guides=max_guides)
    ts = now_iso()
    n = 0
    for g in guides:
        cur.execute(
            """INSERT OR REPLACE INTO crispr_targets (
                 gene_id, guide_sequence, pam, strand, chromosome, start, end,
                 on_target_score, off_target_risk, exon_rank, notes, computed_at
               ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                gid,
                g.guide_sequence,
                g.pam,
                g.strand,
                None,
                g.start,
                g.end,
                g.on_target_score,
                g.off_target_risk,
                g.exon_rank,
                "SpCas9 NGG scan on canonical CDS",
                ts,
            ),
        )
        n += 1
    conn.commit()
    return n
