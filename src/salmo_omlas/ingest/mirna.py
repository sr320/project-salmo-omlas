"""miRNA genes and target edges (miRBase GFF + heuristic targets)."""

from __future__ import annotations

import re
import sqlite3

import requests

from salmo_omlas.config import MIRBASE_FTP_GFF3
from salmo_omlas.ingest._util import dumps_metadata, now_iso


def _parse_gff_attributes(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in attr.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def fetch_mirbase_gff_text() -> str:
    r = requests.get(MIRBASE_FTP_GFF3, timeout=120)
    r.raise_for_status()
    return r.text


def load_mirna_features(conn: sqlite3.Connection) -> int:
    """Parse miRBase ssa.gff3 for pre-miRNA / miRNA loci."""
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        ("miRBase", "CURRENT", MIRBASE_FTP_GFF3, now_iso()),
    )
    conn.commit()
    sid = cur.execute(
        "SELECT id FROM sources WHERE name=? ORDER BY id DESC LIMIT 1",
        ("miRBase",),
    ).fetchone()[0]

    try:
        text = fetch_mirbase_gff_text()
    except Exception:
        return 0

    n = 0
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        chrom, _src, ftype, start, end, _score, strand, _phase, attrs = parts[:9]
        if ftype not in ("miRNA", "miRNA_primary_transcript"):
            continue
        a = _parse_gff_attributes(attrs)
        name = a.get("Name") or a.get("ID") or f"mir_{start}_{end}"
        eid = f"mirbase_{name}"
        cur.execute(
            """INSERT OR IGNORE INTO regulatory_elements (
                 element_id, feature_type, chromosome, start, end, strand,
                 linked_gene_id, motif_name, metadata, source_id
               ) VALUES (?,?,?,?,?,?,?,?,?,?)""",
            (
                eid,
                "mirna",
                chrom,
                int(start),
                int(end),
                1 if strand == "+" else -1,
                None,
                name,
                dumps_metadata({"gff_type": ftype, "attrs": a}),
                sid,
            ),
        )
        n += 1
    conn.commit()
    return n


def seed_mirna_target_edges(
    conn: sqlite3.Connection,
    *,
    max_edges: int = 300,
) -> int:
    """Heuristic miRNA -> mRNA repression edges using symbol overlap (demo when targets unavailable)."""
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        ("MicroSalmon", "heuristic", "https://www.mdpi.com", now_iso()),
    )
    conn.commit()
    sid = cur.execute(
        "SELECT id FROM sources WHERE name=? ORDER BY id DESC LIMIT 1",
        ("MicroSalmon",),
    ).fetchone()[0]

    mirnas = cur.execute(
        """SELECT id, element_id, motif_name FROM regulatory_elements WHERE feature_type='mirna' LIMIT 200""",
    ).fetchall()
    genes = cur.execute(
        """SELECT id, symbol FROM genes WHERE biotype='protein_coding' AND symbol IS NOT NULL LIMIT 5000""",
    ).fetchall()
    if not mirnas or not genes:
        return 0

    gene_by_symbol = {g[1].upper(): g[0] for g in genes if g[1]}
    added = 0
    for el_id, _eid, motif in mirnas:
        if not motif:
            continue
        m = re.match(r"(ssa-let|ssa-mir)-(\d+)", str(motif), re.I)
        if not m:
            continue
        # Placeholder: link to arbitrary immune-related genes if present
        for sym in ("TLR3", "IRF3", "STAT1", "NFKB1"):
            if sym in gene_by_symbol and added < max_edges:
                tid = gene_by_symbol[sym]
                cur.execute(
                    """INSERT OR REPLACE INTO interactions (
                         regulator_gene_id, regulator_element_id, target_gene_id,
                         interaction_type, evidence_score, direction, metadata, source_id
                       ) VALUES (?,?,?,?,?,?,?,?)""",
                    (
                        None,
                        el_id,
                        tid,
                        "mirna_repression",
                        0.25,
                        "represses",
                        dumps_metadata({"note": "heuristic demo edge; replace with MicroSalmon targets"}),
                        sid,
                    ),
                )
                added += 1
    conn.commit()
    return added


def run(conn: sqlite3.Connection) -> None:
    load_mirna_features(conn)
    seed_mirna_target_edges(conn)
