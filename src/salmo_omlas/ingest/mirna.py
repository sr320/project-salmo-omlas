"""miRNA genes and target edges (miRBase GFF + heuristic targets)."""

from __future__ import annotations

import io
import re
import sqlite3

import requests

from salmo_omlas.config import MIRBASE_MATURE_FASTA
from salmo_omlas.ingest._util import dumps_metadata, now_iso


def fetch_mirbase_mature_fasta() -> str:
    r = requests.get(MIRBASE_MATURE_FASTA, timeout=180)
    r.raise_for_status()
    return r.text


def load_mirna_features(conn: sqlite3.Connection) -> int:
    """Load miRBase mature miRNAs (FASTA) for Salmo salar (ssa-).

    miRBase does not consistently provide genome-coordinate GFF3 for salmon,
    so these are stored as regulatory elements without genomic coordinates.
    """
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        ("miRBase", "CURRENT", MIRBASE_MATURE_FASTA, now_iso()),
    )
    conn.commit()
    sid = cur.execute(
        "SELECT id FROM sources WHERE name=? ORDER BY id DESC LIMIT 1",
        ("miRBase",),
    ).fetchone()[0]

    try:
        fasta = fetch_mirbase_mature_fasta()
    except Exception:
        return 0

    n = 0
    name: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        nonlocal n, name, seq_parts
        if not name:
            return
        if not name.lower().startswith("ssa-"):
            name = None
            seq_parts = []
            return
        seq = "".join(seq_parts).strip().upper()
        if not seq:
            name = None
            seq_parts = []
            return
        eid = f"mirbase_{name}"
        cur.execute(
            """INSERT OR IGNORE INTO regulatory_elements (
                 element_id, feature_type, chromosome, start, end, strand,
                 linked_gene_id, motif_name, metadata, source_id
               ) VALUES (?,?,?,?,?,?,?,?,?,?)""",
            (
                eid,
                "mirna",
                None,
                None,
                None,
                None,
                None,
                name,
                dumps_metadata({"sequence": seq, "length": len(seq), "kind": "mature_fasta"}),
                sid,
            ),
        )
        n += 1
        name = None
        seq_parts = []

    for line in io.StringIO(fasta):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            # Take first token after '>' as miRNA id (e.g., ssa-miR-...)
            name = line[1:].split(None, 1)[0]
            continue
        seq_parts.append(line)

    flush()
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
