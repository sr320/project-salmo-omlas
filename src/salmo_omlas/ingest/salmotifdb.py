"""SalMotifDB / motif-derived regulatory features (best-effort ingestion)."""

from __future__ import annotations

import csv
import io
import sqlite3
from pathlib import Path

import requests

from salmo_omlas.config import DATA_RAW, ensure_dirs
from salmo_omlas.ingest._util import dumps_metadata, now_iso

# Public reference for motif metadata (users can place full export under data/raw/salmotif_export.tsv)
SALMOTIF_PLACEHOLDER_URL = None


def load_from_tsv(conn: sqlite3.Connection, path: Path) -> int:
    """Load rows: motif_name, chromosome, start, end, strand, linked_ensembl_gene_id (optional)."""
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        ("SalMotifDB", "custom", "https://salmobase.org/apps/SalMotifDB", now_iso()),
    )
    conn.commit()
    sid = cur.execute(
        "SELECT id FROM sources WHERE name=? ORDER BY id DESC LIMIT 1",
        ("SalMotifDB",),
    ).fetchone()[0]

    gene_pk = dict(cur.execute("SELECT ensembl_gene_id, id FROM genes").fetchall())
    n = 0
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row.get("linked_ensembl_gene_id") or row.get("gene_id")
            if not gid or gid not in gene_pk:
                continue
            eid = row.get("element_id") or f"tfbs_{gid}_{row.get('start')}_{row.get('end')}"
            cur.execute(
                """INSERT OR IGNORE INTO regulatory_elements (
                     element_id, feature_type, chromosome, start, end, strand,
                     linked_gene_id, motif_name, metadata, source_id
                   ) VALUES (?,?,?,?,?,?,?,?,?,?)""",
                (
                    eid,
                    "tfbs",
                    row.get("chromosome") or row.get("chr"),
                    int(row.get("start") or 0),
                    int(row.get("end") or 0),
                    int(row.get("strand") or 0),
                    gene_pk[gid],
                    row.get("motif_name") or row.get("motif"),
                    dumps_metadata({"raw": row}),
                    sid,
                ),
            )
            el_row = cur.execute(
                "SELECT id FROM regulatory_elements WHERE element_id=?",
                (eid,),
            ).fetchone()
            if not el_row:
                continue
            el_id = el_row[0]
            # TF -> target gene edge (motif near gene promoter)
            cur.execute(
                """INSERT OR REPLACE INTO interactions (
                     regulator_gene_id, regulator_element_id, target_gene_id,
                     interaction_type, evidence_score, direction, metadata, source_id
                   ) VALUES (?,?,?,?,?,?,?,?)""",
                (
                    None,
                    el_id,
                    gene_pk[gid],
                    "tf_to_gene",
                    float(row.get("score") or 0.5),
                    "binds",
                    dumps_metadata({"motif": row.get("motif_name")}),
                    sid,
                ),
            )
            n += 1
    conn.commit()
    return n


def seed_placeholder_motifs(conn: sqlite3.Connection, *, max_rows: int = 200) -> int:
    """Create synthetic TFBS-linked rows for genes that lack element data (demo / tests)."""
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        ("SalMotifDB", "placeholder", "https://salmobase.org/apps/SalMotifDB", now_iso()),
    )
    conn.commit()
    sid = cur.execute(
        "SELECT id FROM sources WHERE name=? ORDER BY id DESC LIMIT 1",
        ("SalMotifDB",),
    ).fetchone()[0]

    rows = cur.execute(
        """SELECT id, ensembl_gene_id, chromosome, start, end, strand FROM genes
           WHERE biotype='protein_coding' LIMIT ?""",
        (max_rows,),
    ).fetchall()
    n = 0
    for gid_pk, ens, chrom, gstart, gend, strand in rows:
        # Synthetic promoter-proximal TFBS 500bp upstream of TSS
        if strand >= 0:
            s, e = max(0, gstart - 600), max(0, gstart - 100)
        else:
            s, e = gend + 100, gend + 600
        eid = f"demo_tfbs_{ens}_{s}_{e}"
        cur.execute(
            """INSERT OR IGNORE INTO regulatory_elements (
                 element_id, feature_type, chromosome, start, end, strand,
                 linked_gene_id, motif_name, metadata, source_id
               ) VALUES (?,?,?,?,?,?,?,?,?,?)""",
            (
                eid,
                "tfbs",
                str(chrom),
                int(s),
                int(e),
                int(strand or 1),
                gid_pk,
                "NFKB-like_demo",
                dumps_metadata({"note": "placeholder TFBS for development"}),
                sid,
            ),
        )
        el = cur.execute(
            "SELECT id FROM regulatory_elements WHERE element_id=?",
            (eid,),
        ).fetchone()
        if not el:
            continue
        cur.execute(
            """INSERT OR REPLACE INTO interactions (
                 regulator_gene_id, regulator_element_id, target_gene_id,
                 interaction_type, evidence_score, direction, metadata, source_id
               ) VALUES (?,?,?,?,?,?,?,?)""",
            (
                None,
                el[0],
                gid_pk,
                "tf_to_gene",
                0.35,
                "binds",
                dumps_metadata({"placeholder": True}),
                sid,
            ),
        )
        n += 1
    conn.commit()
    return n


def run(conn: sqlite3.Connection, *, tsv_path: Path | None = None) -> None:
    ensure_dirs()
    path = tsv_path or (DATA_RAW / "salmotif_export.tsv")
    if path.exists():
        load_from_tsv(conn, path)
    else:
        seed_placeholder_motifs(conn)
