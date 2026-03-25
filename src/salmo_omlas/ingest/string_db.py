"""STRING protein association network for Salmo salar."""

from __future__ import annotations

import io
import sqlite3

import pandas as pd
import requests

from salmo_omlas.config import STRING_API, STRING_SPECIES_NCBI_TAXON
from salmo_omlas.ingest._util import dumps_metadata, now_iso


def fetch_network_tsv(
    ensembl_ids: list[str],
    *,
    required_score: int = 400,
) -> pd.DataFrame:
    """Bulk network via STRING tsv/network endpoint."""
    if not ensembl_ids:
        return pd.DataFrame()
    joined = "%0d".join(ensembl_ids)
    url = f"{STRING_API}/tsv/network"
    params = {
        "identifiers": joined,
        "species": STRING_SPECIES_NCBI_TAXON,
        "required_score": required_score,
    }
    r = requests.get(url, params=params, timeout=180)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep="\t")


def _strip_taxon(string_id: str) -> str:
    s = str(string_id)
    if "." in s and "ENS" in s.upper():
        return s.split(".", 1)[-1]
    return s


def load_interactions_for_genes(
    conn: sqlite3.Connection,
    *,
    limit_genes: int | None = 500,
    required_score: int = 400,
) -> None:
    """Add PPI-style edges between salmon genes using STRING."""
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        ("STRING", "v12", "https://string-db.org", now_iso()),
    )
    conn.commit()
    sid = cur.execute(
        "SELECT id FROM sources WHERE name=? ORDER BY id DESC LIMIT 1",
        ("STRING",),
    ).fetchone()[0]

    genes = cur.execute(
        "SELECT ensembl_gene_id FROM genes WHERE biotype='protein_coding' LIMIT ?",
        (limit_genes or 10_000_000,),
    ).fetchall()
    ids = [g[0] for g in genes]
    if not ids:
        return

    gene_pk = dict(cur.execute("SELECT ensembl_gene_id, id FROM genes").fetchall())

    chunk = 50
    for i in range(0, len(ids), chunk):
        batch = ids[i : i + chunk]
        try:
            df = fetch_network_tsv(batch, required_score=required_score)
        except Exception:
            continue
        if df.empty:
            continue
        # Typical columns: node1, node2, combined_score (or stringId_A / stringId_B)
        colmap = {c.lower(): c for c in df.columns}
        c1 = colmap.get("node1") or colmap.get("stringid_a") or colmap.get("protein1")
        c2 = colmap.get("node2") or colmap.get("stringid_b") or colmap.get("protein2")
        sc = colmap.get("combined_score") or colmap.get("score")
        if not c1 or not c2:
            c1, c2 = df.columns[0], df.columns[1]
        if sc is None and len(df.columns) > 2:
            sc = df.columns[2]
        for _, row in df.iterrows():
            id_a = _strip_taxon(row[c1])
            id_b = _strip_taxon(row[c2])
            try:
                score = float(row[sc]) / 1000.0 if sc is not None else None
            except (TypeError, ValueError):
                score = None
            if id_a not in gene_pk or id_b not in gene_pk:
                continue
            cur.execute(
                """INSERT OR REPLACE INTO interactions (
                     regulator_gene_id, regulator_element_id, target_gene_id,
                     interaction_type, evidence_score, direction, metadata, source_id
                   ) VALUES (?,?,?,?,?,?,?,?)""",
                (
                    gene_pk[id_a],
                    None,
                    gene_pk[id_b],
                    "ppi",
                    score,
                    "associated",
                    dumps_metadata({"source": "STRING"}),
                    sid,
                ),
            )
    conn.commit()
