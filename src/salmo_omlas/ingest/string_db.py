"""STRING protein association network for Salmo salar."""

from __future__ import annotations

import io
import sqlite3
import time

import pandas as pd
import requests

from salmo_omlas.config import STRING_API, STRING_SPECIES_NCBI_TAXON
from salmo_omlas.ingest._util import dumps_metadata, now_iso


def fetch_string_ids_tsv(
    identifiers: list[str],
    *,
    species: int = STRING_SPECIES_NCBI_TAXON,
) -> pd.DataFrame:
    """Resolve gene/protein names to STRING identifiers."""
    if not identifiers:
        return pd.DataFrame()
    joined = "%0d".join(identifiers)
    url = f"{STRING_API}/tsv/get_string_ids"
    params = {"identifiers": joined, "species": species, "echo_query": 1}
    r = requests.get(url, params=params, timeout=180)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep="\t")


def fetch_network_tsv(
    string_ids: list[str],
    *,
    required_score: int = 400,
) -> pd.DataFrame:
    """Bulk network via STRING tsv/network endpoint."""
    if not string_ids:
        return pd.DataFrame()
    joined = "\r".join(string_ids)
    url = f"{STRING_API}/tsv/network"
    data = {
        "identifiers": joined,
        "species": STRING_SPECIES_NCBI_TAXON,
        "required_score": required_score,
    }
    r = requests.post(url, data=data, timeout=180)
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
    """Add PPI-style edges between salmon genes using STRING.

    Notes:
    - STRING's network endpoint expects identifiers it can resolve (often protein names/STRING ids).
    - Ensembl *gene* IDs generally won't resolve directly, so we map via gene symbols first.
    """
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

    rows = cur.execute(
        """
        SELECT id, symbol
        FROM genes
        WHERE biotype='protein_coding' AND symbol IS NOT NULL AND symbol != ''
        LIMIT ?
        """,
        (limit_genes or 10_000_000,),
    ).fetchall()
    if not rows:
        return

    gene_by_symbol = {str(sym).upper(): gid for gid, sym in rows if sym}

    symbols = [str(sym) for _gid, sym in rows if sym]

    # Resolve symbols -> STRING ids (batched). Be polite to the API.
    string_ids: list[str] = []
    chunk = 200
    for i in range(0, len(symbols), chunk):
        batch = symbols[i : i + chunk]
        try:
            mapped = fetch_string_ids_tsv(batch)
        except Exception:
            continue
        if not mapped.empty and "stringId" in mapped.columns:
            string_ids.extend([str(x) for x in mapped["stringId"].dropna().tolist() if str(x)])
        time.sleep(1)

    # Fetch interactions among the mapped ids.
    if not string_ids:
        return

    net_chunk = 50
    for i in range(0, len(string_ids), net_chunk):
        batch = string_ids[i : i + net_chunk]
        try:
            df = fetch_network_tsv(batch, required_score=required_score)
        except Exception:
            continue
        if df.empty:
            continue

        colmap = {c.lower(): c for c in df.columns}
        a_name = colmap.get("preferredname_a")
        b_name = colmap.get("preferredname_b")
        sc = colmap.get("score") or colmap.get("combined_score")
        if not a_name or not b_name:
            continue

        for _, row in df.iterrows():
            sa = str(row[a_name]).upper()
            sb = str(row[b_name]).upper()
            if sa not in gene_by_symbol or sb not in gene_by_symbol:
                continue
            try:
                score = float(row[sc]) if sc is not None else None
            except (TypeError, ValueError):
                score = None
            cur.execute(
                """INSERT OR REPLACE INTO interactions (
                     regulator_gene_id, regulator_element_id, target_gene_id,
                     interaction_type, evidence_score, direction, metadata, source_id
                   ) VALUES (?,?,?,?,?,?,?,?)""",
                (
                    gene_by_symbol[sa],
                    None,
                    gene_by_symbol[sb],
                    "ppi",
                    score,
                    "associated",
                    dumps_metadata({"source": "STRING"}),
                    sid,
                ),
            )
        conn.commit()
        time.sleep(1)
    conn.commit()
