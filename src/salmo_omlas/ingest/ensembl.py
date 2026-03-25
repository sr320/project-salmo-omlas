"""Ensembl BioMart + REST: genes, GO biological process, orthologs."""

from __future__ import annotations

import io
import sqlite3
import xml.etree.ElementTree as ET
import pandas as pd
import requests

from salmo_omlas.config import BIOMART_URL
from salmo_omlas.ingest._util import now_iso, save_raw


SALMON_TAXON = 8030

# Ensembl BioMart dataset name for Salmo salar (check Ensembl if this changes)
DATASET = "ssalar_gene_ensembl"


def _biomart_xml(attributes: list[str], filters: list[tuple[str, str]] | None = None) -> str:
    root = ET.Element("Query")
    root.set("virtualSchemaName", "default")
    root.set("formatter", "TSV")
    root.set("header", "1")
    root.set("uniqueRows", "0")
    root.set("count", "")
    root.set("datasetConfigVersion", "0.6")
    ds = ET.SubElement(root, "Dataset")
    ds.set("name", DATASET)
    ds.set("interface", "default")
    for name, value in filters or []:
        flt = ET.SubElement(ds, "Filter")
        flt.set("name", name)
        flt.set("value", value)
    for attr in attributes:
        a = ET.SubElement(ds, "Attribute")
        a.set("name", attr)
    # BioMart expects full XML document
    return '<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE Query>\n' + ET.tostring(
        root, encoding="unicode"
    )


def biomart_tsv(query_xml: str, attributes: list[str] | None = None) -> pd.DataFrame:
    """POST BioMart query; return DataFrame.

    If *attributes* is given, override the BioMart display-name headers with
    the internal attribute names (column order matches the query order).
    """
    r = requests.post(
        BIOMART_URL,
        data={"query": query_xml},
        timeout=600,
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    r.raise_for_status()
    raw = r.text
    if not raw.strip() or raw.strip().lower().startswith("error"):
        raise RuntimeError(f"BioMart error: {raw[:500]}")
    save_raw(f"biomart_{now_iso().replace(':', '-')}.tsv", raw)
    df = pd.read_csv(io.StringIO(raw), sep="\t", dtype=str, na_filter=False)
    if attributes and len(attributes) == len(df.columns):
        df.columns = attributes
    return df


def fetch_genes(limit: int | None = None) -> pd.DataFrame:
    attrs = [
        "ensembl_gene_id",
        "external_gene_name",
        "gene_biotype",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "description",
    ]
    xml = _biomart_xml(attrs)
    df = biomart_tsv(xml, attributes=attrs)
    df.columns = [str(c).strip() for c in df.columns]
    df = df.rename(
        columns={
            "external_gene_name": "symbol",
            "gene_biotype": "biotype",
            "chromosome_name": "chromosome",
            "start_position": "start",
            "end_position": "end",
        }
    )
    for c in ("start", "end", "strand"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    if limit:
        df = df.head(limit)
    return df


def fetch_gene_go() -> pd.DataFrame:
    attrs = [
        "ensembl_gene_id",
        "external_gene_name",
        "go_id",
        "name_1006",
        "namespace_1003",
    ]
    xml = _biomart_xml(attrs)
    df = biomart_tsv(xml, attributes=attrs)
    df.columns = [c.strip() for c in df.columns]
    rename = {
        "external_gene_name": "symbol",
        "name_1006": "go_name",
        "namespace_1003": "go_namespace",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    df = df[df["go_id"].str.startswith("GO:", na=False)]
    df = df[df["go_namespace"].str.contains("biological_process", case=False, na=False)]
    return df


def fetch_orthologs() -> pd.DataFrame:
    """One-to-many: multiple homolog columns per species."""
    attrs = [
        "ensembl_gene_id",
        "external_gene_name",
        "drerio_homolog_ensembl_gene",
        "drerio_homolog_associated_gene_name",
        "mmusculus_homolog_ensembl_gene",
        "mmusculus_homolog_associated_gene_name",
        "hsapiens_homolog_ensembl_gene",
        "hsapiens_homolog_associated_gene_name",
    ]
    xml = _biomart_xml(attrs)
    df = biomart_tsv(xml, attributes=attrs)
    df.columns = [c.strip() for c in df.columns]
    return df


def load_into_sqlite(
    conn: sqlite3.Connection,
    *,
    genes_df: pd.DataFrame | None = None,
    go_df: pd.DataFrame | None = None,
    ortho_df: pd.DataFrame | None = None,
    source_name: str = "Ensembl",
    source_version: str = "BioMart",
) -> None:
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        (source_name, source_version, "https://www.ensembl.org", now_iso()),
    )
    conn.commit()
    sid = cur.execute(
        "SELECT id FROM sources WHERE name=? AND version=? ORDER BY id DESC LIMIT 1",
        (source_name, source_version),
    ).fetchone()[0]

    cur.execute(
        """INSERT OR IGNORE INTO species (taxon_id, scientific_name, common_name, assembly)
           VALUES (?, ?, ?, ?)""",
        (SALMON_TAXON, "Salmo salar", "Atlantic salmon", "Ssal_v3.1"),
    )
    conn.commit()
    sp_row = cur.execute("SELECT id FROM species WHERE taxon_id=?", (SALMON_TAXON,)).fetchone()
    species_id = sp_row[0]

    if genes_df is not None and not genes_df.empty:
        for _, row in genes_df.iterrows():
            gid = row.get("ensembl_gene_id")
            if not gid:
                continue
            strand = int(row.get("strand", 0) or 0)
            cur.execute(
                """INSERT OR REPLACE INTO genes (
                     ensembl_gene_id, symbol, name, biotype, chromosome, start, end, strand,
                     description, species_id, source_id
                   ) VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                (
                    gid,
                    row.get("symbol") or None,
                    row.get("symbol") or None,
                    row.get("biotype") or None,
                    str(row.get("chromosome") or ""),
                    int(row.get("start") or 0),
                    int(row.get("end") or 0),
                    strand,
                    row.get("description") or None,
                    species_id,
                    sid,
                ),
            )
        conn.commit()

    gene_pk = dict(cur.execute("SELECT ensembl_gene_id, id FROM genes").fetchall())

    if go_df is not None and not go_df.empty:
        for _, row in go_df.iterrows():
            gid = row.get("ensembl_gene_id")
            go_id = row.get("go_id")
            if not gid or not go_id or gid not in gene_pk:
                continue
            name = row.get("go_name") or go_id
            cur.execute(
                """INSERT OR IGNORE INTO physiological_processes (term_id, name, namespace, parent_term_id, definition)
                   VALUES (?,?,?,?,?)""",
                (go_id, name, "go_bp", None, None),
            )
            conn.commit()
            pid = cur.execute(
                "SELECT id FROM physiological_processes WHERE term_id=?",
                (go_id,),
            ).fetchone()[0]
            cur.execute(
                """INSERT OR REPLACE INTO gene_processes (gene_id, process_id, evidence) VALUES (?,?,?)""",
                (gene_pk[gid], pid, "IEA"),
            )
        conn.commit()

    if ortho_df is not None and not ortho_df.empty:
        colmap = {
            "drerio_homolog_ensembl_gene": 7955,
            "mmusculus_homolog_ensembl_gene": 10090,
            "hsapiens_homolog_ensembl_gene": 9606,
        }
        symcols = {
            "drerio_homolog_ensembl_gene": "drerio_homolog_associated_gene_name",
            "mmusculus_homolog_ensembl_gene": "mmusculus_homolog_associated_gene_name",
            "hsapiens_homolog_ensembl_gene": "hsapiens_homolog_associated_gene_name",
        }
        for _, row in ortho_df.iterrows():
            gid = row.get("ensembl_gene_id")
            if not gid or gid not in gene_pk:
                continue
            for hcol, tax in colmap.items():
                oid = row.get(hcol)
                if not oid or oid == "":
                    continue
                sym = row.get(symcols.get(hcol, ""), "") or None
                cur.execute(
                    """INSERT OR REPLACE INTO orthologs (
                         gene_id, ortholog_ensembl_id, ortholog_symbol, ortholog_name,
                         species_taxon_id, ortholog_type, identity_percent, confidence, source_id
                       ) VALUES (?,?,?,?,?,?,?,?,?)""",
                    (
                        gene_pk[gid],
                        oid,
                        sym,
                        sym,
                        tax,
                        "ortholog_one2one",
                        None,
                        0.9,
                        sid,
                    ),
                )
        conn.commit()


def run(
    conn: sqlite3.Connection,
    *,
    limit_genes: int | None = None,
    skip_genes: bool = False,
    skip_go: bool = False,
    skip_orthologs: bool = False,
) -> None:
    """Fetch from BioMart and load."""
    genes_df = None
    go_df = None
    ortho_df = None
    if not skip_genes:
        genes_df = fetch_genes(limit=limit_genes)
    if not skip_go:
        go_df = fetch_gene_go()
        if limit_genes and genes_df is not None:
            allowed = set(genes_df["ensembl_gene_id"])
            go_df = go_df[go_df["ensembl_gene_id"].isin(allowed)]
    if not skip_orthologs:
        ortho_df = fetch_orthologs()
        if limit_genes and genes_df is not None:
            allowed = set(genes_df["ensembl_gene_id"])
            ortho_df = ortho_df[ortho_df["ensembl_gene_id"].isin(allowed)]

    load_into_sqlite(conn, genes_df=genes_df, go_df=go_df, ortho_df=ortho_df)
