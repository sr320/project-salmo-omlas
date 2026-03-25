"""KEGG pathways for Salmo salar (sasa)."""

from __future__ import annotations

import re
import sqlite3

import pandas as pd
import requests

from salmo_omlas.config import IMMUNE_KEGG_PATHWAY_IDS, KEGG_BASE
from salmo_omlas.ingest._util import http_get, now_iso


def list_pathways(org: str = "sasa") -> pd.DataFrame:
    """path:map<tab>name lines."""
    url = f"{KEGG_BASE}/list/pathway/{org}"
    r = http_get(url)
    lines = [ln for ln in r.text.strip().splitlines() if ln.strip()]
    rows = []
    for ln in lines:
        parts = ln.split("\t", 1)
        if len(parts) == 2:
            pid, name = parts[0], parts[1]
            m = re.search(r"path:(\w+)", pid)
            term = m.group(1) if m else pid.replace("path:", "")
            rows.append({"term_id": term, "name": name.strip()})
    return pd.DataFrame(rows)


def link_genes_pathway(pathway_term: str, org: str = "sasa") -> pd.DataFrame:
    """Genes in pathway via KEGG link (path:sasaXXXXX)."""
    path_id = pathway_term if pathway_term.startswith("path:") else f"path:{pathway_term}"
    url = f"{KEGG_BASE}/link/{org}/{path_id}"
    r = http_get(url)
    rows = []
    for ln in r.text.strip().splitlines():
        if not ln.strip():
            continue
        parts = ln.split("\t", 1)
        if len(parts) != 2:
            continue
        a, b = parts
        rows.append({"kegg_gene": a.strip(), "pathway": b.strip()})
    return pd.DataFrame(rows)


def _ensembl_id_from_kegg_conv_line(line: str) -> str | None:
    """Parse Ensembl gene id from KEGG conv line."""
    parts = line.split("\t", 1)
    if len(parts) != 2:
        return None
    right = parts[1]
    m = re.search(r"(ENS[A-Z0-9]+)", right)
    return m.group(1) if m else None


def load_kegg_conv_map() -> dict[str, str]:
    """Map KEGG gene id (sasa:nnn) -> Ensembl gene id."""
    r = requests.get(f"{KEGG_BASE}/conv/sasa/ensembl", timeout=180)
    r.raise_for_status()
    m: dict[str, str] = {}
    for ln in r.text.strip().splitlines():
        if not ln.strip():
            continue
        parts = ln.split("\t", 1)
        if len(parts) != 2:
            continue
        kg = parts[0].strip()
        ens = _ensembl_id_from_kegg_conv_line(ln)
        if ens:
            m[kg] = ens
    return m


def load_pathways_and_members(
    conn: sqlite3.Connection,
    *,
    immune_only: bool = True,
) -> None:
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        ("KEGG", "REST", KEGG_BASE, now_iso()),
    )
    conn.commit()

    pathways = list_pathways("sasa")
    if immune_only and not pathways.empty:
        pathways = pathways[
            pathways["term_id"].apply(
                lambda t: any(str(t).endswith(suf) for suf in IMMUNE_KEGG_PATHWAY_IDS)
            )
        ]

    for _, prow in pathways.iterrows():
        tid = str(prow["term_id"])
        name = str(prow["name"])
        cur.execute(
            """INSERT OR IGNORE INTO physiological_processes (term_id, name, namespace, parent_term_id, definition)
               VALUES (?,?,?,?,?)""",
            (tid, name, "kegg_pathway", None, None),
        )
    conn.commit()

    gene_pk = dict(cur.execute("SELECT ensembl_gene_id, id FROM genes").fetchall())
    try:
        kegg_to_ens = load_kegg_conv_map()
    except Exception:
        kegg_to_ens = {}

    for _, prow in pathways.iterrows():
        tid = str(prow["term_id"])
        pid_row = cur.execute(
            "SELECT id FROM physiological_processes WHERE term_id=?",
            (tid,),
        ).fetchone()
        if not pid_row:
            continue
        process_id = pid_row[0]
        try:
            df = link_genes_pathway(tid, "sasa")
        except Exception:
            continue
        for _, row in df.iterrows():
            kg = str(row.get("kegg_gene", "")).strip()
            ens = kegg_to_ens.get(kg)
            if ens and ens in gene_pk:
                cur.execute(
                    """INSERT OR REPLACE INTO gene_processes (gene_id, process_id, evidence) VALUES (?,?,?)""",
                    (gene_pk[ens], process_id, "KEGG"),
                )
    conn.commit()
