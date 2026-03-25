"""CRISPR-oriented queries."""

from __future__ import annotations

import sqlite3

import pandas as pd


def list_crispr_for_gene(conn: sqlite3.Connection, identifier: str) -> pd.DataFrame:
    ident = identifier.strip()
    if ident.upper().startswith("ENS"):
        where = "g.ensembl_gene_id = ?"
        params: tuple[str, ...] = (ident,)
    else:
        where = "(g.symbol = ? COLLATE NOCASE OR g.name = ? COLLATE NOCASE)"
        params = (ident, ident)
    return pd.read_sql_query(
        f"""SELECT ct.*, g.symbol, g.ensembl_gene_id
            FROM crispr_targets ct
            JOIN genes g ON g.id = ct.gene_id
            WHERE {where}""",
        conn,
        params=params,
    )


def rank_genes_for_pathway_crispr(conn: sqlite3.Connection, pathway_query: str) -> pd.DataFrame:
    """Join process search with best available CRISPR guide scores."""
    q = f"%{pathway_query.strip()}%"
    return pd.read_sql_query(
        """SELECT g.symbol, g.ensembl_gene_id, pp.term_id, pp.name AS pathway,
                  MAX(ct.on_target_score) AS best_on_target,
                  MIN(ct.off_target_risk) AS best_off_target_risk,
                  COUNT(ct.id) AS n_guides
           FROM physiological_processes pp
           JOIN gene_processes gp ON gp.process_id = pp.id
           JOIN genes g ON g.id = gp.gene_id
           LEFT JOIN crispr_targets ct ON ct.gene_id = g.id
           WHERE pp.term_id LIKE ? OR pp.name LIKE ?
           GROUP BY g.id, pp.id
           ORDER BY best_on_target DESC, best_off_target_risk ASC""",
        conn,
        params=(q, q),
    )
