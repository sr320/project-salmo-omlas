"""Core relational queries and light network traversal."""

from __future__ import annotations

import sqlite3
from typing import Any

import networkx as nx
import pandas as pd


def _gene_where_clause(identifier: str) -> tuple[str, tuple[Any, ...]]:
    ident = identifier.strip()
    if ident.startswith("ENS"):
        return "g.ensembl_gene_id = ?", (ident,)
    return "(g.symbol = ? COLLATE NOCASE OR g.name = ? COLLATE NOCASE)", (ident, ident)


def query_gene(conn: sqlite3.Connection, identifier: str) -> dict[str, pd.DataFrame]:
    """Return tables for gene metadata, processes, orthologs, elements, interactions, CRISPR."""
    clause, params = _gene_where_clause(identifier)
    gene_sql = f"""SELECT g.*, sp.scientific_name, sp.assembly
                   FROM genes g JOIN species sp ON sp.id = g.species_id
                   WHERE {clause} LIMIT 5"""
    gene_df = pd.read_sql_query(gene_sql, conn, params=params)
    if gene_df.empty:
        return {
            "gene": gene_df,
            "processes": pd.DataFrame(),
            "orthologs": pd.DataFrame(),
            "regulatory_elements": pd.DataFrame(),
            "interactions_out": pd.DataFrame(),
            "interactions_in": pd.DataFrame(),
            "crispr_targets": pd.DataFrame(),
        }

    gid = int(gene_df.iloc[0]["id"])
    ens = str(gene_df.iloc[0]["ensembl_gene_id"])

    processes = pd.read_sql_query(
        """SELECT pp.term_id, pp.name, pp.namespace, gp.evidence
           FROM gene_processes gp
           JOIN physiological_processes pp ON pp.id = gp.process_id
           WHERE gp.gene_id = ?""",
        conn,
        params=(gid,),
    )
    orthologs = pd.read_sql_query(
        """SELECT o.* FROM orthologs o WHERE o.gene_id = ?""",
        conn,
        params=(gid,),
    )
    elements = pd.read_sql_query(
        """SELECT DISTINCT re.*
           FROM regulatory_elements re
           WHERE re.linked_gene_id = ? OR re.motif_name LIKE ?
           UNION
           SELECT DISTINCT re2.*
           FROM interactions i
           JOIN regulatory_elements re2 ON re2.id = i.regulator_element_id
           WHERE i.target_gene_id = ?""",
        conn,
        params=(gid, f"%{identifier}%", gid),
    )
    inter_out = pd.read_sql_query(
        """SELECT i.*, tg.symbol AS target_symbol, tg.ensembl_gene_id AS target_ensembl
           FROM interactions i
           JOIN genes tg ON tg.id = i.target_gene_id
           WHERE i.regulator_gene_id = ?""",
        conn,
        params=(gid,),
    )
    inter_in = pd.read_sql_query(
        """SELECT i.*, rg.symbol AS regulator_symbol, rg.ensembl_gene_id AS regulator_ensembl,
                  re.element_id AS regulator_element
           FROM interactions i
           LEFT JOIN genes rg ON rg.id = i.regulator_gene_id
           LEFT JOIN regulatory_elements re ON re.id = i.regulator_element_id
           WHERE i.target_gene_id = ?""",
        conn,
        params=(gid,),
    )
    crispr = pd.read_sql_query(
        """SELECT * FROM crispr_targets WHERE gene_id = ? ORDER BY on_target_score DESC""",
        conn,
        params=(gid,),
    )
    return {
        "gene": gene_df,
        "processes": processes,
        "orthologs": orthologs,
        "regulatory_elements": elements,
        "interactions_out": inter_out,
        "interactions_in": inter_in,
        "crispr_targets": crispr,
        "_internal_gene_id": gid,
        "_ensembl_gene_id": ens,
    }


def query_process(conn: sqlite3.Connection, text: str) -> dict[str, pd.DataFrame]:
    """Match physiological processes by name or term id; return member genes."""
    q = f"%{text.strip()}%"
    proc_df = pd.read_sql_query(
        """SELECT * FROM physiological_processes
           WHERE term_id LIKE ? OR name LIKE ?""",
        conn,
        params=(q, q),
    )
    if proc_df.empty:
        return {"processes": proc_df, "genes": pd.DataFrame()}
    ids = tuple(proc_df["id"].tolist())
    placeholders = ",".join("?" * len(ids))
    genes = pd.read_sql_query(
        f"""SELECT g.*, pp.term_id, pp.name AS process_name, gp.evidence
            FROM gene_processes gp
            JOIN genes g ON g.id = gp.gene_id
            JOIN physiological_processes pp ON pp.id = gp.process_id
            WHERE gp.process_id IN ({placeholders})""",
        conn,
        params=ids,
    )
    return {"processes": proc_df, "genes": genes}


def query_element(conn: sqlite3.Connection, element_id: str) -> dict[str, pd.DataFrame]:
    """Lookup regulatory element by element_id."""
    el = pd.read_sql_query(
        """SELECT re.*, g.symbol, g.ensembl_gene_id
           FROM regulatory_elements re
           LEFT JOIN genes g ON g.id = re.linked_gene_id
           WHERE re.element_id = ? OR re.motif_name = ?""",
        conn,
        params=(element_id, element_id),
    )
    if el.empty:
        return {"element": el, "interactions": pd.DataFrame()}
    eid = int(el.iloc[0]["id"])
    inter = pd.read_sql_query(
        """SELECT i.*, g.symbol AS target_symbol FROM interactions i
           JOIN genes g ON g.id = i.target_gene_id
           WHERE i.regulator_element_id = ?""",
        conn,
        params=(eid,),
    )
    return {"element": el, "interactions": inter}


def _gene_label(conn: sqlite3.Connection, gene_id: int) -> str:
    sym = conn.execute(
        "SELECT symbol, ensembl_gene_id FROM genes WHERE id=?",
        (gene_id,),
    ).fetchone()
    if not sym:
        return str(gene_id)
    return str(sym[0] or sym[1])


def regulatory_subgraph(
    conn: sqlite3.Connection,
    identifier: str,
    *,
    depth: int = 1,
) -> nx.DiGraph:
    """BFS over interactions from a seed gene (by symbol or Ensembl id)."""
    G = nx.DiGraph()
    res = query_gene(conn, identifier)
    if res["gene"].empty:
        return G
    seed = int(res["gene"].iloc[0]["id"])
    seed_label = str(res["gene"].iloc[0].get("symbol") or res["gene"].iloc[0]["ensembl_gene_id"])
    G.add_node(seed, label=seed_label, kind="gene")

    frontier: set[int] = {seed}
    gene_seen: set[int] = {seed}
    for _ in range(max(1, depth)):
        next_genes: set[int] = set()
        for gid in frontier:
            rows = conn.execute(
                """SELECT i.target_gene_id, i.regulator_gene_id, i.regulator_element_id,
                          i.interaction_type, i.evidence_score
                   FROM interactions i
                   WHERE i.regulator_gene_id = ? OR i.target_gene_id = ?""",
                (gid, gid),
            ).fetchall()
            for tgt_id, reg_g, reg_el, itype, escore in rows:
                if reg_g is not None:
                    if reg_g not in gene_seen:
                        G.add_node(reg_g, label=_gene_label(conn, reg_g), kind="gene")
                        gene_seen.add(reg_g)
                        next_genes.add(reg_g)
                    if tgt_id not in gene_seen:
                        G.add_node(tgt_id, label=_gene_label(conn, tgt_id), kind="gene")
                        gene_seen.add(tgt_id)
                        next_genes.add(tgt_id)
                    G.add_edge(reg_g, tgt_id, type=itype, score=escore)
                elif reg_el is not None:
                    elab = conn.execute(
                        "SELECT element_id, feature_type FROM regulatory_elements WHERE id=?",
                        (reg_el,),
                    ).fetchone()
                    node_id = f"el_{reg_el}"
                    if node_id not in G:
                        G.add_node(
                            node_id,
                            label=elab[0] if elab else str(reg_el),
                            kind=str(elab[1]) if elab else "element",
                        )
                    if tgt_id not in gene_seen:
                        G.add_node(tgt_id, label=_gene_label(conn, tgt_id), kind="gene")
                        gene_seen.add(tgt_id)
                        next_genes.add(tgt_id)
                    G.add_edge(node_id, tgt_id, type=itype, score=escore)
        frontier = next_genes
    return G
