"""Command-line interface."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click
import pandas as pd

from salmo_omlas import __version__
from salmo_omlas.config import get_db_path
from salmo_omlas.db import connect, init_schema, table_exists
from salmo_omlas.demo_seed import seed as demo_seed
from salmo_omlas.ingest import crispr as ingest_crispr
from salmo_omlas.ingest import ensembl, mirna, salmotifdb, string_db
from salmo_omlas.ingest import kegg as kegg_ingest
from salmo_omlas.query.core import query_gene, query_process, regulatory_subgraph
from salmo_omlas.query.crispr_queries import list_crispr_for_gene, rank_genes_for_pathway_crispr


@click.group()
@click.version_option(__version__, prog_name="salmo-omlas")
def main() -> None:
    """Salmon gene regulatory interaction database CLI."""


@main.command("init-db")
@click.option("--db", "db_path", type=click.Path(), default=None, help="SQLite path (default from config)")
@click.option("--force", is_flag=True, help="Drop all tables and recreate (destructive)")
def init_db(db_path: str | None, force: bool) -> None:
    """Create database file and apply schema."""
    path = Path(db_path) if db_path else get_db_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = connect(path)
    if force:
        conn.executescript(
            """
            PRAGMA foreign_keys = OFF;
            DROP TABLE IF EXISTS crispr_targets;
            DROP TABLE IF EXISTS interactions;
            DROP TABLE IF EXISTS regulatory_elements;
            DROP TABLE IF EXISTS orthologs;
            DROP TABLE IF EXISTS gene_processes;
            DROP TABLE IF EXISTS physiological_processes;
            DROP TABLE IF EXISTS genes;
            DROP TABLE IF EXISTS species;
            DROP TABLE IF EXISTS sources;
            PRAGMA foreign_keys = ON;
            """
        )
        conn.commit()
    init_schema(conn)
    click.echo(f"Initialized schema at {path}")


@main.command()
@click.option("--db", "db_path", type=click.Path(), default=None)
def demo(db_path: str | None) -> None:
    """Load minimal demo rows (offline)."""
    conn = connect(db_path)
    if not table_exists(conn, "genes"):
        init_schema(conn)
    demo_seed(conn)
    click.echo("Demo data loaded.")


@main.command()
@click.option("--db", "db_path", type=click.Path(), default=None)
@click.option("--all", "ingest_all", is_flag=True, help="Run all ingest steps")
@click.option("--ensembl", "ensembl_flag", is_flag=True)
@click.option("--kegg", "kegg_flag", is_flag=True)
@click.option("--string", "string_flag", is_flag=True)
@click.option("--salmotif", "salmotif_flag", is_flag=True)
@click.option("--mirna", "mirna_flag", is_flag=True)
@click.option("--limit-genes", type=int, default=None, help="Limit Ensembl genes (testing)")
@click.option("--string-limit", type=int, default=500, help="Max protein-coding genes for STRING")
def ingest(
    db_path: str | None,
    ingest_all: bool,
    ensembl_flag: bool,
    kegg_flag: bool,
    string_flag: bool,
    salmotif_flag: bool,
    mirna_flag: bool,
    limit_genes: int | None,
    string_limit: int,
) -> None:
    """Download / transform upstream data into SQLite."""
    conn = connect(db_path)
    if not table_exists(conn, "genes"):
        init_schema(conn)

    do_all = ingest_all or not any(
        [ensembl_flag, kegg_flag, string_flag, salmotif_flag, mirna_flag]
    )
    if do_all:
        ensembl_flag = kegg_flag = string_flag = salmotif_flag = mirna_flag = True

    if ensembl_flag:
        click.echo("Ingesting Ensembl (BioMart)...")
        ensembl.run(conn, limit_genes=limit_genes)
    if kegg_flag:
        click.echo("Ingesting KEGG pathways...")
        kegg_ingest.load_pathways_and_members(conn, immune_only=True)
    if string_flag:
        click.echo("Ingesting STRING interactions...")
        string_db.load_interactions_for_genes(conn, limit_genes=string_limit)
    if salmotif_flag:
        click.echo("Ingesting SalMotif / placeholder TFBS...")
        salmotifdb.run(conn)
    if mirna_flag:
        click.echo("Ingesting miRNA features...")
        mirna.run(conn)
    click.echo("Ingest complete.")


@main.command()
@click.argument("identifier")
@click.option("--db", "db_path", type=click.Path(), default=None)
@click.option("--json-out", is_flag=True)
def gene(identifier: str, db_path: str | None, json_out: bool) -> None:
    """Query a gene by symbol or Ensembl id."""
    conn = connect(db_path)
    res = query_gene(conn, identifier)
    if json_out:
        out = {k: v.to_dict(orient="records") for k, v in res.items() if isinstance(v, pd.DataFrame)}
        click.echo(json.dumps(out, default=str, indent=2))
        return
    for name, df in res.items():
        if name.startswith("_"):
            continue
        click.echo(f"\n== {name} ==")
        if hasattr(df, "empty") and df.empty:
            click.echo("(empty)")
        else:
            click.echo(df.to_string(index=False))


@main.command()
@click.argument("text")
@click.option("--db", "db_path", type=click.Path(), default=None)
def process(text: str, db_path: str | None) -> None:
    """Search physiological processes / pathways by name or id."""
    conn = connect(db_path)
    res = query_process(conn, text)
    click.echo(res["processes"].to_string(index=False))
    click.echo("\n== genes ==")
    click.echo(res["genes"].to_string(index=False))


@main.command()
@click.option("--gene-id", "ensembl_gene_id", required=True, help="Ensembl gene id")
@click.option("--db", "db_path", type=click.Path(), default=None)
@click.option("--limit", type=int, default=15)
def crispr(ensembl_gene_id: str, db_path: str | None, limit: int) -> None:
    """Design and store SpCas9 guides for a gene (requires Ensembl REST)."""
    conn = connect(db_path)
    n = ingest_crispr.save_guides_for_gene(conn, ensembl_gene_id, max_guides=limit)
    click.echo(f"Stored {n} guide candidates for {ensembl_gene_id}")


@main.command("crispr-list")
@click.argument("identifier")
@click.option("--db", "db_path", type=click.Path(), default=None)
def crispr_list(identifier: str, db_path: str | None) -> None:
    """List stored CRISPR guides for a gene symbol or id."""
    conn = connect(db_path)
    df = list_crispr_for_gene(conn, identifier)
    click.echo(df.to_string(index=False))


@main.command("crispr-pathway")
@click.argument("pathway_query")
@click.option("--db", "db_path", type=click.Path(), default=None)
def crispr_pathway(pathway_query: str, db_path: str | None) -> None:
    """Rank genes in matching pathways by stored CRISPR metrics."""
    conn = connect(db_path)
    df = rank_genes_for_pathway_crispr(conn, pathway_query)
    click.echo(df.to_string(index=False))


@main.command("subgraph")
@click.argument("identifier")
@click.option("--db", "db_path", type=click.Path(), default=None)
@click.option("--depth", type=int, default=1)
def subgraph_cmd(identifier: str, db_path: str | None, depth: int) -> None:
    """Print regulatory subgraph nodes/edges (text)."""
    conn = connect(db_path)
    G = regulatory_subgraph(conn, identifier, depth=depth)
    click.echo(f"Nodes: {G.number_of_nodes()} Edges: {G.number_of_edges()}")
    for u, v, data in G.edges(data=True):
        click.echo(f"{u} -> {v} {data}")


if __name__ == "__main__":
    main()
