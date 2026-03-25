"""Tests for query helpers (offline demo DB)."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from salmo_omlas.db import connect, init_schema
from salmo_omlas.demo_seed import seed
from salmo_omlas.query.core import query_gene, query_process, regulatory_subgraph
from salmo_omlas.query.crispr_queries import list_crispr_for_gene, rank_genes_for_pathway_crispr


@pytest.fixture()
def demo_conn():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "t.db"
        conn = connect(p)
        init_schema(conn)
        seed(conn)
        yield conn
        conn.close()


def test_query_gene_tlr3(demo_conn):
    res = query_gene(demo_conn, "TLR3")
    assert not res["gene"].empty
    assert res["gene"].iloc[0]["symbol"] == "TLR3"
    assert not res["processes"].empty


def test_query_process_immune(demo_conn):
    res = query_process(demo_conn, "immune")
    assert not res["processes"].empty
    assert not res["genes"].empty


def test_subgraph(demo_conn):
    G = regulatory_subgraph(demo_conn, "STAT1", depth=2)
    assert G.number_of_nodes() >= 1


def test_crispr_list(demo_conn):
    df = list_crispr_for_gene(demo_conn, "TLR3")
    assert not df.empty
