"""Streamlit web UI for Salmo-OMLAS."""

from __future__ import annotations

import sys
from pathlib import Path

import networkx as nx
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# Allow running as `streamlit run src/salmo_omlas/web.py`
_SRC = Path(__file__).resolve().parents[1]
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from salmo_omlas.db import connect, init_schema, table_exists
from salmo_omlas.demo_seed import seed as demo_seed
from salmo_omlas.query.core import query_gene, query_process, regulatory_subgraph
from salmo_omlas.query.crispr_queries import list_crispr_for_gene, rank_genes_for_pathway_crispr


def _layout_graph(G: nx.DiGraph) -> dict:
    return nx.spring_layout(G, seed=42, k=0.5)


def plot_network(G: nx.DiGraph) -> go.Figure:
    if G.number_of_nodes() == 0:
        fig = go.Figure()
        fig.update_layout(title="No nodes", height=400)
        return fig
    pos = _layout_graph(G)
    edge_x: list[float] = []
    edge_y: list[float] = []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        line=dict(width=1, color="#888"),
        hoverinfo="none",
        mode="lines",
    )
    node_x = [pos[n][0] for n in G.nodes()]
    node_y = [pos[n][1] for n in G.nodes()]
    node_text = [str(G.nodes[n].get("label", n)) for n in G.nodes()]
    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=node_text,
        textposition="top center",
        hovertext=[f"{n}: {G.nodes[n]}" for n in G.nodes()],
        hoverinfo="text",
        marker=dict(size=16, color="#1f77b4"),
    )
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        showlegend=False,
        margin=dict(l=10, r=10, t=30, b=10),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=500,
        title="Regulatory subgraph",
    )
    return fig


def main() -> None:
    st.set_page_config(page_title="Salmo-OMLAS", layout="wide")
    st.title("Salmo-OMLAS — Atlantic salmon regulatory database")

    db_default = Path(__file__).resolve().parents[2] / "db" / "salmo_reg.db"
    db_path = st.sidebar.text_input("Database path", value=str(db_default))

    conn = connect(db_path)
    if not table_exists(conn, "genes"):
        if st.sidebar.button("Initialize empty DB + demo data"):
            init_schema(conn)
            demo_seed(conn)
            st.success("Schema + demo loaded.")
            st.rerun()
    else:
        if st.sidebar.button("Reload demo seed (additive)"):
            demo_seed(conn)
            st.success("Demo seed applied.")
            st.rerun()

    tab_gene, tab_proc, tab_crispr = st.tabs(["Gene", "Physiological process", "CRISPR"])

    with tab_gene:
        gid = st.text_input("Gene symbol or Ensembl id", value="TLR3")
        depth = st.slider("Subgraph depth", 1, 3, 1)
        if st.button("Search gene"):
            res = query_gene(conn, gid)
            st.subheader("Gene record")
            st.dataframe(res["gene"], use_container_width=True)
            st.subheader("Processes")
            st.dataframe(res["processes"], use_container_width=True)
            st.subheader("Orthologs")
            st.dataframe(res["orthologs"], use_container_width=True)
            st.subheader("Regulatory elements (linked)")
            st.dataframe(res["regulatory_elements"], use_container_width=True)
            c1, c2 = st.columns(2)
            with c1:
                st.markdown("**Interactions out (as regulator)**")
                st.dataframe(res["interactions_out"], use_container_width=True)
            with c2:
                st.markdown("**Interactions in (as target)**")
                st.dataframe(res["interactions_in"], use_container_width=True)
            st.subheader("Stored CRISPR guides")
            st.dataframe(res["crispr_targets"], use_container_width=True)

            G = regulatory_subgraph(conn, gid, depth=depth)
            st.plotly_chart(plot_network(G), use_container_width=True)

    with tab_proc:
        pq = st.text_input("Process / pathway keyword", value="immune")
        if st.button("Search processes"):
            res = query_process(conn, pq)
            st.dataframe(res["processes"], use_container_width=True)
            st.dataframe(res["genes"], use_container_width=True)

    with tab_crispr:
        st.caption("List precomputed guides or rank pathway genes by stored scores.")
        cg = st.text_input("Gene for guide list", value="TLR3")
        if st.button("List CRISPR rows"):
            st.dataframe(list_crispr_for_gene(conn, cg), use_container_width=True)
        pw = st.text_input("Pathway query for ranking", value="immune")
        if st.button("Rank by pathway"):
            st.dataframe(rank_genes_for_pathway_crispr(conn, pw), use_container_width=True)


if __name__ == "__main__":
    main()
