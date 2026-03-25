# Salmon Gene Regulatory Interaction Database (Salmo-OMLAS)

A Python + SQLite database of gene regulatory interactions for Atlantic salmon (*Salmo salar*), integrating coding and non-coding elements, physiological process annotations (GO, KEGG), cross-species orthologs (zebrafish, mouse, human), and CRISPR-oriented target metadata.

## Features

- **Query by gene** — gene metadata, processes, interactions, orthologs, regulatory elements
- **Query by physiological process** — GO biological process and KEGG pathway membership
- **CRISPR use case** — candidate SpCas9 gRNAs from exon sequences with simple on-target scoring
- **Interfaces** — Jupyter notebooks, CLI (`salmo-omlas`), Streamlit web app

## Setup

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
pip install -e .
```

Initialize the database schema:

```bash
pip install -e .
salmo-omlas init-db
# or: PYTHONPATH=src python -m salmo_omlas.cli init-db
```

Load offline demo rows (no network):

```bash
salmo-omlas demo
```

Run ingestion (requires network; may take a while for full Ensembl/STRING):

```bash
salmo-omlas ingest --all
# Or step by step:
salmo-omlas ingest --ensembl --kegg --string
```

## Usage

### CLI

```bash
salmo-omlas gene STAT1
salmo-omlas process "interferon"
salmo-omlas crispr --gene-id ENSG00000123456 --limit 10
```

### Streamlit

```bash
streamlit run src/salmo_omlas/web.py
# from repo root; ensure `pip install -e .` or PYTHONPATH=src
```

### Python

```python
from salmo_omlas.query.core import query_gene, query_process
from salmo_omlas.db import connect

conn = connect()
df = query_gene(conn, "STAT1")
```

## Data sources

| Source | Content |
|--------|---------|
| Ensembl | Genes, GO, orthologs |
| KEGG | Pathways (organism `sasa`) |
| STRING | Protein–protein associations |
| SalMotifDB / motifs | TF motifs (placeholder or downloaded) |
| miRBase | miRNA metadata |

## License

Project data are subject to each upstream database’s terms of use.
