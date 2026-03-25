"""Paths, API endpoints, and species constants."""

from __future__ import annotations

import os
from pathlib import Path

# Project root (parent of src/)
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DB_DIR = PROJECT_ROOT / "db"
SCHEMA_PATH = DB_DIR / "schema.sql"
DEFAULT_DB_PATH = DB_DIR / "salmo_reg.db"

# Ensembl
ENSEMBL_REST = "https://rest.ensembl.org"
ENSEMBL_SPECIES = "salmo_salar"
ENSEMBL_ASSEMBLY = "Ssal_v3.1"
BIOMART_URL = "https://www.ensembl.org/biomart/martservice"

# Cross-species ortholog targets (Ensembl species names)
ORTHOLOG_SPECIES = {
    "danio_rerio": {"taxon_id": 7955, "display": "Zebrafish"},
    "mus_musculus": {"taxon_id": 10090, "display": "Mouse"},
    "homo_sapiens": {"taxon_id": 9606, "display": "Human"},
}

# KEGG
KEGG_BASE = "https://rest.kegg.jp"
KEGG_ORG_SASA = "sasa"  # Salmo salar

# STRING
STRING_API = "https://string-db.org/api"
STRING_SPECIES_NCBI_TAXON = 8030  # Salmo salar

# miRBase
# miRBase (FASTA). Genome-coordinate GFF3 for salmon is not available here,
# so ingestion stores miRNAs without genomic coordinates.
MIRBASE_MATURE_FASTA = "https://www.mirbase.org/download/mature.fa"

# Immune-focused GO roots (optional filter)
IMMUNE_GO_ROOTS = (
    "GO:0002376",  # immune system process
    "GO:0006955",  # immune response
    "GO:0045087",  # innate immune response
    "GO:0006954",  # inflammatory response
    "GO:0051607",  # defense response to virus
)

# KEGG immune-related pathway id suffixes (sasa + id)
IMMUNE_KEGG_PATHWAY_IDS = (
    "04620",
    "04630",
    "04060",
    "04621",
    "04622",
)


def get_db_path() -> Path:
    return Path(os.environ.get("SALMO_OMLAS_DB", DEFAULT_DB_PATH))


def ensure_dirs() -> None:
    DATA_RAW.mkdir(parents=True, exist_ok=True)
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    DB_DIR.mkdir(parents=True, exist_ok=True)
