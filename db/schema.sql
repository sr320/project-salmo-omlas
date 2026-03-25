-- Salmon Gene Regulatory Interaction Database schema
-- SQLite 3

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS species (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  taxon_id INTEGER NOT NULL UNIQUE,
  scientific_name TEXT NOT NULL,
  common_name TEXT,
  assembly TEXT
);

CREATE TABLE IF NOT EXISTS sources (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  name TEXT NOT NULL,
  version TEXT,
  url TEXT,
  downloaded_at TEXT
);

CREATE TABLE IF NOT EXISTS genes (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  ensembl_gene_id TEXT NOT NULL UNIQUE,
  symbol TEXT,
  name TEXT,
  biotype TEXT,
  chromosome TEXT,
  start INTEGER,
  end INTEGER,
  strand INTEGER,
  description TEXT,
  species_id INTEGER NOT NULL REFERENCES species (id),
  source_id INTEGER REFERENCES sources (id)
);

CREATE INDEX IF NOT EXISTS idx_genes_symbol ON genes (symbol);
CREATE INDEX IF NOT EXISTS idx_genes_biotype ON genes (biotype);
CREATE INDEX IF NOT EXISTS idx_genes_chr_range ON genes (chromosome, start, end);

CREATE TABLE IF NOT EXISTS physiological_processes (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  term_id TEXT NOT NULL UNIQUE,
  name TEXT NOT NULL,
  namespace TEXT NOT NULL, -- go_bp, kegg_pathway, custom
  parent_term_id TEXT,
  definition TEXT
);

CREATE INDEX IF NOT EXISTS idx_processes_namespace ON physiological_processes (namespace);

CREATE TABLE IF NOT EXISTS gene_processes (
  gene_id INTEGER NOT NULL REFERENCES genes (id) ON DELETE CASCADE,
  process_id INTEGER NOT NULL REFERENCES physiological_processes (id) ON DELETE CASCADE,
  evidence TEXT,
  PRIMARY KEY (gene_id, process_id)
);

CREATE INDEX IF NOT EXISTS idx_gene_processes_process ON gene_processes (process_id);

CREATE TABLE IF NOT EXISTS orthologs (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  gene_id INTEGER NOT NULL REFERENCES genes (id) ON DELETE CASCADE,
  ortholog_ensembl_id TEXT NOT NULL,
  ortholog_symbol TEXT,
  ortholog_name TEXT,
  species_taxon_id INTEGER NOT NULL,
  ortholog_type TEXT,
  identity_percent REAL,
  confidence REAL,
  source_id INTEGER REFERENCES sources (id),
  UNIQUE (gene_id, ortholog_ensembl_id, species_taxon_id)
);

CREATE INDEX IF NOT EXISTS idx_orthologs_taxon ON orthologs (species_taxon_id);
CREATE INDEX IF NOT EXISTS idx_orthologs_ensembl ON orthologs (ortholog_ensembl_id);

CREATE TABLE IF NOT EXISTS regulatory_elements (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  element_id TEXT NOT NULL UNIQUE,
  feature_type TEXT NOT NULL, -- mirna, lncrna, enhancer, promoter, tfbs, other
  chromosome TEXT,
  start INTEGER,
  end INTEGER,
  strand INTEGER,
  linked_gene_id INTEGER REFERENCES genes (id) ON DELETE SET NULL,
  motif_name TEXT,
  metadata TEXT, -- JSON as text
  source_id INTEGER REFERENCES sources (id)
);

CREATE INDEX IF NOT EXISTS idx_regel_type ON regulatory_elements (feature_type);
CREATE INDEX IF NOT EXISTS idx_regel_chr ON regulatory_elements (chromosome, start, end);

CREATE TABLE IF NOT EXISTS interactions (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  regulator_gene_id INTEGER REFERENCES genes (id) ON DELETE CASCADE,
  regulator_element_id INTEGER REFERENCES regulatory_elements (id) ON DELETE CASCADE,
  target_gene_id INTEGER NOT NULL REFERENCES genes (id) ON DELETE CASCADE,
  interaction_type TEXT NOT NULL, -- tf_to_gene, mirna_repression, ppi, enhancer_to_gene, other
  evidence_score REAL,
  direction TEXT, -- activates, represses, binds, associated
  metadata TEXT,
  source_id INTEGER REFERENCES sources (id),
  CHECK (
    (regulator_gene_id IS NOT NULL AND regulator_element_id IS NULL)
    OR (regulator_gene_id IS NULL AND regulator_element_id IS NOT NULL)
  )
);

CREATE INDEX IF NOT EXISTS idx_inter_target ON interactions (target_gene_id);
CREATE INDEX IF NOT EXISTS idx_inter_regulator_gene ON interactions (regulator_gene_id);
CREATE INDEX IF NOT EXISTS idx_inter_type ON interactions (interaction_type);

CREATE TABLE IF NOT EXISTS crispr_targets (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  gene_id INTEGER NOT NULL REFERENCES genes (id) ON DELETE CASCADE,
  guide_sequence TEXT NOT NULL,
  pam TEXT NOT NULL DEFAULT 'NGG',
  strand TEXT,
  chromosome TEXT,
  start INTEGER,
  end INTEGER,
  on_target_score REAL,
  off_target_risk REAL, -- higher = more risk (heuristic)
  exon_rank INTEGER,
  notes TEXT,
  computed_at TEXT,
  UNIQUE (gene_id, guide_sequence, start, end)
);

CREATE INDEX IF NOT EXISTS idx_crispr_gene ON crispr_targets (gene_id);
