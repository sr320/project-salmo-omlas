"""Minimal deterministic data for offline tests and demos."""

from __future__ import annotations

import sqlite3

from salmo_omlas.ingest._util import dumps_metadata, now_iso


def seed(conn: sqlite3.Connection) -> None:
    cur = conn.cursor()
    cur.execute(
        "INSERT OR IGNORE INTO sources (name, version, url, downloaded_at) VALUES (?, ?, ?, ?)",
        ("demo", "0", "local", now_iso()),
    )
    conn.commit()
    sid = cur.execute("SELECT id FROM sources WHERE name='demo' ORDER BY id DESC LIMIT 1").fetchone()[0]

    cur.execute(
        """INSERT OR IGNORE INTO species (taxon_id, scientific_name, common_name, assembly)
           VALUES (8030, 'Salmo salar', 'Atlantic salmon', 'Ssal_v3.1')"""
    )
    conn.commit()
    sp = cur.execute("SELECT id FROM species WHERE taxon_id=8030").fetchone()[0]

    demo_genes = [
        (
            "ENSSSAG00000000001",
            "TLR3",
            "TLR3",
            "protein_coding",
            "1",
            1000,
            5000,
            1,
            "Toll-like receptor 3",
            sp,
            sid,
        ),
        (
            "ENSSSAG00000000002",
            "STAT1",
            "STAT1",
            "protein_coding",
            "1",
            6000,
            12000,
            1,
            "Signal transducer and activator of transcription 1",
            sp,
            sid,
        ),
        (
            "ENSSSAG00000000003",
            "IRF3",
            "IRF3",
            "protein_coding",
            "2",
            2000,
            8000,
            -1,
            "Interferon regulatory factor 3",
            sp,
            sid,
        ),
    ]
    for row in demo_genes:
        cur.execute(
            """INSERT OR IGNORE INTO genes (
                 ensembl_gene_id, symbol, name, biotype, chromosome, start, end, strand,
                 description, species_id, source_id
               ) VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
            row,
        )
    conn.commit()

    cur.execute(
        """INSERT OR IGNORE INTO physiological_processes (term_id, name, namespace, parent_term_id, definition)
           VALUES ('GO:0006955', 'immune response', 'go_bp', NULL, NULL)"""
    )
    cur.execute(
        """INSERT OR IGNORE INTO physiological_processes (term_id, name, namespace, parent_term_id, definition)
           VALUES ('sasa04620', 'Toll-like receptor signaling pathway', 'kegg_pathway', NULL, NULL)"""
    )
    conn.commit()

    gmap = dict(cur.execute("SELECT ensembl_gene_id, id FROM genes").fetchall())
    p_immune = cur.execute(
        "SELECT id FROM physiological_processes WHERE term_id='GO:0006955'"
    ).fetchone()[0]
    p_kegg = cur.execute(
        "SELECT id FROM physiological_processes WHERE term_id='sasa04620'"
    ).fetchone()[0]

    for ens in demo_genes:
        gid = gmap[ens[0]]
        cur.execute(
            "INSERT OR REPLACE INTO gene_processes (gene_id, process_id, evidence) VALUES (?,?,?)",
            (gid, p_immune, "demo"),
        )
        if ens[0] in ("ENSSSAG00000000001", "ENSSSAG00000000003"):
            cur.execute(
                "INSERT OR REPLACE INTO gene_processes (gene_id, process_id, evidence) VALUES (?,?,?)",
                (gid, p_kegg, "demo"),
            )

    cur.execute(
        """INSERT OR REPLACE INTO orthologs (
             gene_id, ortholog_ensembl_id, ortholog_symbol, ortholog_name,
             species_taxon_id, ortholog_type, identity_percent, confidence, source_id
           ) VALUES (?,?,?,?,?,?,?,?,?)""",
        (
            gmap["ENSSSAG00000000001"],
            "ENSDARG00000000001",
            "tlr3",
            "tlr3",
            7955,
            "ortholog_one2one",
            72.0,
            0.8,
            sid,
        ),
    )
    conn.commit()

    cur.execute(
        """INSERT OR IGNORE INTO regulatory_elements (
             element_id, feature_type, chromosome, start, end, strand,
             linked_gene_id, motif_name, metadata, source_id
           ) VALUES (?,?,?,?,?,?,?,?,?,?)""",
        (
            "demo_mir_ssa-mir-1",
            "mirna",
            "1",
            500,
            580,
            1,
            gmap["ENSSSAG00000000001"],
            "ssa-mir-1",
            dumps_metadata({}),
            sid,
        ),
    )
    conn.commit()
    el_id = cur.execute(
        "SELECT id FROM regulatory_elements WHERE element_id='demo_mir_ssa-mir-1'"
    ).fetchone()[0]

    cur.execute(
        """INSERT OR REPLACE INTO interactions (
             regulator_gene_id, regulator_element_id, target_gene_id,
             interaction_type, evidence_score, direction, metadata, source_id
           ) VALUES (?,?,?,?,?,?,?,?)""",
        (
            None,
            el_id,
            gmap["ENSSSAG00000000001"],
            "mirna_repression",
            0.4,
            "represses",
            dumps_metadata({}),
            sid,
        ),
    )
    cur.execute(
        """INSERT OR REPLACE INTO interactions (
             regulator_gene_id, regulator_element_id, target_gene_id,
             interaction_type, evidence_score, direction, metadata, source_id
           ) VALUES (?,?,?,?,?,?,?,?)""",
        (
            gmap["ENSSSAG00000000003"],
            None,
            gmap["ENSSSAG00000000002"],
            "tf_to_gene",
            0.55,
            "activates",
            dumps_metadata({}),
            sid,
        ),
    )
    conn.commit()

    cur.execute(
        """INSERT OR REPLACE INTO crispr_targets (
             gene_id, guide_sequence, pam, strand, chromosome, start, end,
             on_target_score, off_target_risk, exon_rank, notes, computed_at
           ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""",
        (
            gmap["ENSSSAG00000000001"],
            "GCTAGCTAGCTAGCTAGCTA",
            "TGG",
            "+",
            "1",
            10,
            33,
            0.82,
            0.1,
            1,
            "demo guide",
            now_iso(),
        ),
    )
    conn.commit()
