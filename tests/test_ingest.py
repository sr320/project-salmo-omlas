"""Light tests for pure helpers (no network)."""

from __future__ import annotations

from salmo_omlas.ingest.crispr import design_spcas9_guides


def test_design_guides_finds_ngg():
    cds = "A" * 10 + "GCTAGCTAGCTAGCTAGCTA" + "TGG" + "CCC" * 20
    guides = design_spcas9_guides(cds, max_guides=5)
    assert guides
    assert guides[0].pam.endswith("GG")
    assert len(guides[0].guide_sequence) == 20
