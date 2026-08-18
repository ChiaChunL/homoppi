import logging
import math

import pandas as pd
import pytest

from homoppi.ppidb import PPIIndex, build_ppi_component, lib_proteins_by_taxid

logger = logging.getLogger("test")


def write_lib(tmp_path, text):
    path = tmp_path / "ppis.tsv"
    path.write_text(text)
    return path


def test_canonicalization_and_dedup(tmp_path):
    lib = write_lib(
        tmp_path,
        "protein_a\tprotein_b\ttaxid\tscore\n"
        "P2\tP1\t9606\t0.5\n"      # unsorted -> (P1, P2)
        "P1\tP2\t9606\t0.8\n"      # duplicate -> max score kept
        "P3\tP3\t9606\t0.4\n"      # homodimer template kept
        "Q1\tQ2\t10090\t\n",       # unscored -> NaN
    )
    out = tmp_path / "db" / "ppi"
    stats = build_ppi_component(lib, out, logger)
    assert stats["template_pairs"] == 3
    assert stats["homodimer_templates"] == 1

    df = pd.read_csv(out / "templates.tsv", sep="\t")
    row = df[(df.protein_a == "P1") & (df.protein_b == "P2")].iloc[0]
    assert row.score == pytest.approx(0.8)


def test_score_out_of_range_rejected(tmp_path):
    lib = write_lib(tmp_path, "protein_a\tprotein_b\ttaxid\tscore\nP1\tP2\t9606\t1.5\n")
    with pytest.raises(ValueError, match="outside"):
        build_ppi_component(lib, tmp_path / "db" / "ppi", logger)


def test_missing_column_rejected(tmp_path):
    lib = write_lib(tmp_path, "protein_a\tprotein_b\nP1\tP2\n")
    with pytest.raises(ValueError, match="missing required columns"):
        build_ppi_component(lib, tmp_path / "db" / "ppi", logger)


@pytest.fixture()
def built_db(tmp_path):
    lib = write_lib(
        tmp_path,
        "protein_a\tprotein_b\ttaxid\tscore\n"
        "P1\tP2\t9606\t0.8\n"
        "Q1\tQ2\t10090\t0.6\n"
        "P3\tP3\t9606\t\n",
    )
    build_ppi_component(lib, tmp_path / "db" / "ppi", logger)
    return tmp_path / "db"


def test_index_load_and_adjacency(built_db):
    index = PPIIndex.load(built_db)
    assert index.pairs[("P1", "P2")] == (0.8, 9606)
    assert index.adj["P1"] == {"P2"}
    assert index.adj["P3"] == {"P3"}  # homodimer self-adjacency
    raw, _ = index.pairs[("P3", "P3")]
    assert math.isnan(raw)
    assert index.score_of(("P3", "P3"), default=0.3) == pytest.approx(0.3)


def test_index_taxid_filter(built_db):
    index = PPIIndex.load(built_db, taxids={10090})
    assert set(index.pairs) == {("Q1", "Q2")}


def test_lib_proteins_by_taxid(built_db):
    proteins = lib_proteins_by_taxid(built_db / "ppi")
    assert proteins[9606] == {"P1", "P2", "P3"}
    assert proteins[10090] == {"Q1", "Q2"}
