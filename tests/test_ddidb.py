import logging

import pandas as pd
import pytest

from homoppi.ddidb import DDIIndex, build_ddi_component, parse_3did_flat

logger = logging.getLogger("test")

THREE_DID = """\
#=ID\t1-cysPrx_C\t1-cysPrx_C\t(PF10417.9@Pfam\tPF10417.9@Pfam)
#=3D\t1prx\tA:40-50\tB:40-50\t1.0\t2\t0:0
//
#=ID\tPkinase\tSH3_1\t(PF00069.28@Pfam\tPF00018.28@Pfam)
//
"""


def test_parse_3did_flat(tmp_path):
    path = tmp_path / "3did_flat"
    path.write_text(THREE_DID)
    pairs = parse_3did_flat(path)
    assert pairs == {("PF10417", "PF10417"), ("PF00018", "PF00069")}


def write_em(tmp_path, rows):
    path = tmp_path / "em.tsv"
    path.write_text("pfam_a\tpfam_b\tem_score\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    return path


def test_build_scores_follow_paper_formula(tmp_path):
    three_did = tmp_path / "3did_flat"
    three_did.write_text(THREE_DID)
    em = write_em(tmp_path, [("PF00069", "PF00018", "0.6"), ("PF11111", "PF22222", "0.4")])

    build_ddi_component(tmp_path / "db" / "ddi", logger, three_did=three_did, em_scores=em)
    df = pd.read_csv(tmp_path / "db" / "ddi" / "templates.tsv", sep="\t").set_index(["pfam_a", "pfam_b"])

    # 3did only: 0.5*(1+0); both: 0.5*(1+0.6); EM only: 0.5*(0+0.4)
    assert df.loc[("PF10417", "PF10417")].score == pytest.approx(0.5)
    assert df.loc[("PF10417", "PF10417")].source == "3did"
    assert df.loc[("PF00018", "PF00069")].score == pytest.approx(0.8)
    assert df.loc[("PF00018", "PF00069")].source == "3did+em"
    assert df.loc[("PF11111", "PF22222")].score == pytest.approx(0.2)
    assert df.loc[("PF11111", "PF22222")].source == "em"


def test_build_plain_library(tmp_path):
    plain = tmp_path / "ddis.tsv"
    plain.write_text("pfam_a\tpfam_b\tscore\nPF00002\tPF00001\t0.9\n")
    build_ddi_component(tmp_path / "db" / "ddi", logger, plain=plain)
    index = DDIIndex.load(tmp_path / "db")
    assert index.pairs[("PF00001", "PF00002")] == (0.9, "custom")
    assert index.adj["PF00001"] == {"PF00002"}


def test_plain_and_3did_are_mutually_exclusive(tmp_path):
    plain = tmp_path / "ddis.tsv"
    plain.write_text("pfam_a\tpfam_b\nPF1\tPF2\n")
    with pytest.raises(ValueError, match="cannot be combined"):
        build_ddi_component(tmp_path / "db" / "ddi", logger, three_did=plain, plain=plain)


def test_em_score_out_of_range_rejected(tmp_path):
    em = write_em(tmp_path, [("PF1", "PF2", "1.5")])
    with pytest.raises(ValueError, match="outside"):
        build_ddi_component(tmp_path / "db" / "ddi", logger, em_scores=em)
