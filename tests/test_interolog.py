import logging
import math

import pandas as pd
import pytest

from homoppi.interolog import Hit, infer_all_mode, infer_pairs_mode
from homoppi.ppidb import PPIIndex

logger = logging.getLogger("test")

HIT = Hit(pident=95.0, qcov=90.0, evalue=1e-50)


@pytest.fixture()
def index():
    return PPIIndex(
        pairs={
            ("P1", "P2"): (0.8, 9606),
            ("Q1", "Q2"): (0.6, 10090),
            ("P3", "P3"): (0.5, 9606),  # homodimer template
        },
        adj={
            "P1": {"P2"}, "P2": {"P1"},
            "Q1": {"Q2"}, "Q2": {"Q1"},
            "P3": {"P3"},
        },
    )


@pytest.fixture()
def homologs():
    return {
        "A": {"P1": HIT, "Q1": HIT, "P3": HIT},
        "B": {"P2": HIT, "Q2": HIT, "P3": HIT},
        "C": {"P1": HIT},
    }


def run_pairs(pairs, homologs, index, tmp_path, **kwargs):
    summary = tmp_path / "summary.tsv"
    evidence = tmp_path / "evidence.tsv"
    kwargs.setdefault("include_self", False)
    kwargs.setdefault("default_score", 0.0)
    infer_pairs_mode(pairs, homologs, index, summary, evidence, logger=logger, **kwargs)
    return (
        pd.read_csv(summary, sep="\t"),
        pd.read_csv(evidence, sep="\t"),
    )


def test_pairs_mode_scores_and_evidence(homologs, index, tmp_path):
    summary, evidence = run_pairs([("A", "B"), ("A", "C")], homologs, index, tmp_path)

    row_ab = summary[(summary.query_a == "A") & (summary.query_b == "B")].iloc[0]
    # templates: P1-P2 (0.8), Q1-Q2 (0.6), P3-P3 (0.5) -> 1 - 0.2*0.4*0.5
    assert row_ab.n_templates == 3
    assert row_ab.s_im == pytest.approx(1 - 0.2 * 0.4 * 0.5, abs=1e-4)
    assert row_ab.best_template_a == "P1"
    assert row_ab.best_template_score == pytest.approx(0.8)

    row_ac = summary[(summary.query_a == "A") & (summary.query_b == "C")].iloc[0]
    assert row_ac.n_templates == 0
    assert row_ac.s_im == 0.0
    assert math.isnan(row_ac.best_template_taxid) or row_ac.best_template_taxid == ""

    assert len(evidence[(evidence.query_a == "A") & (evidence.query_b == "B")]) == 3


def test_self_pairs_skipped_by_default(homologs, index, tmp_path):
    summary, _ = run_pairs([("A", "A"), ("A", "B")], homologs, index, tmp_path)
    assert summary.query_a.tolist() == ["A"]
    assert summary.query_b.tolist() == ["B"]


def test_include_self_scores_homodimer(homologs, index, tmp_path):
    summary, _ = run_pairs([("A", "A")], homologs, index, tmp_path, include_self=True)
    row = summary.iloc[0]
    # only the P3-P3 homodimer template applies to A-A
    assert row.n_templates == 1
    assert row.s_im == pytest.approx(0.5, abs=1e-4)


def test_default_score_for_unscored_template(homologs, tmp_path):
    index = PPIIndex(pairs={("P1", "P2"): (math.nan, 9606)}, adj={"P1": {"P2"}, "P2": {"P1"}})
    summary, _ = run_pairs([("A", "B")], homologs, index, tmp_path, default_score=0.3)
    assert summary.iloc[0].s_im == pytest.approx(0.3, abs=1e-4)


def test_taxid_semantics_via_filtered_index(homologs, tmp_path):
    # Simulates --taxids 9606: only human templates remain in the index.
    index = PPIIndex(
        pairs={("P1", "P2"): (0.8, 9606), ("P3", "P3"): (0.5, 9606)},
        adj={"P1": {"P2"}, "P2": {"P1"}, "P3": {"P3"}},
    )
    summary, _ = run_pairs([("A", "B")], homologs, index, tmp_path)
    assert summary.iloc[0].n_templates == 2
    assert summary.iloc[0].s_im == pytest.approx(1 - 0.2 * 0.5, abs=1e-4)


def test_all_mode_enumerates_and_deduplicates(homologs, index, tmp_path):
    summary_path = tmp_path / "summary.tsv"
    evidence_path = tmp_path / "evidence.tsv"
    infer_all_mode(
        homologs, index, summary_path, evidence_path,
        include_self=False, default_score=0.0, min_score=0.0, logger=logger,
    )
    summary = pd.read_csv(summary_path, sep="\t").set_index(["query_a", "query_b"])

    # A-B: all three templates; B-C: via P1-P2 (C~P1, B~P2)
    assert summary.loc[("A", "B")].n_templates == 3
    assert summary.loc[("A", "B")].s_im == pytest.approx(1 - 0.2 * 0.4 * 0.5, abs=1e-4)
    assert summary.loc[("B", "C")].n_templates == 1
    assert summary.loc[("B", "C")].s_im == pytest.approx(0.8, abs=1e-4)
    # no self pairs (A-A via P3-P3 must be absent)
    assert ("A", "A") not in summary.index


def test_all_mode_min_score_filter(homologs, index, tmp_path):
    summary_path = tmp_path / "summary.tsv"
    infer_all_mode(
        homologs, index, summary_path, None,
        include_self=False, default_score=0.0, min_score=0.9, logger=logger,
    )
    summary = pd.read_csv(summary_path, sep="\t")
    assert summary.query_b.tolist() == ["B"]  # only A-B (0.96) survives
