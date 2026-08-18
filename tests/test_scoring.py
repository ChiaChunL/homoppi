import pytest

from homoppi.scoring import bayes_integration


def test_empty_evidence_scores_zero():
    assert bayes_integration([]) == 0.0


def test_single_template():
    assert bayes_integration([0.8]) == pytest.approx(0.8)


def test_bayesian_integration():
    assert bayes_integration([0.8, 0.6]) == pytest.approx(1 - 0.2 * 0.4)


def test_certain_template_dominates():
    assert bayes_integration([1.0, 0.1]) == pytest.approx(1.0)
