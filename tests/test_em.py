import logging

import pytest

from homoppi.em import EMParams, derive_false_pos_rate, load_domain_annotations, run_ddi_em

logger = logging.getLogger("test")


def test_derive_false_pos_rate_positive():
    fp = derive_false_pos_rate(n_ppis=1000, n_proteins=500, false_neg_rate=0.8)
    assert 0 < fp < 1


def test_load_domain_annotations_merges_files(tmp_path):
    header = "query_id\tpfam_acc\tdomain_name\tcevalue\tali_from\tali_to\n"
    f1 = tmp_path / "a.tsv"
    f1.write_text(header + "P1\tPF00001\tX\t1e-20\t1\t50\n")
    f2 = tmp_path / "b.tsv"
    f2.write_text(header + "P1\tPF00002\tY\t1e-20\t60\t100\nP2\tPF00003\tZ\t1e-20\t1\t50\n")
    annotations = load_domain_annotations([f1, f2])
    assert annotations == {"P1": {"PF00001", "PF00002"}, "P2": {"PF00003"}}


@pytest.fixture()
def toy_universe():
    taxid = {p: 9606 for p in ("P1", "P2", "P3", "P4", "P5", "P6")}
    domains = {
        "P1": {"D1"}, "P3": {"D1"},
        "P2": {"D2"}, "P4": {"D2"},
        "P5": {"D3"}, "P6": {"D4"},
    }
    # (D1, D2) supported by 2 of 4 containing pairs; (D3, D4) by its only pair.
    ppis = {("P1", "P2"), ("P3", "P4"), ("P5", "P6")}
    return ppis, taxid, domains


def test_em_scores_are_probabilities_and_ordered(toy_universe):
    ppis, taxid, domains = toy_universe
    result = run_ddi_em(ppis, taxid, domains, EMParams(false_neg_rate=0.1, max_iter=30), logger)
    scores = result.set_index(["pfam_a", "pfam_b"]).em_score

    assert ((scores >= 0) & (scores <= 1)).all()
    # a domain pair whose every containing protein pair interacts outranks a diluted one
    assert scores[("D3", "D4")] > scores[("D1", "D2")]

    row = result.set_index(["pfam_a", "pfam_b"]).loc[("D1", "D2")]
    assert row.n_protein_pairs == 4  # P1P2, P1P4, P3P2, P3P4 (same species)
    assert row.n_interacting_pairs == 2


def test_em_promiscuous_domain_excluded(toy_universe):
    ppis, taxid, domains = toy_universe
    params = EMParams(false_neg_rate=0.1, max_iter=5, max_proteins_per_domain=1)
    result = run_ddi_em(ppis, taxid, domains, params, logger)
    # D1 and D2 each occur in 2 proteins > cap 1 -> (D1, D2) excluded
    keys = set(map(tuple, result[["pfam_a", "pfam_b"]].values))
    assert ("D1", "D2") not in keys
    assert ("D3", "D4") in keys


def test_em_requires_annotated_ppis():
    with pytest.raises(ValueError, match="no known PPIs"):
        run_ddi_em({("A", "B")}, {}, {}, EMParams(), logger)
