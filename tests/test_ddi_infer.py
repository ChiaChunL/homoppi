import logging

import pandas as pd
import pytest

from homoppi.ddi import DomHit, infer_pairs_mode, load_domains
from homoppi.ddidb import DDIIndex

logger = logging.getLogger("test")


def test_load_domains_keeps_best_instance(tmp_path):
    path = tmp_path / "domains.tsv"
    path.write_text(
        "query_id\tpfam_acc\tdomain_name\tcevalue\tali_from\tali_to\n"
        "X\tPF00069\tPkinase\t1e-20\t10\t90\n"
        "X\tPF00069\tPkinase\t1e-30\t200\t280\n"
        "Y\tPF00018\tSH3_1\t1e-15\t5\t60\n"
    )
    domains = load_domains(path)
    assert domains["X"]["PF00069"].cevalue == pytest.approx(1e-30)
    assert set(domains["Y"]) == {"PF00018"}


def test_ddi_pairs_mode(tmp_path):
    index = DDIIndex(
        pairs={("PF00018", "PF00069"): (0.8, "3did"), ("PF10417", "PF10417"): (0.5, "3did")},
        adj={"PF00018": {"PF00069"}, "PF00069": {"PF00018"}, "PF10417": {"PF10417"}},
    )
    domains = {
        "X": {"PF00069": DomHit(1e-30), "PF10417": DomHit(1e-12)},
        "Y": {"PF00018": DomHit(1e-15), "PF10417": DomHit(1e-9)},
        "Z": {"PF00019": DomHit(1e-10)},
    }
    summary_path = tmp_path / "s.tsv"
    evidence_path = tmp_path / "e.tsv"
    infer_pairs_mode(
        [("X", "Y"), ("X", "Z")], domains, index, summary_path, evidence_path,
        include_self=False, default_score=0.0, logger=logger,
    )
    summary = pd.read_csv(summary_path, sep="\t").set_index(["query_a", "query_b"])

    # X-Y: PF00069-PF00018 (0.8) and homotypic PF10417-PF10417 (0.5)
    assert summary.loc[("X", "Y")].n_templates == 2
    assert summary.loc[("X", "Y")].s_ddi == pytest.approx(1 - 0.2 * 0.5, abs=1e-4)
    assert summary.loc[("X", "Y")].best_template_source == "3did"
    assert summary.loc[("X", "Z")].n_templates == 0

    evidence = pd.read_csv(evidence_path, sep="\t")
    assert set(evidence.columns) >= {"template_a", "template_b", "source", "cevalue_a", "cevalue_b"}
    assert len(evidence) == 2
