import logging

import pandas as pd
import pytest

from homoppi.pipeline import DDI_SUMMARY, IM_SUMMARY, merge_stage
from homoppi.workdir import Workdir

logger = logging.getLogger("test")


@pytest.fixture()
def wd_with_summaries(tmp_path):
    wd = Workdir(tmp_path / "run")
    (wd.results_dir / IM_SUMMARY).write_text(
        "query_a\tquery_b\tn_templates\ts_im\tbest_template_a\tbest_template_b\tbest_template_taxid\tbest_template_score\n"
        "A\tB\t3\t0.96\tP1\tP2\t9606\t0.8\n"
        "A\tC\t1\t0.5\tP1\tP2\t9606\t0.5\n"
    )
    (wd.results_dir / DDI_SUMMARY).write_text(
        "query_a\tquery_b\tn_templates\ts_ddi\tbest_template_a\tbest_template_b\tbest_template_source\tbest_template_score\n"
        "A\tB\t2\t0.9\tPF1\tPF2\t3did\t0.8\n"
        "B\tC\t1\t0.4\tPF1\tPF2\tem\t0.4\n"
    )
    return wd


def test_merge_outer_join_with_zero_fill(wd_with_summaries):
    out = merge_stage(wd_with_summaries, logger)
    df = pd.read_csv(out, sep="\t").set_index(["query_a", "query_b"])

    assert df.loc[("A", "B")].s_im == pytest.approx(0.96)
    assert df.loc[("A", "B")].s_ddi == pytest.approx(0.9)
    assert df.loc[("A", "C")].s_ddi == 0.0
    assert df.loc[("A", "C")].n_ddi_templates == 0
    assert df.loc[("B", "C")].s_im == 0.0
    assert "s_fused" not in df.columns


def test_merge_fused_column(wd_with_summaries):
    out = merge_stage(wd_with_summaries, logger, fused=True)
    df = pd.read_csv(out, sep="\t").set_index(["query_a", "query_b"])
    assert df.loc[("A", "B")].s_fused == pytest.approx(1 - (1 - 0.96) * (1 - 0.9), abs=1e-4)
    assert df.loc[("A", "C")].s_fused == pytest.approx(0.5, abs=1e-4)


def test_merge_missing_summary_raises(tmp_path):
    wd = Workdir(tmp_path / "empty")
    with pytest.raises(FileNotFoundError, match="summary missing"):
        merge_stage(wd, logger)
