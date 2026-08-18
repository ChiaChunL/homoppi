"""End-to-end DDI test using real hmmbuild/hmmpress/hmmscan; skipped when HMMER is absent."""

import logging
import shutil
import subprocess

import pandas as pd
import pytest

from homoppi.config import HmmscanParams
from homoppi.ddi import infer_pairs_mode, load_domains
from homoppi.ddidb import DDIIndex
from homoppi.hmmer import ensure_hmmscan_stage
from homoppi.makedb import build_database
from homoppi.workdir import Workdir

pytestmark = pytest.mark.skipif(
    any(shutil.which(tool) is None for tool in ("hmmbuild", "hmmpress", "hmmscan")),
    reason="HMMER not installed",
)

logger = logging.getLogger("test")

SEQ_A = "MKVLLAGGSTRRAAEELGVSQPAVSKWLNGGSVPSAENLLALSKLLGVSLDELVFGNRKT"
SEQ_B = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRG"


def build_hmm(tmp_path, name, seq):
    fasta = tmp_path / f"{name}.fasta"
    fasta.write_text(f">{name}\n{seq}\n")
    hmm = tmp_path / f"{name}.hmm"
    subprocess.run(
        ["hmmbuild", "-n", name, "--amino", str(hmm), str(fasta)],
        check=True, capture_output=True,
    )
    return hmm.read_text()


@pytest.fixture()
def db_and_query(tmp_path):
    pfam = tmp_path / "toy_pfam.hmm"
    pfam.write_text(build_hmm(tmp_path, "DOMA", SEQ_A) + build_hmm(tmp_path, "DOMB", SEQ_B))

    ddis = tmp_path / "ddis.tsv"
    ddis.write_text("pfam_a\tpfam_b\tscore\nDOMA\tDOMB\t0.7\n")

    db_dir = tmp_path / "db"
    build_database(db_dir, None, None, logger, ddi_plain=ddis, pfam_hmm=pfam)

    query = tmp_path / "query.fasta"
    query.write_text(f">X\n{SEQ_A}\n>Y\n{SEQ_B}\n")
    return db_dir, query


def test_end_to_end_ddi(db_and_query, tmp_path):
    db_dir, query = db_and_query
    wd = Workdir(tmp_path / "run")
    domains_path = ensure_hmmscan_stage(wd, db_dir, query, HmmscanParams(), logger)

    domains = load_domains(domains_path)
    assert "DOMA" in domains["X"] and "DOMB" in domains["Y"]

    index = DDIIndex.load(db_dir)
    summary_path = tmp_path / "summary.tsv"
    infer_pairs_mode(
        [("X", "Y")], domains, index, summary_path, None,
        include_self=False, default_score=0.0, logger=logger,
    )
    df = pd.read_csv(summary_path, sep="\t")
    assert df.iloc[0].s_ddi == pytest.approx(0.7, abs=1e-4)
    assert df.iloc[0].best_template_a == "DOMA"


def test_hmmscan_stage_resume(db_and_query, tmp_path):
    db_dir, query = db_and_query
    wd = Workdir(tmp_path / "run")
    params = HmmscanParams()
    first = ensure_hmmscan_stage(wd, db_dir, query, params, logger)
    mtime = first.stat().st_mtime_ns
    second = ensure_hmmscan_stage(wd, db_dir, query, params, logger)  # should be a no-op
    assert second.stat().st_mtime_ns == mtime
