"""End-to-end test using real makeblastdb/blastp; skipped when BLAST+ is not installed."""

import logging
import shutil

import pandas as pd
import pytest

from homoppi.blast import ensure_blast_stage
from homoppi.config import BlastParams
from homoppi.interolog import infer_pairs_mode, load_homologs
from homoppi.makedb import build_database
from homoppi.ppidb import PPIIndex
from homoppi.workdir import Workdir

pytestmark = pytest.mark.skipif(
    shutil.which("blastp") is None or shutil.which("makeblastdb") is None,
    reason="BLAST+ not installed",
)

logger = logging.getLogger("test")

# Two 90-aa template proteins and queries identical to them.
SEQ1 = "MKVLLAGGSTRRAAEELGVSQPAVSKWLNGGSVPSAENLLALSKLLGVSLDELVFGNRKTGDLLEQVRALPEDKQEEVLDYIDFLRQKR"
SEQ2 = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNE"


@pytest.fixture()
def db_and_query(tmp_path):
    ppi = tmp_path / "ppis.tsv"
    ppi.write_text("protein_a\tprotein_b\ttaxid\tscore\nT1\tT2\t9606\t0.8\n")
    templates = tmp_path / "human.fasta"
    templates.write_text(f">sp|T1|X_HUMAN\n{SEQ1}\n>sp|T2|Y_HUMAN\n{SEQ2}\n")
    db_dir = tmp_path / "db"
    build_database(db_dir, ppi, [f"9606={templates}"], logger)

    query = tmp_path / "query.fasta"
    query.write_text(f">A\n{SEQ1}\n>B\n{SEQ2}\n")
    return db_dir, query


def test_end_to_end_pairs(db_and_query, tmp_path):
    db_dir, query = db_and_query
    wd = Workdir(tmp_path / "run")
    homologs_path = ensure_blast_stage(wd, db_dir, query, BlastParams(), logger)

    homologs = load_homologs(homologs_path)
    assert "T1" in homologs["A"] and "T2" in homologs["B"]

    index = PPIIndex.load(db_dir)
    summary = tmp_path / "summary.tsv"
    infer_pairs_mode(
        [("A", "B")], homologs, index, summary, None,
        include_self=False, default_score=0.0, logger=logger,
    )
    df = pd.read_csv(summary, sep="\t")
    assert df.iloc[0].s_im == pytest.approx(0.8, abs=1e-4)
    assert df.iloc[0].best_template_a == "T1"


def test_blast_stage_resume(db_and_query, tmp_path):
    db_dir, query = db_and_query
    wd = Workdir(tmp_path / "run")
    params = BlastParams()
    first = ensure_blast_stage(wd, db_dir, query, params, logger)
    mtime = first.stat().st_mtime_ns
    second = ensure_blast_stage(wd, db_dir, query, params, logger)  # should be a no-op
    assert second.stat().st_mtime_ns == mtime


def test_blast_tightening_refilters_cached_raw(db_and_query, tmp_path):
    db_dir, query = db_and_query
    wd = Workdir(tmp_path / "run")
    ensure_blast_stage(wd, db_dir, query, BlastParams(), logger)
    raw_mtime = (wd.blast_dir / "raw.fmt6.tsv").stat().st_mtime_ns

    tight = BlastParams(evalue=1e-30, min_identity=90.0)
    out = ensure_blast_stage(wd, db_dir, query, tight, logger)

    assert (wd.blast_dir / "raw.fmt6.tsv").stat().st_mtime_ns == raw_mtime  # no blastp rerun
    homologs = pd.read_csv(out, sep="\t")
    assert not homologs.empty  # identical sequences still pass the strict thresholds
    assert (homologs.pident >= 90).all()
