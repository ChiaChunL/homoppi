"""Refilter-vs-rerun logic of the cached blast/hmmscan stages (no external tools needed)."""

import logging

import pandas as pd
import pytest
from test_hmmer import domtblout_line

from homoppi.blast import ensure_blast_stage
from homoppi.config import BlastParams, HmmscanParams
from homoppi.external import ExternalToolError
from homoppi.hmmer import ensure_hmmscan_stage
from homoppi.workdir import Workdir

logger = logging.getLogger("test")


def fake_blast_db(tmp_path):
    db = tmp_path / "db"
    (db / "blastdb").mkdir(parents=True)
    (db / "blastdb" / "templates.pin").write_text("")
    (db / "blastdb" / "protein2taxid.tsv").write_text("protein_id\ttaxid\nP1\t9606\n")
    return db


@pytest.fixture()
def blast_workdir(tmp_path):
    wd = Workdir(tmp_path / "run")
    fasta = tmp_path / "q.fasta"
    fasta.write_text(">A\nMKVLLA\n")
    raw = wd.blast_dir / "raw.fmt6.tsv"
    raw.write_text(
        "A\tP1\t90.0\t100\t90.0\t1e-30\t300.0\t200\t1\t180\n"
        "A\tP1\t80.0\t100\t90.0\t1e-06\t100.0\t200\t1\t150\n"
    )
    wd.mark_stage(
        "blast",
        {"search": {"evalue": 1e-5, "max_target_seqs": 500},
         "filter": BlastParams(evalue=1e-5).filter_params()},
        {"query_fasta": fasta},
    )
    return wd, fasta, raw


def test_blast_tighter_evalue_refilters_without_rerun(blast_workdir, tmp_path):
    wd, fasta, raw = blast_workdir
    db = fake_blast_db(tmp_path)
    raw_mtime = raw.stat().st_mtime_ns

    out = ensure_blast_stage(wd, db, fasta, BlastParams(evalue=1e-20), logger)

    assert raw.stat().st_mtime_ns == raw_mtime  # blastp did not rerun
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 1 and df.iloc[0].evalue == pytest.approx(1e-30)
    # the recorded search E-value stays that of the raw file, allowing further refilters
    assert wd.get_stage("blast")["params"]["search"]["evalue"] == 1e-5


def test_blast_loosening_forces_rerun(blast_workdir, tmp_path):
    wd, fasta, _ = blast_workdir
    db = fake_blast_db(tmp_path)
    # 1e-3 is looser than the raw run's 1e-5: raw is incomplete, blastp must rerun,
    # which fails here (missing binary or fake database) - proving the rerun path.
    with pytest.raises(ExternalToolError):
        ensure_blast_stage(wd, db, fasta, BlastParams(evalue=1e-3), logger)


def test_blast_max_target_seqs_change_forces_rerun(blast_workdir, tmp_path):
    wd, fasta, _ = blast_workdir
    db = fake_blast_db(tmp_path)
    with pytest.raises(ExternalToolError):
        ensure_blast_stage(wd, db, fasta, BlastParams(evalue=1e-20, max_target_seqs=100), logger)


@pytest.fixture()
def hmmscan_workdir(tmp_path):
    wd = Workdir(tmp_path / "run")
    fasta = tmp_path / "q.fasta"
    fasta.write_text(">PROT1\nMKVLLA\n")
    raw = wd.hmmscan_dir / "raw.domtblout"
    raw.write_text(
        domtblout_line("DomA", "PF00001.1", "PROT1", "1e-25", seq_evalue="1e-30") + "\n"
        + domtblout_line("DomB", "PF00002.1", "PROT1", "1e-12", seq_evalue="1e-15") + "\n"
    )
    wd.mark_stage(
        "hmmscan",
        {"search": {"cevalue": 1e-10, "cut_tc": False},
         "filter": HmmscanParams(cevalue=1e-10).filter_params()},
        {"query_fasta": fasta},
    )
    return wd, fasta, raw


def test_hmmscan_tighter_cevalue_refilters_without_rerun(hmmscan_workdir, tmp_path):
    wd, fasta, raw = hmmscan_workdir
    db = tmp_path / "db"  # no pfam component: a rerun attempt would fail
    db.mkdir(exist_ok=True)
    raw_mtime = raw.stat().st_mtime_ns

    out = ensure_hmmscan_stage(wd, db, fasta, HmmscanParams(cevalue=1e-20), logger)

    assert raw.stat().st_mtime_ns == raw_mtime  # hmmscan did not rerun
    df = pd.read_csv(out, sep="\t")
    assert df.pfam_acc.tolist() == ["PF00001"]
    assert wd.get_stage("hmmscan")["params"]["search"]["cevalue"] == 1e-10


def test_hmmscan_cut_tc_switch_forces_rerun(hmmscan_workdir, tmp_path):
    wd, fasta, _ = hmmscan_workdir
    db = tmp_path / "db"
    db.mkdir(exist_ok=True)
    with pytest.raises(FileNotFoundError, match="no Pfam component"):
        ensure_hmmscan_stage(wd, db, fasta, HmmscanParams(cut_tc=True), logger)
