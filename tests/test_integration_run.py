"""Full-pipeline test: makedb (all components) -> `homoppi run` -> combined summary.

Uses real makeblastdb/blastp/hmmbuild/hmmpress/hmmscan; skipped when missing.
"""

import logging
import shutil
import subprocess

import pandas as pd
import pytest
from typer.testing import CliRunner

from homoppi.cli import app
from homoppi.makedb import build_database

TOOLS = ("makeblastdb", "blastp", "hmmbuild", "hmmpress", "hmmscan")
pytestmark = pytest.mark.skipif(
    any(shutil.which(tool) is None for tool in TOOLS),
    reason="BLAST+/HMMER not installed",
)

logger = logging.getLogger("test")
runner = CliRunner()

SEQ_1 = "MKVLLAGGSTRRAAEELGVSQPAVSKWLNGGSVPSAENLLALSKLLGVSLDELVFGNRKTGDLLEQVRALPEDKQEEVLDYIDFLRQKR"
SEQ_2 = "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNE"


def build_hmm(tmp_path, name, seq):
    fasta = tmp_path / f"{name}.fasta"
    fasta.write_text(f">{name}\n{seq}\n")
    hmm = tmp_path / f"{name}.hmm"
    subprocess.run(["hmmbuild", "-n", name, "--amino", str(hmm), str(fasta)], check=True, capture_output=True)
    return hmm.read_text()


def test_run_produces_combined_summary(tmp_path):
    # Database: one scored PPI template (T1-T2) and one scored DDI template (DOMA-DOMB).
    ppi = tmp_path / "ppis.tsv"
    ppi.write_text("protein_a\tprotein_b\ttaxid\tscore\nT1\tT2\t9606\t0.8\n")
    templates = tmp_path / "human.fasta"
    templates.write_text(f">sp|T1|X_HUMAN\n{SEQ_1}\n>sp|T2|Y_HUMAN\n{SEQ_2}\n")
    ddis = tmp_path / "ddis.tsv"
    ddis.write_text("pfam_a\tpfam_b\tscore\nDOMA\tDOMB\t0.7\n")
    pfam = tmp_path / "toy_pfam.hmm"
    pfam.write_text(build_hmm(tmp_path, "DOMA", SEQ_1) + build_hmm(tmp_path, "DOMB", SEQ_2))

    db_dir = tmp_path / "db"
    build_database(db_dir, ppi, [f"9606={templates}"], logger, ddi_plain=ddis, pfam_hmm=pfam)

    query = tmp_path / "query.fasta"
    query.write_text(f">X\n{SEQ_1}\n>Y\n{SEQ_2}\n")
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("X\tY\n")

    workdir = tmp_path / "run"
    result = runner.invoke(
        app,
        ["run", "--db", str(db_dir), "--workdir", str(workdir),
         "--fasta", str(query), "--pairs", str(pairs), "--fused"],
    )
    assert result.exit_code == 0, result.output

    combined = pd.read_csv(workdir / "results" / "combined.summary.tsv", sep="\t")
    row = combined.iloc[0]
    assert (row.query_a, row.query_b) == ("X", "Y")
    assert row.s_im == pytest.approx(0.8, abs=1e-4)
    assert row.s_ddi == pytest.approx(0.7, abs=1e-4)
    assert row.s_fused == pytest.approx(1 - 0.2 * 0.3, abs=1e-4)

    # Second invocation must reuse both cached stages (resume behavior).
    result2 = runner.invoke(
        app,
        ["run", "--db", str(db_dir), "--workdir", str(workdir),
         "--fasta", str(query), "--pairs", str(pairs), "--fused"],
    )
    assert result2.exit_code == 0, result2.output
