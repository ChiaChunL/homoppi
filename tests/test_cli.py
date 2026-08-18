from typer.testing import CliRunner

from homoppi import __version__
from homoppi.cli import app

runner = CliRunner()


def test_version():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_help_lists_commands():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    for command in ("makedb", "blast", "domainanno", "interolog", "ddi", "ddi-em", "run"):
        assert command in result.output


def test_run_rejects_skipping_both_methods(tmp_path):
    db = tmp_path / "db"
    db.mkdir()
    fasta = tmp_path / "q.fasta"
    fasta.write_text(">A\nMKV\n")
    result = runner.invoke(
        app,
        ["run", "--db", str(db), "--workdir", str(tmp_path / "wd"),
         "--fasta", str(fasta), "--skip-im", "--skip-ddi"],
    )
    assert result.exit_code != 0


def test_interolog_requires_valid_db(tmp_path):
    result = runner.invoke(app, ["interolog", "--db", str(tmp_path / "nope"), "--workdir", str(tmp_path / "wd")])
    assert result.exit_code != 0
