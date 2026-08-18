import logging

import pytest

from homoppi.pairs import read_pairs

logger = logging.getLogger("test")


def test_read_pairs_with_header(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text("protein_a\tprotein_b\nA\tB\nC\tD\n")
    assert read_pairs(path, logger) == [("A", "B"), ("C", "D")]


def test_read_pairs_without_header(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text("A\tB\nC\tD\n")
    assert read_pairs(path, logger) == [("A", "B"), ("C", "D")]


def test_duplicates_collapsed_either_orientation(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text("A\tB\nB\tA\nA\tB\nC\tD\n")
    assert read_pairs(path, logger) == [("A", "B"), ("C", "D")]


def test_extra_columns_ignored(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text("A\tB\textra\tstuff\n")
    assert read_pairs(path, logger) == [("A", "B")]


def test_empty_id_raises(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text("A\t\n")
    with pytest.raises(ValueError, match="empty protein ID"):
        read_pairs(path, logger)
