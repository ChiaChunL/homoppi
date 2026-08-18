import re

from homoppi.fasta import count_sequences, extract_id, read_fasta, write_fasta


def test_read_fasta_multiline(tmp_path):
    path = tmp_path / "test.fasta"
    path.write_text(">seq1 description here\nMKV\nLLA\n>seq2\nAAA\n")
    records = list(read_fasta(path))
    assert records == [("seq1 description here", "MKVLLA"), ("seq2", "AAA")]


def test_extract_id_uniprot_style():
    assert extract_id("sp|P12345|NAME_HUMAN some description") == "P12345"
    assert extract_id("tr|A0A087WXM9|A0A087WXM9_HUMAN") == "A0A087WXM9"


def test_extract_id_plain_first_token():
    assert extract_id("ENSANP00000012345 pep chromosome:1") == "ENSANP00000012345"


def test_extract_id_custom_regex_wins():
    regex = re.compile(r"GN=(\w+)")
    assert extract_id("sp|P12345|NAME_HUMAN GN=TP53 PE=1", regex) == "TP53"


def test_write_and_count_roundtrip(tmp_path):
    path = tmp_path / "out.fasta"
    n = write_fasta([("a", "MKV" * 50), ("b", "AAA")], path)
    assert n == 2
    assert count_sequences(path) == 2
    assert dict(read_fasta(path)) == {"a": "MKV" * 50, "b": "AAA"}
