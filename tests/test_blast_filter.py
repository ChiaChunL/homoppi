import logging
import math

import pytest

from homoppi.blast import FMT6_FIELDS, LEGACY_FMT6_FIELDS, filter_hits
from homoppi.config import BlastParams

logger = logging.getLogger("test")


def write_raw(tmp_path, rows, n_fields=None):
    n_fields = n_fields if n_fields is not None else len(FMT6_FIELDS)
    path = tmp_path / "raw.fmt6.tsv"
    lines = ["\t".join(str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n")
    assert all(len(row) == n_fields for row in rows)
    return path


P2TAXID = {"P1": 9606, "P2": 9606, "Q1": 10090}
PARAMS = BlastParams(evalue=1e-10, min_identity=30.0, min_coverage=40.0)


def test_thresholds_are_inclusive(tmp_path):
    raw = write_raw(
        tmp_path,
        [
            # qseqid sseqid pident length qcovs evalue bitscore slen sstart send
            ["A", "P1", 30.0, 100, 40.0, 1e-10, 200.0, 200, 1, 100],  # exactly on thresholds -> kept
            ["A", "P2", 29.9, 100, 40.0, 1e-10, 200.0, 200, 1, 100],  # identity below -> dropped
            ["A", "Q1", 30.0, 100, 39.9, 1e-10, 200.0, 200, 1, 100],  # coverage below -> dropped
            ["B", "P1", 90.0, 100, 90.0, 1e-05, 200.0, 200, 1, 100],  # evalue above -> dropped
        ],
    )
    df = filter_hits(raw, P2TAXID, PARAMS, n_query_sequences=2, logger=logger)
    assert df[["query_id", "template_id"]].values.tolist() == [["A", "P1"]]


def test_best_hit_kept_per_query_template(tmp_path):
    raw = write_raw(
        tmp_path,
        [
            ["A", "P1", 50.0, 100, 80.0, 1e-20, 100.0, 200, 1, 100],
            ["A", "P1", 95.0, 100, 90.0, 1e-50, 300.0, 200, 1, 100],  # better evalue -> kept
        ],
    )
    df = filter_hits(raw, P2TAXID, PARAMS, n_query_sequences=1, logger=logger)
    assert len(df) == 1
    assert df.iloc[0].pident == 95.0


def test_subject_coverage_is_hsp_union(tmp_path):
    raw = write_raw(
        tmp_path,
        [
            # two HSPs on subject intervals 1-50 and 41-100 of slen 200 -> union 100 -> 50%
            ["A", "P1", 90.0, 60, 90.0, 1e-50, 300.0, 200, 1, 50],
            ["A", "P1", 85.0, 60, 90.0, 1e-40, 250.0, 200, 41, 100],
        ],
    )
    df = filter_hits(raw, P2TAXID, PARAMS, n_query_sequences=1, logger=logger)
    assert df.iloc[0].scov == pytest.approx(50.0)

    strict = BlastParams(evalue=1e-10, min_subject_coverage=60.0)
    df = filter_hits(raw, P2TAXID, strict, n_query_sequences=1, logger=logger)
    assert df.empty

    loose = BlastParams(evalue=1e-10, min_subject_coverage=40.0)
    df = filter_hits(raw, P2TAXID, loose, n_query_sequences=1, logger=logger)
    assert len(df) == 1


def test_unmapped_template_dropped(tmp_path):
    raw = write_raw(tmp_path, [["A", "UNKNOWN", 90.0, 100, 90.0, 1e-50, 300.0, 200, 1, 100]])
    df = filter_hits(raw, P2TAXID, PARAMS, n_query_sequences=1, logger=logger)
    assert df.empty


def test_legacy_raw_without_subject_columns(tmp_path):
    raw = write_raw(
        tmp_path,
        [["A", "P1", 90.0, 100, 90.0, 1e-50, 300.0]],
        n_fields=len(LEGACY_FMT6_FIELDS),
    )
    df = filter_hits(raw, P2TAXID, PARAMS, n_query_sequences=1, logger=logger)
    assert len(df) == 1
    assert math.isnan(df.iloc[0].scov)

    with pytest.raises(ValueError, match="predates subject-coverage"):
        filter_hits(raw, P2TAXID, BlastParams(min_subject_coverage=50.0), 1, logger)
