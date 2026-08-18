"""blastp execution and hit filtering (the evidence-preparation step for IM)."""

from __future__ import annotations

import logging
import math
from pathlib import Path

import pandas as pd

from .config import BlastParams
from .external import find_binary, run_command
from .fasta import count_sequences
from .workdir import Workdir

FMT6_FIELDS = ["qseqid", "sseqid", "pident", "length", "qcovs", "evalue", "bitscore", "slen", "sstart", "send"]
LEGACY_FMT6_FIELDS = FMT6_FIELDS[:7]  # raw files written before subject columns were added
OUTFMT = "6 " + " ".join(FMT6_FIELDS)

RAW_TSV = "raw.fmt6.tsv"
HOMOLOGS_TSV = "homologs.tsv"
PROTEIN2TAXID_TSV = "protein2taxid.tsv"


def blastdb_prefix(db_dir: Path) -> Path:
    return db_dir / "blastdb" / "templates"


def check_blast_component(db_dir: Path) -> None:
    prefix = blastdb_prefix(db_dir)
    if not prefix.with_suffix(".pin").exists():
        raise FileNotFoundError(f"no BLAST component in database {db_dir}; run `homoppi makedb` first.")
    if not (db_dir / "blastdb" / PROTEIN2TAXID_TSV).exists():
        raise FileNotFoundError(f"{PROTEIN2TAXID_TSV} missing in {db_dir}/blastdb; rebuild with `homoppi makedb`.")


def load_protein2taxid(db_dir: Path) -> dict[str, int]:
    path = db_dir / "blastdb" / PROTEIN2TAXID_TSV
    mapping: dict[str, int] = {}
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            protein, taxid = line.rstrip("\n").split("\t")
            mapping[protein] = int(taxid)
    return mapping


def run_blastp(
    query_fasta: Path,
    db_dir: Path,
    out_path: Path,
    params: BlastParams,
    log_path: Path,
    logger: logging.Logger,
    blastp_bin: Path | None = None,
) -> None:
    blastp = find_binary("blastp", blastp_bin)
    cmd = [
        blastp,
        "-query", str(query_fasta),
        "-db", str(blastdb_prefix(db_dir)),
        "-out", str(out_path),
        "-outfmt", OUTFMT,
        "-evalue", str(params.evalue),
        "-max_target_seqs", str(params.max_target_seqs),
        "-num_threads", str(params.threads),
    ]
    logger.info("blastp does not report progress; follow along with: tail -f %s", log_path)
    run_command(cmd, log_path, logger)


def _read_raw(raw_path: Path) -> pd.DataFrame:
    """Read a raw fmt6 file, tolerating the legacy 7-column layout."""
    df = pd.read_csv(raw_path, sep="\t", header=None, dtype={0: str, 1: str})
    if df.shape[1] == len(FMT6_FIELDS):
        df.columns = FMT6_FIELDS
    elif df.shape[1] == len(LEGACY_FMT6_FIELDS):
        df.columns = LEGACY_FMT6_FIELDS
        df["slen"] = math.nan
        df["sstart"] = math.nan
        df["send"] = math.nan
    else:
        raise ValueError(f"unexpected column count ({df.shape[1]}) in {raw_path}")
    return df


def raw_has_subject_columns(raw_path: Path) -> bool:
    with open(raw_path) as fh:
        first = fh.readline()
    return len(first.rstrip("\n").split("\t")) == len(FMT6_FIELDS)


def _subject_coverage(df: pd.DataFrame) -> pd.Series:
    """Union subject coverage (%) per (query, subject), mirroring qcovs semantics.

    Merges the subject intervals of all HSPs of a pair; computed on the raw
    hits so weaker secondary HSPs still contribute, as blastp does for qcovs.
    """
    sub = df[["qseqid", "sseqid", "slen", "sstart", "send"]].sort_values(
        ["qseqid", "sseqid", "sstart"], kind="stable"
    )
    coverage: dict[tuple[str, str], float] = {}
    key = None
    covered = 0
    current_start = current_end = 0
    slen = 1.0

    def flush() -> None:
        if key is not None:
            coverage[key] = 100.0 * (covered + current_end - current_start + 1) / slen

    for row in sub.itertuples(index=False):
        row_key = (row.qseqid, row.sseqid)
        if row_key != key:
            flush()
            key, slen = row_key, float(row.slen)
            covered, current_start, current_end = 0, int(row.sstart), int(row.send)
            continue
        if row.sstart > current_end + 1:
            covered += current_end - current_start + 1
            current_start, current_end = int(row.sstart), int(row.send)
        else:
            current_end = max(current_end, int(row.send))
    flush()
    return df.apply(lambda r: coverage[(r.qseqid, r.sseqid)], axis=1)


def filter_hits(
    raw_path: Path,
    protein2taxid: dict[str, int],
    params: BlastParams,
    n_query_sequences: int,
    logger: logging.Logger,
) -> pd.DataFrame:
    """Filter raw fmt6 hits by the IM thresholds; keep the best hit per (query, template)."""
    df = _read_raw(raw_path)
    n_raw = len(df)

    has_subject_cols = df["slen"].notna().any()
    if has_subject_cols:
        df["scov"] = _subject_coverage(df)
    else:
        df["scov"] = math.nan
        if params.min_subject_coverage > 0:
            raise ValueError(
                f"{raw_path} predates subject-coverage columns; rerun the blast stage with --force "
                "to use --min-subject-coverage."
            )

    df = df[
        (df["pident"] >= params.min_identity)
        & (df["qcovs"] >= params.min_coverage)
        & (df["evalue"] <= params.evalue)
    ]
    if params.min_subject_coverage > 0:
        df = df[df["scov"] >= params.min_subject_coverage]
    df = df.sort_values("evalue", kind="stable").drop_duplicates(subset=["qseqid", "sseqid"])
    df["taxid"] = df["sseqid"].map(protein2taxid)
    n_unmapped = int(df["taxid"].isna().sum())
    if n_unmapped:
        logger.warning("%s hits target proteins missing from protein2taxid map; dropped.", f"{n_unmapped:,}")
        df = df.dropna(subset=["taxid"])
    df["taxid"] = df["taxid"].astype(int)

    df = df.rename(columns={"qseqid": "query_id", "sseqid": "template_id", "qcovs": "qcov"})
    df = df[["query_id", "template_id", "taxid", "pident", "qcov", "scov", "evalue", "bitscore"]]

    scov_note = f", subject-cov>={params.min_subject_coverage}%" if params.min_subject_coverage > 0 else ""
    logger.info(
        "blastp hits: %s raw -> %s after filtering (identity>=%s%%, qcov>=%s%%, E<=%s%s)",
        f"{n_raw:,}", f"{len(df):,}", params.min_identity, params.min_coverage, params.evalue, scov_note,
    )
    logger.info(
        "%s / %s query proteins have at least one template homolog",
        f"{df['query_id'].nunique():,}", f"{n_query_sequences:,}",
    )
    return df


def _raw_covers(old_search: dict, params: BlastParams, raw_path: Path) -> bool:
    """True if the recorded raw blastp output is a superset of what the new params need."""
    if old_search.get("max_target_seqs") != params.max_target_seqs:
        return False
    if old_search.get("evalue", 0.0) < params.evalue:
        return False  # looser E-value requested: raw is missing hits
    if params.min_subject_coverage > 0 and not raw_has_subject_columns(raw_path):
        return False
    return True


def ensure_blast_stage(
    wd: Workdir,
    db_dir: Path,
    query_fasta: Path | None,
    params: BlastParams,
    logger: logging.Logger,
    force: bool = False,
    blastp_bin: Path | None = None,
) -> Path:
    """Return the path to homologs.tsv, running or re-filtering the blast stage as needed.

    Tightening thresholds re-filters the cached raw output without rerunning
    blastp; loosening them (or changing the query) reruns the search.
    """
    homologs_path = wd.blast_dir / HOMOLOGS_TSV
    raw_path = wd.blast_dir / RAW_TSV

    if query_fasta is None:
        if homologs_path.exists():
            logger.info("Reusing existing homolog table %s", homologs_path)
            return homologs_path
        raise FileNotFoundError(
            f"no homolog table in {wd.path}; provide --fasta so the blast stage can run, or run `homoppi blast` first."
        )

    stage_params = {"search": params.search_params(), "filter": params.filter_params()}
    inputs = {"query_fasta": query_fasta}
    if not force and homologs_path.exists() and wd.is_stage_current("blast", stage_params, inputs):
        logger.info("blast stage is up to date; skipping (use --force to rerun).")
        return homologs_path

    record = wd.get_stage("blast")
    search_used: dict | None = None
    if (
        not force
        and raw_path.exists()
        and record is not None
        and wd.inputs_current(record, inputs)
        and _raw_covers(record.get("params", {}).get("search", {}), params, raw_path)
    ):
        search_used = record["params"]["search"]
        logger.info("Cached raw blastp output covers the new thresholds; re-filtering without rerunning blastp.")

    if search_used is None:
        check_blast_component(db_dir)
        run_blastp(query_fasta, db_dir, raw_path, params, wd.logs_dir / "blastp.log", logger, blastp_bin)
        search_used = params.search_params()

    n_query = count_sequences(query_fasta)
    hits = filter_hits(raw_path, load_protein2taxid(db_dir), params, n_query, logger)
    hits.to_csv(homologs_path, sep="\t", index=False, float_format="%.6g")
    logger.info("Homolog table written to %s", homologs_path)

    wd.mark_stage("blast", {"search": search_used, "filter": params.filter_params()}, inputs)
    return homologs_path
