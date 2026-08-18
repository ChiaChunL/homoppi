"""hmmscan execution and domtblout filtering (the evidence-preparation step for DDI)."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pandas as pd

from .config import HmmscanParams
from .external import find_binary, run_command
from .fasta import count_sequences
from .workdir import Workdir

RAW_DOMTBLOUT = "raw.domtblout"
DOMAINS_TSV = "domains.tsv"
PFAM_HMM = "pfam.hmm"
CLANS_TSV = "clans.tsv"

DOMAINS_COLUMNS = ["query_id", "pfam_acc", "domain_name", "cevalue", "hmm_cov", "ali_from", "ali_to"]


def pfam_path(db_dir: Path) -> Path:
    return db_dir / "pfam" / PFAM_HMM


def clans_path(db_dir: Path) -> Path:
    return db_dir / "pfam" / CLANS_TSV


def check_pfam_component(db_dir: Path) -> None:
    hmm = pfam_path(db_dir)
    if not hmm.exists():
        raise FileNotFoundError(f"no Pfam component in database {db_dir}; run `homoppi makedb --pfam-hmm ...` first.")
    if not hmm.with_name(hmm.name + ".h3m").exists():
        raise FileNotFoundError(f"Pfam HMM in {db_dir} is not pressed; rebuild with `homoppi makedb --pfam-hmm ...`.")


def parse_pfam_dat(path: Path) -> dict[str, str]:
    """Family -> clan mapping from Pfam-A.hmm.dat (#=GF AC / #=GF CL stanzas)."""
    clans: dict[str, str] = {}
    acc: str | None = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#=GF AC"):
                acc = line.split()[2].split(".")[0]
            elif line.startswith("#=GF CL") and acc is not None:
                clans[acc] = line.split()[2]
            elif line.startswith("//"):
                acc = None
    return clans


def load_clans(db_dir: Path) -> dict[str, str]:
    path = clans_path(db_dir)
    if not path.exists():
        raise FileNotFoundError(
            f"no clan mapping in database {db_dir}; --resolve-clan-overlap needs "
            "`homoppi makedb --pfam-dat Pfam-A.hmm.dat` first."
        )
    with open(path) as fh:
        next(fh)  # header
        return dict(line.rstrip("\n").split("\t") for line in fh)


def run_hmmscan(
    query_fasta: Path,
    db_dir: Path,
    out_path: Path,
    params: HmmscanParams,
    log_path: Path,
    logger: logging.Logger,
    hmmscan_bin: Path | None = None,
) -> None:
    hmmscan = find_binary("hmmscan", hmmscan_bin)
    cmd = [hmmscan, "-o", os.devnull, "--noali", "--domtblout", str(out_path), "--cpu", str(params.threads)]
    if params.cut_tc:
        cmd.append("--cut_tc")
    else:
        cmd += ["-E", str(params.cevalue)]
    cmd += [str(pfam_path(db_dir)), str(query_fasta)]
    logger.info("hmmscan does not report progress; follow along with the domtblout file size.")
    run_command(cmd, log_path, logger)


def parse_domtblout(path: Path) -> pd.DataFrame:
    """Parse hmmscan --domtblout output into one row per domain hit."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split(maxsplit=22)
            target_name, taccession, tlen = fields[0], fields[1], int(fields[2])
            acc = taccession if taccession not in ("-", "") else target_name
            hmm_from, hmm_to = int(fields[15]), int(fields[16])
            rows.append(
                {
                    "query_id": fields[3],
                    "pfam_acc": acc.split(".")[0],  # strip Pfam version suffix
                    "domain_name": target_name,
                    "seq_evalue": float(fields[6]),
                    "cevalue": float(fields[11]),
                    "hmm_cov": round(100.0 * (hmm_to - hmm_from + 1) / tlen, 1),
                    "ali_from": int(fields[17]),
                    "ali_to": int(fields[18]),
                }
            )
    return pd.DataFrame(rows, columns=["seq_evalue", *DOMAINS_COLUMNS])


def resolve_clan_overlaps(df: pd.DataFrame, clans: dict[str, str], logger: logging.Logger) -> pd.DataFrame:
    """Drop hits overlapping a better hit from the same Pfam clan (pfam_scan-style).

    Within each query, instances are visited by ascending c-Evalue; an instance
    is dropped when it overlaps an already-accepted instance whose family
    belongs to the same clan. Families without a clan never conflict, and
    overlaps across different clans are kept (genuinely distinct annotations).
    """
    keep: list[int] = []
    n_dropped = 0
    for _, group in df.groupby("query_id", sort=False):
        accepted: list[tuple[int, int, str]] = []  # (ali_from, ali_to, clan)
        for row in group.sort_values("cevalue", kind="stable").itertuples():
            clan = clans.get(row.pfam_acc)
            if clan is not None and any(
                row.ali_from <= to and frm <= row.ali_to for frm, to, c in accepted if c == clan
            ):
                n_dropped += 1
                continue
            keep.append(row.Index)
            if clan is not None:
                accepted.append((row.ali_from, row.ali_to, clan))
    if n_dropped:
        logger.info("Clan overlap resolution dropped %s redundant domain hits.", f"{n_dropped:,}")
    return df.loc[sorted(keep)]


def filter_domains(
    df: pd.DataFrame,
    params: HmmscanParams,
    n_query_sequences: int,
    logger: logging.Logger,
    clans: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Apply the configured filters and keep the best instance per (query, domain).

    The full-sequence E-value filter mirrors hmmscan's own -E gate, so
    re-filtering a cached raw file with a tighter threshold reproduces what a
    fresh, tighter run would report.
    """
    n_raw = len(df)
    if not params.cut_tc:
        df = df[(df["seq_evalue"] <= params.cevalue) & (df["cevalue"] <= params.cevalue)]
    if params.min_hmm_coverage > 0:
        df = df[df["hmm_cov"] >= params.min_hmm_coverage]
    if params.resolve_clan_overlap:
        if clans is None:
            raise ValueError("clan overlap resolution requested but no clan mapping was provided.")
        df = resolve_clan_overlaps(df, clans, logger)
    df = df.sort_values("cevalue", kind="stable").drop_duplicates(subset=["query_id", "pfam_acc"])
    df = df.sort_values(["query_id", "ali_from"], kind="stable")[DOMAINS_COLUMNS]

    criteria = ["Pfam trusted cutoffs" if params.cut_tc else f"c-Evalue<={params.cevalue}"]
    if params.min_hmm_coverage > 0:
        criteria.append(f"hmm-cov>={params.min_hmm_coverage}%")
    if params.resolve_clan_overlap:
        criteria.append("clan overlaps resolved")
    logger.info("hmmscan hits: %s raw -> %s after filtering (%s)", f"{n_raw:,}", f"{len(df):,}", ", ".join(criteria))
    logger.info(
        "%s / %s query proteins have at least one Pfam domain",
        f"{df['query_id'].nunique():,}", f"{n_query_sequences:,}",
    )
    return df


def _raw_covers(old_search: dict, params: HmmscanParams) -> bool:
    """True if the recorded raw domtblout is a superset of what the new params need."""
    if old_search.get("cut_tc") != params.cut_tc:
        return False
    if params.cut_tc:
        return True
    return old_search.get("cevalue", 0.0) >= params.cevalue


def ensure_hmmscan_stage(
    wd: Workdir,
    db_dir: Path,
    query_fasta: Path | None,
    params: HmmscanParams,
    logger: logging.Logger,
    force: bool = False,
    hmmscan_bin: Path | None = None,
) -> Path:
    """Return the path to domains.tsv, running or re-filtering the hmmscan stage as needed."""
    domains_path = wd.hmmscan_dir / DOMAINS_TSV
    raw_path = wd.hmmscan_dir / RAW_DOMTBLOUT

    if query_fasta is None:
        if domains_path.exists():
            logger.info("Reusing existing domain table %s", domains_path)
            return domains_path
        raise FileNotFoundError(
            f"no domain table in {wd.path}; provide --fasta so the hmmscan stage can run, "
            "or run `homoppi domainanno` first."
        )

    clans = load_clans(db_dir) if params.resolve_clan_overlap else None

    stage_params = {"search": params.search_params(), "filter": params.filter_params()}
    inputs = {"query_fasta": query_fasta}
    if not force and domains_path.exists() and wd.is_stage_current("hmmscan", stage_params, inputs):
        logger.info("hmmscan stage is up to date; skipping (use --force to rerun).")
        return domains_path

    record = wd.get_stage("hmmscan")
    search_used: dict | None = None
    if (
        not force
        and raw_path.exists()
        and record is not None
        and wd.inputs_current(record, inputs)
        and _raw_covers(record.get("params", {}).get("search", {}), params)
    ):
        search_used = record["params"]["search"]
        logger.info("Cached raw hmmscan output covers the new thresholds; re-filtering without rerunning hmmscan.")

    if search_used is None:
        check_pfam_component(db_dir)
        run_hmmscan(query_fasta, db_dir, raw_path, params, wd.logs_dir / "hmmscan.log", logger, hmmscan_bin)
        search_used = params.search_params()

    n_query = count_sequences(query_fasta)
    domains = filter_domains(parse_domtblout(raw_path), params, n_query, logger, clans=clans)
    domains.to_csv(domains_path, sep="\t", index=False, float_format="%.6g")
    logger.info("Domain table written to %s", domains_path)

    wd.mark_stage("hmmscan", {"search": search_used, "filter": params.filter_params()}, inputs)
    return domains_path
