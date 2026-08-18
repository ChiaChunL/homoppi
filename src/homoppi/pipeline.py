"""Reusable stage runners shared by the individual commands and `homoppi run`."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from .config import BlastParams, HmmscanParams, snapshot_params
from .workdir import Workdir

IM_SUMMARY = "interolog.summary.tsv"
IM_EVIDENCE = "interolog.evidence.tsv"
DDI_SUMMARY = "ddi.summary.tsv"
DDI_EVIDENCE = "ddi.evidence.tsv"
COMBINED_SUMMARY = "combined.summary.tsv"


def interolog_stage(
    wd: Workdir,
    db_dir: Path,
    logger: logging.Logger,
    *,
    pairs_path: Path | None,
    fasta: Path | None,
    taxid_set: set[int] | None,
    include_self: bool,
    default_score: float,
    min_score: float,
    no_evidence: bool,
    blast_params: BlastParams,
    blastp_bin: Path | None = None,
    force: bool = False,
) -> dict:
    from .blast import ensure_blast_stage
    from .interolog import infer_all_mode, infer_pairs_mode, load_homologs
    from .pairs import read_pairs
    from .ppidb import PPIIndex

    homologs_path = ensure_blast_stage(wd, db_dir, fasta, blast_params, logger, force=force, blastp_bin=blastp_bin)

    logger.info("Loading PPI template index from %s", db_dir)
    index = PPIIndex.load(db_dir, taxids=taxid_set)
    logger.info("Template index: %s pairs, %s proteins", f"{len(index.pairs):,}", f"{len(index.adj):,}")
    homologs = load_homologs(homologs_path, taxids=taxid_set)
    logger.info("Homolog table: %s query proteins with hits", f"{len(homologs):,}")

    summary_path = wd.results_dir / IM_SUMMARY
    evidence_path = None if no_evidence else wd.results_dir / IM_EVIDENCE

    if pairs_path is not None:
        query_pairs = read_pairs(pairs_path, logger)
        stats = infer_pairs_mode(
            query_pairs, homologs, index, summary_path, evidence_path, include_self, default_score, logger
        )
    else:
        logger.info("No --pairs given: running proteome-wide all-vs-all mode (interolog).")
        stats = infer_all_mode(
            homologs, index, summary_path, evidence_path, include_self, default_score, min_score, logger
        )

    snapshot_params(
        wd.results_dir / "interolog.params.json",
        {
            "mode": "pairs" if pairs_path is not None else "all-vs-all",
            "db": str(db_dir.resolve()),
            "pairs": str(pairs_path.resolve()) if pairs_path else None,
            "taxids": sorted(taxid_set) if taxid_set else None,
            "include_self": include_self,
            "default_template_score": default_score,
            "min_score": min_score,
            "blast": blast_params.to_dict(),
            "stats": stats,
        },
    )
    return stats


def ddi_stage(
    wd: Workdir,
    db_dir: Path,
    logger: logging.Logger,
    *,
    pairs_path: Path | None,
    fasta: Path | None,
    include_self: bool,
    default_score: float,
    min_score: float,
    no_evidence: bool,
    hmm_params: HmmscanParams,
    hmmscan_bin: Path | None = None,
    force: bool = False,
) -> dict:
    from .ddi import infer_all_mode, infer_pairs_mode, load_domains
    from .ddidb import DDIIndex
    from .hmmer import ensure_hmmscan_stage
    from .pairs import read_pairs

    domains_path = ensure_hmmscan_stage(wd, db_dir, fasta, hmm_params, logger, force=force, hmmscan_bin=hmmscan_bin)

    logger.info("Loading DDI template index from %s", db_dir)
    index = DDIIndex.load(db_dir)
    logger.info("DDI index: %s pairs, %s domains", f"{len(index.pairs):,}", f"{len(index.adj):,}")
    domains = load_domains(domains_path)
    logger.info("Domain table: %s query proteins with domains", f"{len(domains):,}")

    summary_path = wd.results_dir / DDI_SUMMARY
    evidence_path = None if no_evidence else wd.results_dir / DDI_EVIDENCE

    if pairs_path is not None:
        query_pairs = read_pairs(pairs_path, logger)
        stats = infer_pairs_mode(
            query_pairs, domains, index, summary_path, evidence_path, include_self, default_score, logger
        )
    else:
        logger.info("No --pairs given: running proteome-wide all-vs-all mode (DDI).")
        stats = infer_all_mode(
            domains, index, summary_path, evidence_path, include_self, default_score, min_score, logger
        )

    snapshot_params(
        wd.results_dir / "ddi.params.json",
        {
            "mode": "pairs" if pairs_path is not None else "all-vs-all",
            "db": str(db_dir.resolve()),
            "pairs": str(pairs_path.resolve()) if pairs_path else None,
            "include_self": include_self,
            "default_template_score": default_score,
            "min_score": min_score,
            "hmmscan": hmm_params.to_dict(),
            "stats": stats,
        },
    )
    return stats


def merge_stage(wd: Workdir, logger: logging.Logger, *, fused: bool = False) -> Path:
    """Join the IM and DDI summaries into one table; pairs missing from one side get 0."""
    im_path = wd.results_dir / IM_SUMMARY
    ddi_path = wd.results_dir / DDI_SUMMARY
    for path, name in ((im_path, "interolog"), (ddi_path, "ddi")):
        if not path.exists():
            raise FileNotFoundError(f"cannot merge: {name} summary missing at {path}")

    im = pd.read_csv(im_path, sep="\t", dtype={"query_a": str, "query_b": str})
    ddi = pd.read_csv(ddi_path, sep="\t", dtype={"query_a": str, "query_b": str})
    im = im[["query_a", "query_b", "n_templates", "s_im"]].rename(columns={"n_templates": "n_im_templates"})
    ddi = ddi[["query_a", "query_b", "n_templates", "s_ddi"]].rename(columns={"n_templates": "n_ddi_templates"})

    merged = im.merge(ddi, on=["query_a", "query_b"], how="outer")
    for col in ("n_im_templates", "n_ddi_templates"):
        merged[col] = merged[col].fillna(0).astype(int)
    for col in ("s_im", "s_ddi"):
        merged[col] = merged[col].fillna(0.0)
    if fused:
        merged["s_fused"] = 1.0 - (1.0 - merged["s_im"]) * (1.0 - merged["s_ddi"])

    out_path = wd.results_dir / COMBINED_SUMMARY
    merged.to_csv(out_path, sep="\t", index=False, float_format="%.4f")
    logger.info("Combined summary (%s pairs) written to %s", f"{len(merged):,}", out_path)
    return out_path
