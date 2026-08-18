"""Interolog-mapping inference: IM-specific wiring around the generic engine.

A query pair (A, B) is supported by template pair (A', B') when A has a
filtered blastp hit to A', B to B', and (A', B') is in the PPI template
library.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from . import inference
from .inference import MethodSpec
from .ppidb import PPIIndex


@dataclass(frozen=True)
class Hit:
    pident: float
    qcov: float
    evalue: float
    scov: float = math.nan  # subject (template) coverage; NaN for legacy homolog tables


def _fmt(value: float, spec: str = ".4g") -> str:
    return "" if math.isnan(value) else format(value, spec)


SPEC = MethodSpec(
    label="interolog",
    score_column="s_im",
    meta_columns=("taxid",),
    hit_stat_columns=("pident", "qcov", "scov", "evalue"),
    hit_stats=lambda h: (_fmt(h.pident), _fmt(h.qcov), _fmt(h.scov), _fmt(h.evalue, ".3g")),
)


def load_homologs(path: Path, taxids: set[int] | None = None) -> dict[str, dict[str, Hit]]:
    """homologs.tsv -> {query_id: {template_id: Hit}}."""
    df = pd.read_csv(path, sep="\t", dtype={"query_id": str, "template_id": str})
    if "scov" not in df.columns:  # legacy homolog table
        df["scov"] = math.nan
    if taxids is not None:
        df = df[df["taxid"].isin(taxids)]
    homologs: dict[str, dict[str, Hit]] = {}
    for row in df.itertuples(index=False):
        homologs.setdefault(row.query_id, {})[row.template_id] = Hit(
            pident=float(row.pident), qcov=float(row.qcov), evalue=float(row.evalue), scov=float(row.scov)
        )
    return homologs


def infer_pairs_mode(
    pairs: list[tuple[str, str]],
    homologs: dict[str, dict[str, Hit]],
    index: PPIIndex,
    summary_path: Path,
    evidence_path: Path | None,
    include_self: bool,
    default_score: float,
    logger: logging.Logger,
) -> dict:
    return inference.infer_pairs_mode(
        pairs, homologs, index, SPEC, summary_path, evidence_path,
        include_self=include_self, default_score=default_score, logger=logger,
    )


def infer_all_mode(
    homologs: dict[str, dict[str, Hit]],
    index: PPIIndex,
    summary_path: Path,
    evidence_path: Path | None,
    include_self: bool,
    default_score: float,
    min_score: float,
    logger: logging.Logger,
) -> dict:
    return inference.infer_all_mode(
        homologs, index, SPEC, summary_path, evidence_path,
        include_self=include_self, default_score=default_score, min_score=min_score, logger=logger,
    )
