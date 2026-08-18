"""DDI inference: DDI-specific wiring around the generic engine.

A query pair (A, B) is supported by DDI template (Da, Db) when A carries
domain Da, B carries Db, and (Da, Db) is in the DDI library.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from . import inference
from .ddidb import DDIIndex
from .inference import MethodSpec


@dataclass(frozen=True)
class DomHit:
    cevalue: float


SPEC = MethodSpec(
    label="DDI",
    score_column="s_ddi",
    meta_columns=("source",),
    hit_stat_columns=("cevalue",),
    hit_stats=lambda h: (f"{h.cevalue:.3g}",),
)


def load_domains(path: Path) -> dict[str, dict[str, DomHit]]:
    """domains.tsv -> {query_id: {pfam_acc: DomHit}} (best instance per domain)."""
    df = pd.read_csv(path, sep="\t", dtype={"query_id": str, "pfam_acc": str})
    domains: dict[str, dict[str, DomHit]] = {}
    for row in df.itertuples(index=False):
        current = domains.setdefault(row.query_id, {})
        hit = current.get(row.pfam_acc)
        if hit is None or float(row.cevalue) < hit.cevalue:
            current[row.pfam_acc] = DomHit(cevalue=float(row.cevalue))
    return domains


def infer_pairs_mode(
    pairs: list[tuple[str, str]],
    domains: dict[str, dict[str, DomHit]],
    index: DDIIndex,
    summary_path: Path,
    evidence_path: Path | None,
    include_self: bool,
    default_score: float,
    logger: logging.Logger,
) -> dict:
    return inference.infer_pairs_mode(
        pairs, domains, index, SPEC, summary_path, evidence_path,
        include_self=include_self, default_score=default_score, logger=logger,
    )


def infer_all_mode(
    domains: dict[str, dict[str, DomHit]],
    index: DDIIndex,
    summary_path: Path,
    evidence_path: Path | None,
    include_self: bool,
    default_score: float,
    min_score: float,
    logger: logging.Logger,
) -> dict:
    return inference.infer_all_mode(
        domains, index, SPEC, summary_path, evidence_path,
        include_self=include_self, default_score=default_score, min_score=min_score, logger=logger,
    )
