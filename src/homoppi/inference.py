"""Generic template-mapping inference engine shared by the IM and DDI methods.

Both methods share one structure: each query protein maps to a set of
*features* (homologous template proteins for IM, Pfam domains for DDI); a
query pair (A, B) is supported by a template pair (fa, fb) from a library
when A carries fa, B carries fb, and (fa, fb) is a library pair. Each
distinct template pair counts once per query pair, and evidence is combined
with Bayesian integration.

A template index exposes `pairs` (canonical key -> tuple whose first element
is the score, possibly NaN, and whose remaining elements are method-specific
metadata) and `adj` (feature -> partner features).
"""

from __future__ import annotations

import logging
import math
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, TextIO

from tqdm import tqdm

from .scoring import bayes_integration


class TemplateIndex(Protocol):
    pairs: Mapping[tuple[str, str], tuple]
    adj: Mapping[str, set[str]]


@dataclass(frozen=True)
class MethodSpec:
    """Method-specific naming and evidence formatting."""

    label: str  # for log messages, e.g. "interolog"
    score_column: str  # e.g. "s_im"
    meta_columns: tuple[str, ...]  # names of index metadata fields, e.g. ("taxid",)
    hit_stat_columns: tuple[str, ...]  # per-side hit statistics, e.g. ("pident", "qcov", "evalue")
    hit_stats: Callable[[Any], tuple[str, ...]]  # format a hit into hit_stat_columns values

    @property
    def summary_columns(self) -> list[str]:
        return [
            "query_a", "query_b", "n_templates", self.score_column,
            "best_template_a", "best_template_b",
            *[f"best_template_{m}" for m in self.meta_columns],
            "best_template_score",
        ]

    @property
    def evidence_columns(self) -> list[str]:
        return [
            "query_a", "query_b", "template_a", "template_b",
            *self.meta_columns, "template_score",
            *[f"{c}_a" for c in self.hit_stat_columns],
            *[f"{c}_b" for c in self.hit_stat_columns],
        ]


def _fmt_score(raw: float) -> str:
    return "" if math.isnan(raw) else f"{raw:.4g}"


def _effective(raw: float, default: float) -> float:
    return default if math.isnan(raw) else raw


def _evidence_row(
    spec: MethodSpec, qa: str, qb: str, ta: str, tb: str, value: tuple, hit_a: Any, hit_b: Any
) -> str:
    raw, meta = value[0], value[1:]
    return (
        "\t".join(
            [qa, qb, ta, tb, *[str(m) for m in meta], _fmt_score(raw),
             *spec.hit_stats(hit_a), *spec.hit_stats(hit_b)]
        )
        + "\n"
    )


def infer_pairs_mode(
    pairs: list[tuple[str, str]],
    features: dict[str, dict[str, Any]],
    index: TemplateIndex,
    spec: MethodSpec,
    summary_path: Path,
    evidence_path: Path | None,
    *,
    include_self: bool,
    default_score: float,
    logger: logging.Logger,
) -> dict:
    """Score an explicit list of query pairs; every non-skipped pair gets a summary row."""
    n_self_skipped = 0
    n_supported = 0

    evidence_fh: TextIO | None = None
    with open(summary_path, "w") as summary_fh:
        summary_fh.write("\t".join(spec.summary_columns) + "\n")
        if evidence_path is not None:
            evidence_fh = open(evidence_path, "w")
            evidence_fh.write("\t".join(spec.evidence_columns) + "\n")
        try:
            for qa, qb in tqdm(pairs, desc=f"Scoring pairs ({spec.label})", unit="pair"):
                if qa == qb and not include_self:
                    n_self_skipped += 1
                    continue
                feats_a = features.get(qa, {})
                feats_b = features.get(qb, {})

                # canonical template pair -> (feature matched to A, feature matched to B)
                found: dict[tuple[str, str], tuple[str, str]] = {}
                if feats_a and feats_b:
                    set_b = set(feats_b)
                    for ta in feats_a:
                        partners = index.adj.get(ta)
                        if not partners:
                            continue
                        for tb in partners & set_b:
                            key = (ta, tb) if ta <= tb else (tb, ta)
                            found.setdefault(key, (ta, tb))

                scores = [_effective(index.pairs[key][0], default_score) for key in found]
                s_total = bayes_integration(scores)
                best_fields = [""] * (3 + len(spec.meta_columns))
                if found:
                    n_supported += 1
                    best_key = max(found, key=lambda k: _effective(index.pairs[k][0], default_score))
                    value = index.pairs[best_key]
                    raw, meta = value[0], value[1:]
                    best_fields = [
                        best_key[0], best_key[1], *[str(m) for m in meta],
                        _fmt_score(raw) or f"{default_score:.4g}",
                    ]

                summary_fh.write(
                    "\t".join([qa, qb, str(len(found)), f"{s_total:.4f}", *best_fields]) + "\n"
                )
                if evidence_fh is not None:
                    for key, (ta, tb) in found.items():
                        evidence_fh.write(
                            _evidence_row(spec, qa, qb, ta, tb, index.pairs[key], feats_a[ta], feats_b[tb])
                        )
        finally:
            if evidence_fh is not None:
                evidence_fh.close()

    if n_self_skipped:
        logger.warning(
            "Skipped %s self pairs (A == B); use --include-self to score them.", f"{n_self_skipped:,}"
        )
    stats = {
        "pairs_scored": len(pairs) - n_self_skipped,
        "pairs_with_evidence": n_supported,
        "self_pairs_skipped": n_self_skipped,
    }
    logger.info(
        "%s / %s pairs have %s evidence; summary written to %s",
        f"{n_supported:,}", f"{stats['pairs_scored']:,}", spec.label, summary_path,
    )
    return stats


def infer_all_mode(
    features: dict[str, dict[str, Any]],
    index: TemplateIndex,
    spec: MethodSpec,
    summary_path: Path,
    evidence_path: Path | None,
    *,
    include_self: bool,
    default_score: float,
    min_score: float,
    logger: logging.Logger,
) -> dict:
    """Proteome-wide enumeration: expand every library template pair to query pairs.

    Iterates from the template side and streams evidence to disk; the summary
    keeps one running product per query pair, so memory scales with the number
    of predicted pairs rather than with the evidence count.
    """
    reverse: dict[str, list[str]] = {}
    for query, feats in features.items():
        for feature in feats:
            reverse.setdefault(feature, []).append(query)

    # query pair -> [product(1-s), n_templates, best_eff, best_ta, best_tb, best_value]
    acc: dict[tuple[str, str], list] = {}
    evidence_fh: TextIO | None = None
    if evidence_path is not None:
        evidence_fh = open(evidence_path, "w")
        evidence_fh.write("\t".join(spec.evidence_columns) + "\n")

    try:
        for (ta, tb), value in tqdm(index.pairs.items(), desc=f"Expanding templates ({spec.label})", unit="tmpl"):
            queries_a = reverse.get(ta)
            queries_b = reverse.get(tb)
            if not queries_a or not queries_b:
                continue
            eff = _effective(value[0], default_score)

            emitted: set[tuple[str, str]] = set()
            for qa in queries_a:
                for qb in queries_b:
                    if qa == qb and not include_self:
                        continue
                    key = (qa, qb) if qa <= qb else (qb, qa)
                    if key in emitted:
                        continue
                    emitted.add(key)

                    entry = acc.get(key)
                    if entry is None:
                        entry = acc[key] = [1.0, 0, -1.0, "", "", None]
                    entry[0] *= 1.0 - eff
                    entry[1] += 1
                    if eff > entry[2]:
                        entry[2], entry[3], entry[4], entry[5] = eff, ta, tb, value

                    if evidence_fh is not None:
                        # align templates to the canonical query orientation
                        out_ta, out_tb = (ta, tb) if key == (qa, qb) else (tb, ta)
                        evidence_fh.write(
                            _evidence_row(
                                spec, key[0], key[1], out_ta, out_tb, value,
                                features[key[0]][out_ta], features[key[1]][out_tb],
                            )
                        )
    finally:
        if evidence_fh is not None:
            evidence_fh.close()

    n_written = 0
    with open(summary_path, "w") as summary_fh:
        summary_fh.write("\t".join(spec.summary_columns) + "\n")
        for (qa, qb), (product, n, best_eff, bta, btb, bvalue) in tqdm(
            acc.items(), desc="Writing summary", unit="pair"
        ):
            s_total = 1.0 - product
            if s_total < min_score:
                continue
            n_written += 1
            raw, meta = bvalue[0], bvalue[1:]
            best_score = _fmt_score(raw) or f"{best_eff:.4g}"
            summary_fh.write(
                "\t".join(
                    [qa, qb, str(n), f"{s_total:.4f}", bta, btb, *[str(m) for m in meta], best_score]
                )
                + "\n"
            )

    logger.info(
        "Proteome-wide %s mode predicted %s query pairs (%s written with %s>=%s); summary at %s",
        spec.label, f"{len(acc):,}", f"{n_written:,}", spec.score_column, min_score, summary_path,
    )
    return {"pairs_predicted": len(acc), "pairs_written": n_written}
