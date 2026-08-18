"""PPI template library: building (makedb) and runtime loading."""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

TEMPLATES_TSV = "templates.tsv"
REQUIRED_COLUMNS = {"protein_a", "protein_b", "taxid"}


def build_ppi_component(ppi_tsv: Path, ppi_dir: Path, logger: logging.Logger) -> dict:
    """Validate and canonicalize a user PPI library into <db>/ppi/templates.tsv.

    Input TSV columns: protein_a, protein_b, taxid[, score]. Pairs are stored
    with sorted protein IDs; duplicates keep the maximum score. Homodimer
    templates (protein_a == protein_b) are kept: they are valid evidence.
    """
    df = pd.read_csv(ppi_tsv, sep="\t", dtype={"protein_a": str, "protein_b": str})
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"PPI library {ppi_tsv} is missing required columns: {sorted(missing)}")

    n_raw = len(df)
    if "score" not in df.columns:
        logger.warning(
            "PPI library has no 'score' column; all templates will rely on --default-template-score at prediction time."
        )
        df["score"] = np.nan
    df["score"] = pd.to_numeric(df["score"], errors="coerce")

    out_of_range = df["score"].dropna().pipe(lambda s: ((s < 0) | (s > 1)).sum())
    if out_of_range:
        raise ValueError(f"{out_of_range} scores fall outside [0, 1]; template scores must be probabilities.")

    df["protein_a"] = df["protein_a"].str.strip()
    df["protein_b"] = df["protein_b"].str.strip()
    df = df.dropna(subset=["protein_a", "protein_b", "taxid"])
    df = df[(df["protein_a"] != "") & (df["protein_b"] != "")]
    df["taxid"] = df["taxid"].astype(int)

    # Canonical order within each pair.
    swap = df["protein_a"] > df["protein_b"]
    df.loc[swap, ["protein_a", "protein_b"]] = df.loc[swap, ["protein_b", "protein_a"]].values

    df = (
        df.groupby(["protein_a", "protein_b"], as_index=False)
        .agg(taxid=("taxid", "first"), score=("score", "max"))
    )

    n_self = int((df["protein_a"] == df["protein_b"]).sum())
    n_scored = int(df["score"].notna().sum())

    ppi_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(ppi_dir / TEMPLATES_TSV, sep="\t", index=False, float_format="%.6g")

    stats = {
        "input_rows": n_raw,
        "template_pairs": len(df),
        "homodimer_templates": n_self,
        "scored_fraction": round(n_scored / len(df), 4) if len(df) else 0.0,
        "pairs_per_taxid": df["taxid"].value_counts().to_dict(),
    }
    logger.info(
        "PPI library: %s rows -> %s canonical template pairs (%s homodimers, %.1f%% scored)",
        f"{n_raw:,}", f"{len(df):,}", n_self, 100 * stats["scored_fraction"],
    )
    for taxid, count in sorted(stats["pairs_per_taxid"].items()):
        logger.info("  taxid %-8s %s pairs", taxid, f"{count:,}")
    if stats["scored_fraction"] < 0.5 and n_scored > 0:
        logger.warning(
            "Less than half of the templates carry a score; unscored ones fall back to --default-template-score."
        )
    return stats


def lib_proteins_by_taxid(ppi_dir: Path) -> dict[int, set[str]]:
    """Proteins appearing in the template library, grouped by taxid."""
    df = pd.read_csv(ppi_dir / TEMPLATES_TSV, sep="\t", dtype={"protein_a": str, "protein_b": str})
    out: dict[int, set[str]] = {}
    for taxid, sub in df.groupby("taxid"):
        out[int(taxid)] = set(sub["protein_a"]) | set(sub["protein_b"])
    return out


@dataclass
class PPIIndex:
    """In-memory index of the template library for fast pair lookup."""

    pairs: dict[tuple[str, str], tuple[float, int]]  # (a, b) sorted -> (score or nan, taxid)
    adj: dict[str, set[str]]  # protein -> interacting partners

    @classmethod
    def load(cls, db_dir: Path, taxids: set[int] | None = None) -> PPIIndex:
        path = db_dir / "ppi" / TEMPLATES_TSV
        if not path.exists():
            raise FileNotFoundError(f"no PPI component in database {db_dir}; run `homoppi makedb` first.")
        df = pd.read_csv(path, sep="\t", dtype={"protein_a": str, "protein_b": str})
        if taxids is not None:
            df = df[df["taxid"].isin(taxids)]

        pairs: dict[tuple[str, str], tuple[float, int]] = {}
        adj: dict[str, set[str]] = {}
        for a, b, taxid, score in df.itertuples(index=False):
            score = float(score) if not pd.isna(score) else math.nan
            pairs[(a, b)] = (score, int(taxid))
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
        return cls(pairs=pairs, adj=adj)

    def score_of(self, key: tuple[str, str], default: float) -> float:
        raw, _ = self.pairs[key]
        return default if math.isnan(raw) else raw
