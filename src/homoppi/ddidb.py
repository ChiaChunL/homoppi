"""DDI template library: 3did parsing, building (makedb) and runtime loading.

Library scores follow the published formula
    S_ddi_template = 1/2 * (S_known + S_EM)
where S_known is 1 for DDIs present in 3did (0 otherwise) and S_EM is the
EM-algorithm score in [0, 1] (0 when absent). A user-supplied pre-scored
library bypasses the formula.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

TEMPLATES_TSV = "templates.tsv"


def _canonical(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def parse_3did_flat(path: Path) -> set[tuple[str, str]]:
    """Extract canonical Pfam-accession pairs from a 3did flat file.

    Interaction records start with lines like:
        #=ID  domain1  domain2  (PF00001.24@Pfam  PF00002.10@Pfam)
    """
    pairs: set[tuple[str, str]] = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#=ID"):
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            pfam1 = fields[3].strip("()").split("@")[0].split(".")[0]
            pfam2 = fields[4].strip("()").split("@")[0].split(".")[0]
            if pfam1 and pfam2:
                pairs.add(_canonical(pfam1, pfam2))
    return pairs


def _load_em_scores(path: Path) -> dict[tuple[str, str], float]:
    df = pd.read_csv(path, sep="\t", dtype={"pfam_a": str, "pfam_b": str})
    required = {"pfam_a", "pfam_b", "em_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"EM score file {path} is missing columns: {sorted(missing)}")
    scores: dict[tuple[str, str], float] = {}
    for row in df.itertuples(index=False):
        key = _canonical(str(row.pfam_a), str(row.pfam_b))
        scores[key] = max(scores.get(key, 0.0), float(row.em_score))
    bad = [s for s in scores.values() if not 0.0 <= s <= 1.0]
    if bad:
        raise ValueError(f"{len(bad)} EM scores fall outside [0, 1] in {path}")
    return scores


def build_ddi_component(
    ddi_dir: Path,
    logger: logging.Logger,
    three_did: Path | None = None,
    em_scores: Path | None = None,
    plain: Path | None = None,
) -> dict:
    """Build <db>/ddi/templates.tsv from 3did and/or EM scores, or a pre-scored TSV."""
    if plain is not None and (three_did is not None or em_scores is not None):
        raise ValueError("--ddi (pre-scored) cannot be combined with --ddi-3did/--ddi-em.")
    if plain is None and three_did is None and em_scores is None:
        raise ValueError("DDI component needs --ddi-3did and/or --ddi-em, or a pre-scored --ddi TSV.")

    if plain is not None:
        df = pd.read_csv(plain, sep="\t", dtype={"pfam_a": str, "pfam_b": str})
        missing = {"pfam_a", "pfam_b"} - set(df.columns)
        if missing:
            raise ValueError(f"DDI library {plain} is missing columns: {sorted(missing)}")
        if "score" not in df.columns:
            logger.warning(
                "DDI library has no 'score' column; templates will rely on --default-template-score."
            )
            df["score"] = math.nan
        df["score"] = pd.to_numeric(df["score"], errors="coerce")
        out_of_range = df["score"].dropna().pipe(lambda s: ((s < 0) | (s > 1)).sum())
        if out_of_range:
            raise ValueError(f"{out_of_range} DDI scores fall outside [0, 1].")
        swap = df["pfam_a"] > df["pfam_b"]
        df.loc[swap, ["pfam_a", "pfam_b"]] = df.loc[swap, ["pfam_b", "pfam_a"]].values
        df = df.groupby(["pfam_a", "pfam_b"], as_index=False).agg(score=("score", "max"))
        df["source"] = "custom"
    else:
        known = parse_3did_flat(three_did) if three_did is not None else set()
        em = _load_em_scores(em_scores) if em_scores is not None else {}
        if three_did is not None:
            logger.info("3did: %s canonical DDI pairs parsed from %s", f"{len(known):,}", three_did)
        if em_scores is not None:
            logger.info("EM: %s scored DDI pairs loaded from %s", f"{len(em):,}", em_scores)

        rows = []
        for pair in sorted(known | set(em)):
            s_known = 1.0 if pair in known else 0.0
            s_em = em.get(pair, 0.0)
            source = "3did+em" if pair in known and pair in em else ("3did" if pair in known else "em")
            rows.append(
                {"pfam_a": pair[0], "pfam_b": pair[1], "source": source, "score": 0.5 * (s_known + s_em)}
            )
        df = pd.DataFrame(rows, columns=["pfam_a", "pfam_b", "source", "score"])

    n_self = int((df["pfam_a"] == df["pfam_b"]).sum())
    ddi_dir.mkdir(parents=True, exist_ok=True)
    df = df[["pfam_a", "pfam_b", "source", "score"]]
    df.to_csv(ddi_dir / TEMPLATES_TSV, sep="\t", index=False, float_format="%.6g")

    stats = {
        "template_pairs": len(df),
        "homotypic_templates": n_self,
        "by_source": df["source"].value_counts().to_dict(),
    }
    logger.info(
        "DDI library: %s template pairs (%s homotypic) by source: %s",
        f"{len(df):,}", n_self, stats["by_source"],
    )
    return stats


@dataclass
class DDIIndex:
    """In-memory index of the DDI library; same shape the inference engine expects."""

    pairs: dict[tuple[str, str], tuple[float, str]]  # (a, b) sorted -> (score or nan, source)
    adj: dict[str, set[str]]

    @classmethod
    def load(cls, db_dir: Path) -> DDIIndex:
        path = db_dir / "ddi" / TEMPLATES_TSV
        if not path.exists():
            raise FileNotFoundError(
                f"no DDI component in database {db_dir}; run `homoppi makedb` with DDI inputs first."
            )
        df = pd.read_csv(path, sep="\t", dtype={"pfam_a": str, "pfam_b": str, "source": str})

        pairs: dict[tuple[str, str], tuple[float, str]] = {}
        adj: dict[str, set[str]] = {}
        for a, b, source, score in df.itertuples(index=False):
            score = float(score) if not pd.isna(score) else math.nan
            pairs[(a, b)] = (score, source)
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
        return cls(pairs=pairs, adj=adj)
