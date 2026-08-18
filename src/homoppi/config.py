"""Default parameters and parameter snapshotting.

Defaults follow the published method: blastp hits are kept when
identity >= 30%, query coverage >= 40% and E-value <= 1e-10.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class BlastParams:
    evalue: float = 1e-10
    min_identity: float = 30.0  # percent
    min_coverage: float = 40.0  # percent of the query covered (blastp qcovs)
    min_subject_coverage: float = 0.0  # percent of the template covered; 0 disables (paper behavior)
    max_target_seqs: int = 500
    threads: int = 4

    def to_dict(self) -> dict:
        return asdict(self)

    def search_params(self) -> dict:
        """Parameters baked into the raw blastp output (a change forces a rerun)."""
        return {"evalue": self.evalue, "max_target_seqs": self.max_target_seqs}

    def filter_params(self) -> dict:
        """Pure post-filters; a change only requires re-filtering the raw output."""
        return {
            "evalue": self.evalue,
            "min_identity": self.min_identity,
            "min_coverage": self.min_coverage,
            "min_subject_coverage": self.min_subject_coverage,
        }


@dataclass(frozen=True)
class HmmscanParams:
    cevalue: float = 1e-10  # conditional (domain) E-value cutoff
    cut_tc: bool = False  # use Pfam trusted cutoffs instead of E-values
    min_hmm_coverage: float = 0.0  # percent of the HMM model matched; 0 disables (paper behavior)
    resolve_clan_overlap: bool = False  # drop overlapping hits from the same Pfam clan (keep best)
    threads: int = 4

    def to_dict(self) -> dict:
        return asdict(self)

    def search_params(self) -> dict:
        """Parameters baked into the raw domtblout (a change forces a rerun)."""
        return {"cevalue": self.cevalue, "cut_tc": self.cut_tc}

    def filter_params(self) -> dict:
        """Pure post-filters; a change only requires re-filtering the raw output."""
        return {
            "cevalue": self.cevalue,
            "cut_tc": self.cut_tc,
            "min_hmm_coverage": self.min_hmm_coverage,
            "resolve_clan_overlap": self.resolve_clan_overlap,
        }


DEFAULT_TEMPLATE_SCORE = 0.0


def parse_taxids(value: str | None) -> set[int] | None:
    """Parse a comma-separated taxid list ('9606,10090') into a set."""
    if value is None or not value.strip():
        return None
    try:
        return {int(t) for t in value.replace(" ", "").split(",") if t}
    except ValueError as exc:
        raise ValueError(f"invalid --taxids value: {value!r}") from exc


def snapshot_params(path: Path, params: dict) -> None:
    """Write the effective parameters of a run next to its outputs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(params, indent=2, default=str) + "\n")
