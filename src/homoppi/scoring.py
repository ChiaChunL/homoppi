"""Score integration shared by the IM and DDI methods."""

from __future__ import annotations

from collections.abc import Iterable


def bayes_integration(scores: Iterable[float]) -> float:
    """Bayesian integration of independent evidence: S = 1 - prod(1 - s_i).

    Returns 0.0 for an empty evidence list.
    """
    product = 1.0
    for s in scores:
        product *= 1.0 - s
    return 1.0 - product
