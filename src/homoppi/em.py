"""EM expansion of the DDI library from known PPIs (Deng et al. model).

Candidate domain pairs are those co-occurring across known interacting
proteins. For each candidate (m, n), lambda_mn = P(domains m and n interact)
is estimated by expectation-maximization over all same-species protein pairs
that contain the domain pair, with a fixed false-negative rate and a derived
(or user-supplied) false-positive rate.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


@dataclass(frozen=True)
class EMParams:
    false_neg_rate: float = 0.8
    false_pos_rate: float | None = None  # derived from library counts when None
    max_iter: int = 50
    tol: float = 1e-4
    max_proteins_per_domain: int = 1000


def load_domain_annotations(paths: list[Path]) -> dict[str, set[str]]:
    """Merge domainanno outputs (domains.tsv) into {protein: {pfam_acc}}."""
    annotations: dict[str, set[str]] = {}
    for path in paths:
        df = pd.read_csv(path, sep="\t", dtype={"query_id": str, "pfam_acc": str})
        for row in df.itertuples(index=False):
            annotations.setdefault(row.query_id, set()).add(row.pfam_acc)
    return annotations


def derive_false_pos_rate(n_ppis: int, n_proteins: int, false_neg_rate: float) -> float:
    """False-positive rate consistent with the observed interaction density."""
    avg_partner = n_ppis / n_proteins
    fp = (2 * n_ppis - (1 - false_neg_rate) * n_proteins * avg_partner) / (
        (n_proteins + 1 - avg_partner) * n_proteins
    )
    return max(fp, 1e-12)


def run_ddi_em(
    ppi_pairs: set[tuple[str, str]],
    protein_taxid: dict[str, int],
    protein_domains: dict[str, set[str]],
    params: EMParams,
    logger: logging.Logger,
) -> pd.DataFrame:
    """Estimate lambda for every candidate domain pair; returns pfam_a/pfam_b/em_score/... rows."""
    # Universe: proteins with both a taxid and at least one domain.
    proteins = {p for p in protein_domains if p in protein_taxid}
    ppis = {tuple(sorted(pair)) for pair in ppi_pairs if pair[0] in proteins and pair[1] in proteins}
    logger.info(
        "EM universe: %s annotated proteins, %s known PPIs with both sides annotated",
        f"{len(proteins):,}", f"{len(ppis):,}",
    )
    if not ppis:
        raise ValueError("no known PPIs remain after restricting to domain-annotated proteins.")

    dom2pros: dict[str, list[str]] = {}
    for protein in proteins:
        for domain in protein_domains[protein]:
            dom2pros.setdefault(domain, []).append(protein)
    promiscuous = {d for d, ps in dom2pros.items() if len(ps) > params.max_proteins_per_domain}
    if promiscuous:
        logger.warning(
            "Excluding %s promiscuous domains present in > %s proteins.",
            len(promiscuous), params.max_proteins_per_domain,
        )

    # Candidate domain pairs: co-occurring across known interacting proteins.
    candidates: set[tuple[str, str]] = set()
    for p, q in ppis:
        for dm in protein_domains[p]:
            if dm in promiscuous:
                continue
            for dn in protein_domains[q]:
                if dn in promiscuous:
                    continue
                candidates.add((dm, dn) if dm <= dn else (dn, dm))
    logger.info("Candidate domain pairs: %s", f"{len(candidates):,}")

    # Protein-pair universe per candidate: same-species pairs containing the domain pair.
    pair_index: dict[tuple[str, str], int] = {}
    observed: list[int] = []
    pair_dps: list[list[int]] = []
    dp_pair_idx: list[np.ndarray] = []
    dp_keys: list[tuple[str, str]] = []

    for dm, dn in tqdm(sorted(candidates), desc="Collating protein pairs", unit="dp"):
        pros_m, pros_n = dom2pros[dm], dom2pros[dn]
        pairs_mn: set[tuple[str, str]] = set()
        for i in pros_m:
            taxid_i = protein_taxid[i]
            for j in pros_n:
                if i == j or protein_taxid[j] != taxid_i:
                    continue
                pairs_mn.add((i, j) if i <= j else (j, i))
        if not pairs_mn:
            continue
        dp_idx = len(dp_keys)
        dp_keys.append((dm, dn))
        indices = np.empty(len(pairs_mn), dtype=np.int64)
        for k, pair in enumerate(pairs_mn):
            idx = pair_index.get(pair)
            if idx is None:
                idx = pair_index[pair] = len(observed)
                observed.append(1 if pair in ppis else 0)
                pair_dps.append([])
            pair_dps[idx].append(dp_idx)
            indices[k] = idx
        dp_pair_idx.append(indices)

    n_pairs = len(observed)
    obs = np.array(observed, dtype=np.float64)
    logger.info("Protein-pair universe: %s pairs (%s observed interacting)", f"{n_pairs:,}", f"{int(obs.sum()):,}")

    fn = params.false_neg_rate
    fp = params.false_pos_rate
    if fp is None:
        fp = derive_false_pos_rate(len(ppis), len(proteins), fn)
    logger.info("EM parameters: false_neg_rate=%s, false_pos_rate=%.3g", fn, fp)

    # Initialization: lambda_mn = fraction of containing pairs that are known PPIs.
    lam = np.array([obs[idx].mean() for idx in dp_pair_idx], dtype=np.float64)

    for iteration in range(1, params.max_iter + 1):
        # P(pair interacts) via its domain pairs, then P(pair observed).
        log_not_inter = np.zeros(n_pairs)
        for dp_idx, indices in enumerate(dp_pair_idx):
            log_not_inter[indices] += np.log1p(-min(lam[dp_idx], 1 - 1e-12))
        p_inter = 1.0 - np.exp(log_not_inter)
        p_obs = np.clip((1.0 - p_inter) * fp + p_inter * (1.0 - fn), 1e-12, 1 - 1e-12)

        # E-step weight per protein pair; M-step: lambda = mean expectation.
        weight = np.where(obs == 1.0, (1.0 - fn) / p_obs, fn / (1.0 - p_obs))
        lam_new = np.array(
            [min(lam[dp_idx] * weight[indices].mean(), 1.0) for dp_idx, indices in enumerate(dp_pair_idx)]
        )
        delta = float(np.max(np.abs(lam_new - lam))) if len(lam) else 0.0
        lam = lam_new
        logger.info("EM iteration %s: max |delta lambda| = %.6f", iteration, delta)
        if delta < params.tol:
            logger.info("Converged after %s iterations.", iteration)
            break
    else:
        logger.warning("EM did not converge within %s iterations (last delta %.6f).", params.max_iter, delta)

    rows = [
        {
            "pfam_a": dp_keys[i][0],
            "pfam_b": dp_keys[i][1],
            "em_score": round(float(lam[i]), 6),
            "n_protein_pairs": int(len(dp_pair_idx[i])),
            "n_interacting_pairs": int(obs[dp_pair_idx[i]].sum()),
        }
        for i in range(len(dp_keys))
    ]
    return pd.DataFrame(rows, columns=["pfam_a", "pfam_b", "em_score", "n_protein_pairs", "n_interacting_pairs"])
