"""Reading and validating query pair lists (the batch-inference input)."""

from __future__ import annotations

import logging
from pathlib import Path

HEADER_ALIASES = {("protein_a", "protein_b"), ("query_a", "query_b"), ("proteina", "proteinb")}


def read_pairs(path: Path, logger: logging.Logger) -> list[tuple[str, str]]:
    """Read a two-column TSV of query pairs.

    A header row is optional and detected by name. Only the first two columns
    are used. Duplicate pairs (in either orientation) are collapsed to their
    first occurrence; input order and orientation are preserved in the output.
    """
    pairs: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    n_dupes = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{lineno}: expected at least 2 tab-separated columns")
            a, b = fields[0].strip(), fields[1].strip()
            if lineno == 1 and (a.lower(), b.lower()) in HEADER_ALIASES:
                continue
            if not a or not b:
                raise ValueError(f"{path}:{lineno}: empty protein ID")
            key = (a, b) if a <= b else (b, a)
            if key in seen:
                n_dupes += 1
                continue
            seen.add(key)
            pairs.append((a, b))
    if n_dupes:
        logger.warning("Collapsed %s duplicate pairs from %s", f"{n_dupes:,}", path)
    logger.info("Loaded %s query pairs from %s", f"{len(pairs):,}", path)
    return pairs
