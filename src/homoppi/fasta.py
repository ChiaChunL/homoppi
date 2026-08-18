"""Minimal FASTA reading/writing and sequence-ID extraction."""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from pathlib import Path

UNIPROT_HEADER_RE = re.compile(r"^(?:sp|tr)\|([^|\s]+)\|")


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) tuples; header excludes the leading '>'."""
    header: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def extract_id(header: str, id_regex: re.Pattern[str] | None = None) -> str:
    """Extract a protein ID from a FASTA header.

    Priority: user-supplied regex (first capture group) > UniProt-style
    'sp|ACC|NAME' accession > first whitespace-delimited token.
    """
    if id_regex is not None:
        m = id_regex.search(header)
        if m:
            return m.group(1)
    m = UNIPROT_HEADER_RE.match(header)
    if m:
        return m.group(1)
    return header.split()[0]


def write_fasta(records: Iterable[tuple[str, str]], path: Path, width: int = 60) -> int:
    """Write records as FASTA; return the number of sequences written."""
    n = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")
            n += 1
    return n


def count_sequences(path: Path) -> int:
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))
