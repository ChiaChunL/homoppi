"""Logging setup: console always, plus an optional per-run log file."""

from __future__ import annotations

import logging
import sys
from pathlib import Path


def get_logger(log_file: Path | None = None, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger("homoppi")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s", "%H:%M:%S"))
    logger.addHandler(console)

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file)
        fh.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-7s | %(message)s"))
        logger.addHandler(fh)

    return logger
