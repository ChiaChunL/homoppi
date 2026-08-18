"""Discovery and execution of external tools (blastp, makeblastdb, hmmscan...)."""

from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from pathlib import Path

INSTALL_HINTS = {
    "blastp": "conda install -c bioconda blast",
    "makeblastdb": "conda install -c bioconda blast",
    "hmmscan": "conda install -c bioconda hmmer",
    "hmmpress": "conda install -c bioconda hmmer",
}


class ExternalToolError(RuntimeError):
    pass


def find_binary(name: str, override: Path | None = None) -> str:
    """Resolve an external binary from an explicit override or PATH."""
    if override is not None:
        if not override.exists():
            raise ExternalToolError(f"{name} binary not found at {override}")
        return str(override)
    found = shutil.which(name)
    if found is None:
        hint = INSTALL_HINTS.get(name, "")
        raise ExternalToolError(
            f"'{name}' was not found on PATH. Install it (e.g. `{hint}`) or point to it with --{name}-bin."
        )
    return found


def run_command(cmd: list[str], log_path: Path, logger: logging.Logger) -> None:
    """Run a command, teeing stdout+stderr to a log file; raise with the log tail on failure."""
    logger.info("Running: %s", shlex.join(cmd))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_fh:
        proc = subprocess.run(cmd, stdout=log_fh, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        tail = "".join(log_path.read_text().splitlines(keepends=True)[-20:])
        raise ExternalToolError(
            f"command failed (exit {proc.returncode}): {shlex.join(cmd)}\n"
            f"--- last lines of {log_path} ---\n{tail}"
        )
