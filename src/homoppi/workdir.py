"""Run directory layout and stage bookkeeping (the basis for resume)."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path


class Workdir:
    """A per-run working directory with cached stages.

    Layout:
        <workdir>/state.json      stage completion records
        <workdir>/blast/          raw + filtered blastp results
        <workdir>/hmmscan/        raw + filtered hmmscan results (phase 2)
        <workdir>/results/        final tables + parameter snapshots
        <workdir>/logs/           external tool logs and run logs
    """

    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.mkdir(parents=True, exist_ok=True)

    def _subdir(self, name: str) -> Path:
        d = self.path / name
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def blast_dir(self) -> Path:
        return self._subdir("blast")

    @property
    def hmmscan_dir(self) -> Path:
        return self._subdir("hmmscan")

    @property
    def results_dir(self) -> Path:
        return self._subdir("results")

    @property
    def logs_dir(self) -> Path:
        return self._subdir("logs")

    @property
    def _state_path(self) -> Path:
        return self.path / "state.json"

    def _load_state(self) -> dict:
        if self._state_path.exists():
            return json.loads(self._state_path.read_text())
        return {"stages": {}}

    def _save_state(self, state: dict) -> None:
        self._state_path.write_text(json.dumps(state, indent=2) + "\n")

    @staticmethod
    def _fingerprint(path: Path) -> dict:
        st = path.stat()
        return {"path": str(path.resolve()), "size": st.st_size, "mtime": int(st.st_mtime)}

    def get_stage(self, stage: str) -> dict | None:
        """The completion record of a stage, or None if it never completed."""
        return self._load_state()["stages"].get(stage)

    def inputs_current(self, record: dict, inputs: dict[str, Path]) -> bool:
        """True if every input file still matches the fingerprints in the record."""
        for name, path in inputs.items():
            if not path.exists():
                return False
            if record.get("inputs", {}).get(name) != self._fingerprint(path):
                return False
        return True

    def is_stage_current(self, stage: str, params: dict, inputs: dict[str, Path]) -> bool:
        """True if the stage completed earlier with identical params and unchanged inputs."""
        record = self.get_stage(stage)
        if record is None:
            return False
        if record.get("params") != params:
            return False
        return self.inputs_current(record, inputs)

    def mark_stage(self, stage: str, params: dict, inputs: dict[str, Path]) -> None:
        state = self._load_state()
        state["stages"][stage] = {
            "completed_at": datetime.now(UTC).isoformat(timespec="seconds"),
            "params": params,
            "inputs": {name: self._fingerprint(path) for name, path in inputs.items()},
        }
        self._save_state(state)
