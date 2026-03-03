from __future__ import annotations

from pathlib import Path

from .run import run_orchestrator_v13
from .types import Config, VerifierError, WitnessRaw
from .util import canonical_json_dumps


def verify_replay(raw_log: list[WitnessRaw], config: Config, out_dir: Path) -> list[VerifierError]:
    errors: list[VerifierError] = []

    expected = run_orchestrator_v13(raw_log, config)
    expected_files = {
        "frontier_state.json": canonical_json_dumps(expected.frontier_state),
        "partition_plan.json": canonical_json_dumps(expected.partition_plan),
        "repair_queue.json": canonical_json_dumps(expected.repair_queue),
        "v13_obligations.json": canonical_json_dumps(expected.v13_obligations),
    }

    for name, exp in expected_files.items():
        path = out_dir / name
        if not path.exists():
            errors.append(
                VerifierError(
                    code="E_MISSING_ARTIFACT",
                    path=f"/{name}",
                    expected="present",
                    got="missing",
                )
            )
            continue
        got = path.read_bytes()
        if got != exp:
            errors.append(
                VerifierError(
                    code="E_ARTIFACT_MISMATCH",
                    path=f"/{name}",
                    expected=exp.decode("utf-8"),
                    got=got.decode("utf-8"),
                )
            )

    errors.sort(key=lambda e: (e.code, e.path))
    return errors
