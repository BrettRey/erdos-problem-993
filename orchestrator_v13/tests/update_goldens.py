from __future__ import annotations

from pathlib import Path

from orchestrator_v13 import emit_artifacts, run_orchestrator_v13
from orchestrator_v13.io import load_fixture

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures"
GOLDEN_ROOT = Path(__file__).resolve().parent / "goldens"


def main() -> None:
    GOLDEN_ROOT.mkdir(parents=True, exist_ok=True)
    for fixture_path in sorted(FIXTURE_DIR.glob("T*.json")):
        config, raw_log, _expected = load_fixture(fixture_path)
        out = run_orchestrator_v13(raw_log, config)
        out_dir = GOLDEN_ROOT / fixture_path.stem
        emit_artifacts(out, out_dir)
        print(f"wrote {fixture_path.stem} -> {out_dir}")


if __name__ == "__main__":
    main()
