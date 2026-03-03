from __future__ import annotations

import argparse
import json
from pathlib import Path

from .io import config_from_dict, witness_raw_from_dict
from .replay_verify import verify_replay
from .run import emit_artifacts, run_orchestrator_v13


def _load_input(path: Path):
    obj = json.loads(path.read_text())
    config = config_from_dict(obj["config"])
    raw_log = [witness_raw_from_dict(r) for r in obj["raw_log"]]
    return config, raw_log


def main() -> None:
    parser = argparse.ArgumentParser(description="Run deterministic orchestrator_v13 and emit artifacts")
    parser.add_argument("--input", required=True, help="JSON input with {config, raw_log}")
    parser.add_argument("--out-dir", required=True, help="Output directory for four artifacts")
    parser.add_argument("--verify", action="store_true", help="Run replay verification after emit")
    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.out_dir)

    config, raw_log = _load_input(input_path)
    out = run_orchestrator_v13(raw_log, config)
    emit_artifacts(out, out_dir)

    if args.verify:
        errors = verify_replay(raw_log, config, out_dir)
        if errors:
            for e in errors:
                print(f"{e.code} {e.path}")
            raise SystemExit(1)
        print("verify: OK")


if __name__ == "__main__":
    main()
