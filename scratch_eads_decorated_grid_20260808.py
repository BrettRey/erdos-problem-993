#!/usr/bin/env python3
"""Exact parameter grid around the first two-certificate decorated spider."""

from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

from scratch_eads_decorated_spider_optimizer_20260808 import build, canonical
from scratch_eads_pareto_optimizer_20260808 import (
    encode,
    make_individual,
    scaled_total_key,
    total_key,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-load", type=int, default=1)
    parser.add_argument("--max-load", type=int, default=15)
    parser.add_argument(
        "--template", choices=("one_anchor", "two_endpoint", "mixed_anchor"),
        default="one_anchor",
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("results/eads_decorated_grid_20260808.json"),
    )
    args = parser.parse_args()
    best = None
    scaled_best = None
    eads_counterexample = None
    main_counterexample = None
    evaluated = 0
    hist: dict[int, int] = {}
    templates = {
        "one_anchor": lambda a, b, c, d: [(a, 0, b), (c,), (d, 1)],
        "two_endpoint": lambda a, b, c, d: [(a, 0, b), (c, 1), (d, 1)],
        "mixed_anchor": lambda a, b, c, d: [(a, 1, b), (c,), (d, 1)],
    }
    values = range(args.min_load, args.max_load + 1)
    for a, b, c, d in product(values, repeat=4):
        spec = canonical(0, templates[args.template](a, b, c, d))
        item = make_individual(build(spec), f"a{a}:b{b}:c{c}:d{d}")
        item["spec"] = spec
        evaluated += 1
        ev = item["evaluation"]
        assert isinstance(ev, dict)
        good = int(ev["good_count"])
        hist[good] = hist.get(good, 0) + 1
        if best is None or total_key(item) < total_key(best):
            best = item
        if scaled_best is None or scaled_total_key(item) < scaled_total_key(scaled_best):
            scaled_best = item
        if not ev["unimodal"]:
            main_counterexample = item
            break
        if good == 0:
            eads_counterexample = item
            break
        if evaluated % 5000 == 0:
            print(json.dumps({
                "event": "progress", "evaluated": evaluated,
                "best_good": best["evaluation"]["good_count"],
                "scaled_best_n": scaled_best["n"],
            }), flush=True)
    assert best is not None and scaled_best is not None

    def record(item: dict[str, object] | None):
        return None if item is None else {"spec": item["spec"], **encode(item)}

    result = {
        "claim_scope": "exact finite parameter grid; evidence only",
        "template": args.template,
        "load_range": [args.min_load, args.max_load],
        "evaluated": evaluated,
        "good_count_histogram": {str(k): v for k, v in sorted(hist.items())},
        "erdos_993_counterexample_found": main_counterexample is not None,
        "eads_counterexample_found": eads_counterexample is not None,
        "erdos_993_counterexample": record(main_counterexample),
        "eads_counterexample": record(eads_counterexample),
        "best": record(best),
        "scaled_best": record(scaled_best),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "event": "done", "evaluated": evaluated,
        "best_good": best["evaluation"]["good_count"],
        "eads_counterexample_found": eads_counterexample is not None,
        "erdos_993_counterexample_found": main_counterexample is not None,
        "out": str(args.out),
    }), flush=True)
    raise SystemExit(1 if eads_counterexample or main_counterexample else 0)


if __name__ == "__main__":
    main()
