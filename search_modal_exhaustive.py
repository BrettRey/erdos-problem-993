#!/usr/bin/env python3
"""Run exhaustive unimodality checks on Modal cloud.

This generalizes ``search_modal.py`` so you can run any single n.

Examples:
    modal run --detach search_modal_exhaustive.py --n 28 --workers 1024
    modal run search_modal_exhaustive.py::collect --n 28 --workers 1024
"""

from __future__ import annotations

import json
import time
from typing import Any

import modal

app = modal.App("erdos-993-exhaustive")

image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("build-essential", "curl")
    .run_commands(
        "curl -sL https://pallini.di.uniroma1.it/nauty2_8_9.tar.gz | tar xz",
        "cd nauty2_8_9 && ./configure --quiet && make -j$(nproc) geng",
        "cp nauty2_8_9/geng /usr/local/bin/geng",
        "rm -rf nauty2_8_9",
    )
    .pip_install("numpy")
    .add_local_file("graph6.py", "/root/src/graph6.py")
    .add_local_file("indpoly.py", "/root/src/indpoly.py")
    .add_local_file("trees.py", "/root/src/trees.py")
)


def _dict_name(n: int) -> str:
    return f"erdos-993-n{n}-unimodality-results"


@app.function(image=image, timeout=43200, cpu=1)
def check_partition(n: int, res: int, mod: int, dict_name: str) -> dict[str, Any]:
    """Check one geng partition and write result to a persistent Modal Dict."""
    import sys

    sys.path.insert(0, "/root/src")
    from indpoly import independence_poly, is_unimodal
    from trees import trees_geng

    results_dict = modal.Dict.from_name(dict_name, create_if_missing=True)

    count = 0
    for tree_n, adj in trees_geng(n, res=res, mod=mod):
        poly = independence_poly(tree_n, adj)
        if not is_unimodal(poly):
            result = {
                "count": count + 1,
                "counterexample": {
                    "n": tree_n,
                    "adj": [list(a) for a in adj],
                    "poly": list(poly),
                },
            }
            results_dict[f"{res}/{mod}"] = result
            return result
        count += 1

    result = {"count": count, "counterexample": None}
    results_dict[f"{res}/{mod}"] = result
    return result


@app.function(image=image, timeout=3600, cpu=1)
def launch_partitions(n: int, workers: int, dict_name: str) -> dict[str, Any]:
    """Launch one detached worker task per partition from a server-side function."""
    for res in range(workers):
        check_partition.spawn(n, res, workers, dict_name)
    return {"launched": workers, "n": n, "dict_name": dict_name}


@app.function(image=image, timeout=300)
def collect_results(n: int, workers: int, dict_name: str) -> dict[str, Any]:
    """Collect completed results from the persistent dict."""
    results_dict = modal.Dict.from_name(dict_name, create_if_missing=True)
    completed: dict[str, Any] = {}
    for res in range(workers):
        key = f"{res}/{workers}"
        try:
            completed[key] = results_dict[key]
        except KeyError:
            pass
    return {
        "n": n,
        "completed": len(completed),
        "workers": workers,
        "trees": sum(r["count"] for r in completed.values()),
        "counterexamples": [
            r["counterexample"]
            for r in completed.values()
            if r["counterexample"] is not None
        ],
    }


@app.local_entrypoint()
def main(n: int = 28, workers: int = 1024, expected: int = 0, out_json: str = ""):
    """Launch exhaustive run and wait for completion."""
    if out_json == "":
        out_json = f"results/analysis_n{n}_modal_unimodality.json"

    dict_name = _dict_name(n)

    print(f"Erdos #993 exhaustive unimodality on Modal: n={n}, workers={workers}")
    print(f"Dict: {dict_name}")
    print("Launching partitions...")

    t0 = time.time()
    completed = 0
    total_trees = 0
    counterexamples = []

    for result in check_partition.starmap(
        [(n, res, workers, dict_name) for res in range(workers)]
    ):
        completed += 1
        total_trees += result["count"]
        if result["counterexample"] is not None:
            counterexamples.append(result["counterexample"])
        if completed % 50 == 0 or completed == workers:
            elapsed = time.time() - t0
            print(
                f"  {completed}/{workers} partitions done, "
                f"{total_trees:,} trees, {elapsed:.0f}s elapsed"
            )

    elapsed = time.time() - t0

    summary = {
        "n": n,
        "workers": workers,
        "total_trees": total_trees,
        "expected": expected if expected > 0 else None,
        "match": (total_trees == expected) if expected > 0 else None,
        "counterexamples": len(counterexamples),
        "counterexample_examples": counterexamples[:3],
        "wall_time_s": round(elapsed, 1),
        "platform": "Modal",
        "dict_name": dict_name,
        "date": time.strftime("%Y-%m-%d"),
    }

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\nDone")
    print(json.dumps(summary, indent=2))
    print(f"Saved: {out_json}")


@app.local_entrypoint()
def dispatch(n: int = 28, workers: int = 1024):
    """Fire-and-forget launch: spawn one task per partition and exit."""
    dict_name = _dict_name(n)
    print(
        f"Dispatching unimodality tasks for n={n}, workers={workers}, "
        f"dict={dict_name}"
    )
    for res in range(workers):
        check_partition.spawn(n, res, workers, dict_name)
    print("Dispatch submitted.")
