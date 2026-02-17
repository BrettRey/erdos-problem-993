#!/usr/bin/env python3
"""Run n=27 exhaustive unimodality check on Modal cloud compute.

Usage:
    modal run --detach search_modal.py

Each worker saves its result to a Modal Dict as it finishes,
so completed work is never lost. A separate collector reads
the dict and produces the final summary.
"""

import json
import time

import modal

app = modal.App("erdos-993-n27")

# Persistent dict to store results as they complete
results_dict = modal.Dict.from_name("erdos-993-n27-results", create_if_missing=True)

# Container image: build nauty from source + numpy + source files
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


@app.function(image=image, timeout=7200, cpu=1)
def check_partition(n: int, res: int, mod: int) -> dict:
    """Check one geng partition for unimodality violations.

    Saves result to Modal Dict immediately on completion.
    """
    import sys

    sys.path.insert(0, "/root/src")
    from indpoly import independence_poly, is_unimodal
    from trees import trees_geng

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


@app.local_entrypoint()
def main():
    n = 27
    mod = 1024
    expected = 751_065_460

    print(f"Erdos Problem #993: n={n} exhaustive search on Modal")
    print(f"Launching {mod} parallel workers...")
    print(f"Results saved to Modal Dict 'erdos-993-n27-results' as they complete.")
    print(f"Run 'modal run search_modal.py::collect' to gather results at any time.")
    print()

    t0 = time.time()
    completed = 0
    total_trees = 0
    counterexamples = []

    for result in check_partition.starmap(
        [(n, res, mod) for res in range(mod)]
    ):
        completed += 1
        total_trees += result["count"]
        if result["counterexample"] is not None:
            counterexamples.append(result["counterexample"])
        if completed % 50 == 0 or completed == mod:
            elapsed = time.time() - t0
            print(
                f"  {completed}/{mod} partitions done, "
                f"{total_trees:,} trees so far, "
                f"{elapsed:.0f}s elapsed"
            )

    elapsed = time.time() - t0
    print(f"\nn={n}: {total_trees:,} trees checked in {elapsed:.1f}s")
    if total_trees == expected:
        print(f"Tree count matches OEIS A000055 ({expected:,})")
    else:
        print(f"WARNING: count mismatch! Got {total_trees:,}, expected {expected:,}")

    if counterexamples:
        print(f"\n*** COUNTEREXAMPLE FOUND ***")
        for cx in counterexamples:
            print(json.dumps(cx, indent=2))
    else:
        print(f"All unimodal. Conjecture holds for n <= {n}.")

    print(f"\nWall time: {elapsed:.1f}s ({elapsed / 60:.1f} minutes)")

    summary = {
        "n": n,
        "total_trees": total_trees,
        "expected": expected,
        "match": total_trees == expected,
        "counterexamples": len(counterexamples),
        "wall_time_s": round(elapsed, 1),
        "workers": mod,
    }
    with open(f"results/analysis_n{n}_modal.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved to results/analysis_n{n}_modal.json")


@app.function(image=image, timeout=300)
def collect_results(mod: int = 1024) -> dict:
    """Collect completed results from the Modal Dict."""
    completed = {}
    for res in range(mod):
        key = f"{res}/{mod}"
        try:
            completed[key] = results_dict[key]
        except KeyError:
            pass
    return {
        "completed": len(completed),
        "total": mod,
        "trees": sum(r["count"] for r in completed.values()),
        "counterexamples": [
            r["counterexample"]
            for r in completed.values()
            if r["counterexample"] is not None
        ],
    }
