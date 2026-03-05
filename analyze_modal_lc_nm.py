#!/usr/bin/env python3
"""Run LC + near-miss analysis on Modal for one n.

Outputs worker-merged summary similar to ``analyze.py`` but cloud-scaled.

Examples:
    modal run --detach analyze_modal_lc_nm.py --n 27 --workers 1024 --top-k 200
    modal run analyze_modal_lc_nm.py::collect --n 27 --workers 1024
"""

from __future__ import annotations

import heapq
import json
import time
from typing import Any

import modal

app = modal.App("erdos-993-lc-nm-analysis")

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
    return f"erdos-993-n{n}-lc-nm-results"


@app.function(image=image, timeout=43200, cpu=1)
def analyze_partition(
    n: int,
    res: int,
    mod: int,
    top_k: int,
    lc_top_k: int,
    dict_name: str,
    collect_all_lc: bool,
) -> dict[str, Any]:
    """Analyze one geng partition and persist summary."""
    import sys

    sys.path.insert(0, "/root/src")
    from indpoly import (
        independence_poly,
        is_log_concave,
        is_unimodal,
        log_concavity_ratio,
        near_miss_ratio,
    )
    from trees import trees_geng_raw

    results_dict = modal.Dict.from_name(dict_name, create_if_missing=True)

    count = 0
    non_unimodal_count = 0
    lc_fail_count = 0

    worst_lc_ratio = -1.0
    worst_lc_item: dict[str, Any] | None = None

    lc_heap: list[tuple[float, int, dict[str, Any]]] = []
    nm_heap: list[tuple[float, int, dict[str, Any]]] = []
    all_lc_failures: list[dict[str, Any]] = []

    serial = 0

    for tree_n, adj, g6 in trees_geng_raw(n, res=res, mod=mod):
        count += 1
        poly = independence_poly(tree_n, adj)

        if not is_unimodal(poly):
            non_unimodal_count += 1

        if not is_log_concave(poly):
            lc_fail_count += 1
            lc_ratio, lc_pos = log_concavity_ratio(poly)
            lc_item = {
                "graph6": g6.decode(),
                "n": tree_n,
                "lc_ratio": lc_ratio,
                "lc_pos": lc_pos,
                "poly": poly if collect_all_lc else None,
            }
            if collect_all_lc:
                all_lc_failures.append(lc_item)
            serial += 1
            if len(lc_heap) < lc_top_k:
                heapq.heappush(lc_heap, (lc_ratio, serial, lc_item))
            elif lc_ratio > lc_heap[0][0]:
                heapq.heapreplace(lc_heap, (lc_ratio, serial, lc_item))
            if lc_ratio > worst_lc_ratio:
                worst_lc_ratio = lc_ratio
                worst_lc_item = lc_item

        nm_ratio, nm_pos = near_miss_ratio(poly)
        if nm_ratio > 0:
            nm_item = {
                "graph6": g6.decode(),
                "n": tree_n,
                "nm_ratio": nm_ratio,
                "nm_pos": nm_pos,
                "poly": poly,
            }
            serial += 1
            if len(nm_heap) < top_k:
                heapq.heappush(nm_heap, (nm_ratio, serial, nm_item))
            elif nm_ratio > nm_heap[0][0]:
                heapq.heapreplace(nm_heap, (nm_ratio, serial, nm_item))

    top_lc = [d for _, _, d in sorted(lc_heap, key=lambda x: -x[0])]
    top_nm = [d for _, _, d in sorted(nm_heap, key=lambda x: -x[0])]

    result = {
        "count": count,
        "non_unimodal_count": non_unimodal_count,
        "lc_fail_count": lc_fail_count,
        "worst_lc_ratio": worst_lc_ratio,
        "worst_lc_item": worst_lc_item,
        "top_lc": top_lc,
        "top_nm": top_nm,
        "all_lc_failures": all_lc_failures if collect_all_lc else None,
    }

    results_dict[f"{res}/{mod}"] = result
    return result


@app.function(image=image, timeout=3600, cpu=1)
def launch_partitions(
    n: int,
    workers: int,
    top_k: int,
    lc_top_k: int,
    dict_name: str,
    collect_all_lc: bool,
) -> dict[str, Any]:
    """Launch one detached analysis task per partition from server-side function."""
    for res in range(workers):
        analyze_partition.spawn(
            n, res, workers, top_k, lc_top_k, dict_name, collect_all_lc
        )
    return {
        "launched": workers,
        "n": n,
        "dict_name": dict_name,
        "top_k": top_k,
        "lc_top_k": lc_top_k,
        "collect_all_lc": collect_all_lc,
    }


@app.function(image=image, timeout=300)
def collect_results(n: int, workers: int, dict_name: str) -> dict[str, Any]:
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
        "trees": sum(v["count"] for v in completed.values()),
        "lc_failures": sum(v["lc_fail_count"] for v in completed.values()),
        "non_unimodal": sum(v["non_unimodal_count"] for v in completed.values()),
    }


@app.local_entrypoint()
def main(
    n: int = 27,
    workers: int = 1024,
    top_k: int = 200,
    lc_top_k: int = 200,
    collect_all_lc: bool = False,
    out_json: str = "",
):
    if out_json == "":
        out_json = f"results/analysis_n{n}_modal_lc_nm.json"

    dict_name = _dict_name(n)

    print(
        f"Modal LC/NM analysis: n={n}, workers={workers}, "
        f"top_k={top_k}, lc_top_k={lc_top_k}, collect_all_lc={collect_all_lc}"
    )
    print(f"Dict: {dict_name}")

    t0 = time.time()

    agg_count = 0
    agg_non_unimodal = 0
    agg_lc_fail = 0

    global_worst_lc_ratio = -1.0
    global_worst_lc_item: dict[str, Any] | None = None

    merged_nm: list[tuple[float, int, dict[str, Any]]] = []
    merged_lc: list[tuple[float, int, dict[str, Any]]] = []
    all_lc_failures: list[dict[str, Any]] = []
    serial = 0

    for idx, result in enumerate(
        analyze_partition.starmap(
            [
                (n, res, workers, top_k, lc_top_k, dict_name, collect_all_lc)
                for res in range(workers)
            ]
        ),
        start=1,
    ):
        agg_count += result["count"]
        agg_non_unimodal += result["non_unimodal_count"]
        agg_lc_fail += result["lc_fail_count"]

        if result["worst_lc_ratio"] > global_worst_lc_ratio:
            global_worst_lc_ratio = result["worst_lc_ratio"]
            global_worst_lc_item = result["worst_lc_item"]

        for item in result["top_nm"]:
            serial += 1
            ratio = item["nm_ratio"]
            if len(merged_nm) < top_k:
                heapq.heappush(merged_nm, (ratio, serial, item))
            elif ratio > merged_nm[0][0]:
                heapq.heapreplace(merged_nm, (ratio, serial, item))

        for item in result["top_lc"]:
            serial += 1
            ratio = item["lc_ratio"]
            if len(merged_lc) < lc_top_k:
                heapq.heappush(merged_lc, (ratio, serial, item))
            elif ratio > merged_lc[0][0]:
                heapq.heapreplace(merged_lc, (ratio, serial, item))

        if collect_all_lc and result["all_lc_failures"]:
            all_lc_failures.extend(result["all_lc_failures"])

        if idx % 50 == 0 or idx == workers:
            elapsed = time.time() - t0
            print(
                f"  {idx}/{workers} partitions done, "
                f"trees={agg_count:,}, lc_fail={agg_lc_fail:,}, "
                f"non_unimodal={agg_non_unimodal:,}, {elapsed:.0f}s elapsed"
            )

    elapsed = time.time() - t0

    top_nm = [d for _, _, d in sorted(merged_nm, key=lambda x: -x[0])]
    top_lc = [d for _, _, d in sorted(merged_lc, key=lambda x: -x[0])]

    summary = {
        "n": n,
        "workers": workers,
        "total_trees": agg_count,
        "non_unimodal_count": agg_non_unimodal,
        "lc_failure_count": agg_lc_fail,
        "worst_lc_ratio": global_worst_lc_ratio,
        "worst_lc_item": global_worst_lc_item,
        "top_near_misses": top_nm,
        "top_lc_failures": top_lc,
        "all_lc_failures": all_lc_failures if collect_all_lc else None,
        "wall_time_s": round(elapsed, 1),
        "platform": "Modal",
        "dict_name": dict_name,
        "date": time.strftime("%Y-%m-%d"),
    }

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\nDone")
    print(json.dumps({
        "n": n,
        "workers": workers,
        "total_trees": agg_count,
        "non_unimodal_count": agg_non_unimodal,
        "lc_failure_count": agg_lc_fail,
        "worst_lc_ratio": global_worst_lc_ratio,
        "wall_time_s": round(elapsed, 1),
    }, indent=2))
    print(f"Saved: {out_json}")


@app.local_entrypoint()
def dispatch(
    n: int = 27,
    workers: int = 1024,
    top_k: int = 200,
    lc_top_k: int = 200,
    collect_all_lc: bool = False,
):
    """Fire-and-forget launch: spawn one analysis task per partition."""
    dict_name = _dict_name(n)
    print(
        f"Dispatching LC/NM tasks for n={n}, workers={workers}, "
        f"top_k={top_k}, lc_top_k={lc_top_k}, collect_all_lc={collect_all_lc}, "
        f"dict={dict_name}"
    )
    for res in range(workers):
        analyze_partition.spawn(
            n, res, workers, top_k, lc_top_k, dict_name, collect_all_lc
        )
    print("Dispatch submitted.")
