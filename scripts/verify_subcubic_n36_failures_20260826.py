#!/usr/bin/env python3
"""Independently replay the n=36 subcubic census certificates.

This verifier does not import the census engine. It sums the 32 shard STATS
records, rebuilds each FAIL tree from its parent array, recomputes its
independence polynomial by a separate tree DP, and audits the centre/arm
structure and rooted-arm gap data used in the accompanying note.
"""

import argparse
import glob
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

import networkx as nx


STATS_RE = re.compile(
    r"^STATS n=(\d+) alpha_max=(\d+) trees=(\d+) passed=(\d+) "
    r"failures=(\d+) alarms=(\d+)$",
    re.MULTILINE,
)
FAIL_RE = re.compile(r"^FAIL .+$", re.MULTILINE)


def poly_add(a: list[int], b: list[int]) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return out


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            out[i + j] += left * right
    return out


def tree_states(graph: nx.Graph, root: int) -> tuple[list[int], list[int]]:
    """Return (root free, root excluded) independence polynomials."""
    parent = {root: -1}
    order = [root]
    for vertex in order:
        for neighbour in graph[vertex]:
            if neighbour != parent[vertex]:
                parent[neighbour] = vertex
                order.append(neighbour)
    assert len(order) == len(graph), "tree DP did not reach every vertex"

    included: dict[int, list[int]] = {}
    excluded: dict[int, list[int]] = {}
    for vertex in reversed(order):
        inc = [0, 1]
        exc = [1]
        for child in graph[vertex]:
            if parent.get(child) == vertex:
                inc = poly_mul(inc, excluded[child])
                exc = poly_mul(exc, poly_add(included[child], excluded[child]))
        included[vertex] = inc
        excluded[vertex] = exc
    return poly_add(included[root], excluded[root]), excluded[root]


def rooted_code(graph: nx.Graph, root: int, parent: int = -1) -> str:
    children = sorted(
        rooted_code(graph, child, root)
        for child in graph[root]
        if child != parent
    )
    return "(" + "".join(children) + ")"


def is_unimodal(sequence: list[int]) -> bool:
    descending = False
    for left, right in zip(sequence, sequence[1:]):
        if right < left:
            descending = True
        elif right > left and descending:
            return False
    return True


def shard_number(path: str) -> int:
    return int(re.search(r"_s(\d+)of32", path).group(1))


def parse_failure(line: str) -> tuple[list[int], list[int]]:
    parents = [int(x) for x in re.search(r" par=([^ ]+)", line).group(1).split(",")]
    polynomial = [int(x) for x in re.search(r" poly=([^ ]+)", line).group(1).split(",")]
    return parents, polynomial


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--census-dir",
        default="results/subcubic_census_20260826",
        type=Path,
    )
    parser.add_argument(
        "--output",
        default="results/subcubic_n36_witness_audit_20260826.json",
        type=Path,
    )
    args = parser.parse_args()

    summary_path = args.census_dir / "summary_n36.json"
    summary = json.loads(summary_path.read_text())
    paths = sorted(
        glob.glob(str(args.census_dir / "n36_s*of32.txt")),
        key=shard_number,
    )
    assert len(paths) == 32, f"expected 32 shards, found {len(paths)}"

    totals = Counter()
    failures: list[tuple[str, str]] = []
    alarm_lines: list[str] = []
    for path in paths:
        text = Path(path).read_text()
        stats = STATS_RE.findall(text)
        assert len(stats) == 1, f"{path}: expected one STATS record"
        n, _, trees, passed, failed, alarms = map(int, stats[0])
        assert n == 36
        totals.update(
            trees=trees,
            passed=passed,
            failures=failed,
            alarms=alarms,
        )
        failures.extend((Path(path).name, line) for line in FAIL_RE.findall(text))
        alarm_lines.extend(re.findall(r"^ALARM .+$", text, re.MULTILINE))

    expected = summary["expected_trees_A000672"]
    assert totals["trees"] == expected == summary["trees"]
    assert totals["passed"] == expected
    assert totals["failures"] == len(failures) == summary["failures"] == 15
    assert totals["alarms"] == len(alarm_lines) == summary["alarms"] == 0
    assert sorted(line for _, line in failures) == sorted(summary["fail_lines"])

    records = []
    graphs = []
    for number, (shard, line) in enumerate(failures, 1):
        parents, saved_polynomial = parse_failure(line)
        assert len(parents) == 36 and parents[0] == 0
        graph = nx.Graph()
        graph.add_nodes_from(range(36))
        graph.add_edges_from((i, parents[i] - 1) for i in range(1, 36))
        assert nx.is_tree(graph)
        assert max(dict(graph.degree()).values()) == 3

        polynomial, _ = tree_states(graph, 0)
        assert polynomial == saved_polynomial
        breaks = [
            k
            for k in range(1, len(polynomial) - 1)
            if polynomial[k] ** 2 < polynomial[k - 1] * polynomial[k + 1]
        ]
        alpha = len(polynomial) - 1
        threshold = math.ceil((2 * alpha - 1) / 3)

        centres = nx.center(graph)
        assert len(centres) == 1
        centre = centres[0]
        assert graph.degree(centre) == 3
        without_centre = graph.copy()
        without_centre.remove_node(centre)

        arms = []
        for component in nx.connected_components(without_centre):
            vertices = set(component)
            root = next(v for v in graph[centre] if v in vertices)
            arm_graph = without_centre.subgraph(vertices)
            free, root_excluded = tree_states(arm_graph, root)
            code = rooted_code(arm_graph, root)
            arms.append(
                {
                    "size": len(vertices),
                    "gap": (len(free) - 1) - (len(root_excluded) - 1),
                    "alpha_free": len(free) - 1,
                    "alpha_root_excluded": len(root_excluded) - 1,
                    "shape_id": hashlib.sha256(code.encode()).hexdigest()[:12],
                }
            )
        arms.sort(key=lambda arm: (arm["size"], arm["shape_id"]))

        record = {
            "number": number,
            "shard": shard,
            "centre": centre + 1,
            "arm_sizes": [arm["size"] for arm in arms],
            "arm_gaps": [arm["gap"] for arm in arms],
            "arm_shape_ids": [arm["shape_id"] for arm in arms],
            "distinct_rooted_arm_shapes": len({arm["shape_id"] for arm in arms}),
            "branch_vertices": sum(graph.degree(v) == 3 for v in graph),
            "leaves": sum(graph.degree(v) == 1 for v in graph),
            "alpha": alpha,
            "breaks": breaks,
            "distances": [k - threshold for k in breaks],
            "unimodal": is_unimodal(polynomial),
        }
        assert record["alpha"] == 19
        assert record["breaks"] == [18]
        assert record["distances"] == [5]
        assert record["unimodal"]
        records.append(record)
        graphs.append(graph)

    duplicate_pairs = []
    for right in range(len(graphs)):
        for left in range(right):
            if nx.is_isomorphic(graphs[left], graphs[right]):
                duplicate_pairs.append([left + 1, right + 1])
    assert not duplicate_pairs

    shape_counts = Counter(r["distinct_rooted_arm_shapes"] for r in records)
    gap_counts = Counter(tuple(r["arm_gaps"]) for r in records)
    assert shape_counts == Counter({3: 10, 2: 5})
    assert gap_counts == Counter({(0, 0, 0): 9, (1, 1, 1): 6})

    result = {
        "n": 36,
        "shards": len(paths),
        "trees": totals["trees"],
        "expected_trees_A000672": expected,
        "failures": len(failures),
        "alarms": len(alarm_lines),
        "duplicate_witness_pairs": duplicate_pairs,
        "rooted_shape_count_distribution": dict(sorted(shape_counts.items())),
        "arm_gap_vector_distribution": {
            ",".join(map(str, key)): value for key, value in sorted(gap_counts.items())
        },
        "arm_size_range": [
            min(size for record in records for size in record["arm_sizes"]),
            max(size for record in records for size in record["arm_sizes"]),
        ],
        "records": records,
    }
    args.output.write_text(json.dumps(result, indent=1) + "\n")
    print(
        "verified: 32 shards, 25,664,800,714 trees, 15 distinct failures, "
        "0 alarms"
    )
    print("verified: all polynomials replay; alpha=19, break=[18], dist=[5], unimodal")
    print("rooted arm shapes: 5 witnesses with two shapes; 10 with three")
    print("arm gaps: 9 witnesses with (0,0,0); 6 with (1,1,1)")
    print(f"saved {args.output}")


if __name__ == "__main__":
    main()
