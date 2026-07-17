#!/usr/bin/env python3
"""Exactly replay the archived n=28 LC witnesses and hard-family grids.

The checker deliberately avoids floating-point mode/mean decisions.  It
regenerates every structured tree from the parameter grids stored in the four
summary artifacts, recomputes its independence polynomial, and checks the
archived exact histograms and counts.  It also decodes and recomputes all 19
n=28 log-concavity witnesses retained by the Modal LC/near-miss summary.

The output is a compact deterministic certificate: individual rows are folded
into SHA-256 digests, while exact aggregate counts, order ranges, and
``d_leaf`` histograms remain human-readable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from graph6 import parse_graph6  # noqa: E402
from indpoly import independence_poly, is_log_concave, is_unimodal  # noqa: E402
from scripts.analyze_prufer_corpus import (  # noqa: E402
    GeneratedRecord,
    first_descent,
    graph6_from_adj,
    iter_family_records,
    structural_metrics,
)


FAMILY_ARTIFACTS = (
    "galvin_family_grid_m1-50_t1-50_nle500_summary.json",
    "bautista_ramos_family_grid_m1-8_t1-8_summary.json",
    "li_family_grid_m1-50_n1-50_summary.json",
    "li_star_family_grid_m1-50_n1-50_summary.json",
)
N28_ARTIFACT = "analysis_n28_modal_lc_nm.json"
DEFAULT_OUT = "mode_mean_artifact_replay_certificate.json"


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def histogram(values: Iterable[Any]) -> dict[str, int]:
    return {
        str(key): count
        for key, count in sorted(
            Counter(values).items(),
            key=lambda item: (str(type(item[0])), str(item[0])),
        )
    }


def exact_mean(poly: list[int]) -> Fraction:
    total = sum(poly)
    weighted = sum(k * coefficient for k, coefficient in enumerate(poly))
    return Fraction(weighted, total)


def mode_interval(poly: list[int]) -> tuple[int, int]:
    peak = max(poly)
    modes = [k for k, coefficient in enumerate(poly) if coefficient == peak]
    return modes[0], modes[-1]


def exact_lc_defects(poly: list[int]) -> list[dict[str, int]]:
    defects: list[dict[str, int]] = []
    for k in range(1, len(poly) - 1):
        square = poly[k] * poly[k]
        cross_product = poly[k - 1] * poly[k + 1]
        if cross_product > square:
            defects.append(
                {
                    "k": k,
                    "cross_product": cross_product,
                    "square": square,
                    "defect": cross_product - square,
                }
            )
    return defects


def exact_lc_ratio_position(poly: list[int]) -> tuple[Fraction, int]:
    worst = Fraction(0)
    worst_k = -1
    for k in range(1, len(poly) - 1):
        square = poly[k] * poly[k]
        if square == 0:
            continue
        ratio = Fraction(poly[k - 1] * poly[k + 1], square)
        if ratio > worst:
            worst = ratio
            worst_k = k
    return worst, worst_k


def exact_near_miss_position(poly: list[int]) -> int:
    descent = first_descent(poly)
    if descent is None:
        return -1
    worst = Fraction(0)
    worst_k = -1
    for k in range(descent, len(poly) - 1):
        if poly[k] == 0:
            continue
        ratio = Fraction(poly[k + 1], poly[k])
        if ratio > worst:
            worst = ratio
            worst_k = k
    return worst_k


def exact_target_score(record: GeneratedRecord, poly: list[int]) -> int | None:
    target_k = record.default_target_k
    if target_k is None:
        target_k = len(record.adj) // 2
    if target_k <= 0 or target_k >= len(poly) - 1:
        return None
    return poly[target_k - 1] * poly[target_k + 1] - poly[target_k] * poly[target_k]


def is_connected(adj: list[list[int]]) -> bool:
    if not adj:
        return True
    seen = {0}
    stack = [0]
    while stack:
        vertex = stack.pop()
        for neighbor in adj[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == len(adj)


def check_equal(
    mismatches: list[str],
    label: str,
    actual: Any,
    expected: Any,
) -> None:
    if actual != expected:
        mismatches.append(f"{label}: replayed {actual!r}, archived {expected!r}")


def replay_family(path: Path) -> tuple[dict[str, Any], list[str]]:
    archived = json.loads(path.read_text(encoding="utf-8"))
    source = archived["source"]
    records = iter_family_records(
        family=source["family"],
        m_values=source["m_values"],
        t_values=source["t_values"],
        max_n=source["max_n"],
    )
    limit = source.get("limit")
    if limit is not None:
        records = records[:limit]

    rows: list[dict[str, Any]] = []
    graph6_values: set[str] = set()
    for record in records:
        poly = independence_poly(len(record.adj), record.adj)
        mode_first, mode_last = mode_interval(poly)
        mean = exact_mean(poly)
        ceil_mean = math.ceil(mean)
        defects = exact_lc_defects(poly)
        defect_positions = {row["k"] for row in defects}
        near_miss_position = exact_near_miss_position(poly)
        structure = structural_metrics(record.adj)
        encoded = graph6_from_adj(record.adj)
        graph6_values.add(encoded)
        target_score = exact_target_score(record, poly)
        expected_breaks = tuple(record.expected_lc_breaks)

        rows.append(
            {
                "family": record.family,
                "params": record.params,
                "n": len(record.adj),
                "alpha": len(poly) - 1,
                "mode_first": mode_first,
                "mode_last": mode_last,
                "mean_numerator": mean.numerator,
                "mean_denominator": mean.denominator,
                "ceil_mean": ceil_mean,
                "ceil_mean_minus_mode": ceil_mean - mode_first,
                "first_descent": first_descent(poly),
                "unimodal": is_unimodal(poly),
                "log_concave": is_log_concave(poly),
                "lc_defects": defects,
                "near_miss_position": near_miss_position,
                "target_score": target_score,
                "expected_lc_breaks": list(expected_breaks),
                "expected_lc_breaks_missing": [
                    k for k in expected_breaks if k not in defect_positions
                ],
                "max_d_leaf": structure["max_d_leaf"],
                "d_leaf_le1": structure["d_leaf_le1"],
                "multi_leaf_hubs": structure["multi_leaf_hubs"],
                "mean_lt_n_over_3": 3 * mean.numerator < len(record.adj) * mean.denominator,
                "low_mode": mode_first <= len(record.adj) // 3 + 1,
                "poly_sha256": sha256_json(poly),
            }
        )

    lc_positions = [
        defect["k"]
        for row in rows
        for defect in row["lc_defects"]
    ]
    expected_rows = [row for row in rows if row["expected_lc_breaks"]]
    replayed_counts = {
        "non_unimodal": sum(not row["unimodal"] for row in rows),
        "non_log_concave": sum(not row["log_concave"] for row in rows),
        "low_mode": sum(row["low_mode"] for row in rows),
        "mean_lt_n_over_3": sum(row["mean_lt_n_over_3"] for row in rows),
        "mode_le_ceil_mu": sum(row["mode_first"] <= row["ceil_mean"] for row in rows),
        "lc_defects_all_after_mode": sum(
            bool(row["lc_defects"])
            and all(defect["k"] > row["mode_first"] for defect in row["lc_defects"])
            for row in rows
        ),
        "lc_defects_all_after_near_miss": sum(
            bool(row["lc_defects"])
            and row["near_miss_position"] != -1
            and all(
                defect["k"] > row["near_miss_position"]
                for defect in row["lc_defects"]
            )
            for row in rows
        ),
        "positive_target_score": sum(
            row["target_score"] is not None and row["target_score"] > 0
            for row in rows
        ),
        "expected_lc_breaks_all_hit": sum(
            not row["expected_lc_breaks_missing"] for row in expected_rows
        ),
        "expected_lc_breaks_any_missing": sum(
            bool(row["expected_lc_breaks_missing"]) for row in expected_rows
        ),
        "d_leaf_le1": sum(row["d_leaf_le1"] for row in rows),
        "multi_leaf_hub": sum(row["multi_leaf_hubs"] > 0 for row in rows),
        "unique_graph6": len(graph6_values),
        "unique_labelled_prufer": None,
        "unique_canonical_graph6": None,
    }
    replayed_exact_fields = {
        "processed": len(rows),
        "n_hist": histogram(row["n"] for row in rows),
        "alpha_hist": histogram(row["alpha"] for row in rows),
        "family_hist": histogram(row["family"] for row in rows),
        "mode_first_hist": histogram(row["mode_first"] for row in rows),
        "first_descent_hist": histogram(row["first_descent"] for row in rows),
        "lc_position_hist": histogram(lc_positions),
        "lc_defect_count_hist": histogram(len(row["lc_defects"]) for row in rows),
        "counts": replayed_counts,
    }

    mismatches: list[str] = []
    for field, actual in replayed_exact_fields.items():
        check_equal(mismatches, field, actual, archived[field])

    exact_means = [Fraction(row["mean_numerator"], row["mean_denominator"]) for row in rows]
    certificate = {
        "source_artifact": str(path.relative_to(REPO_ROOT)),
        "source_sha256": sha256_file(path),
        "family": source["family"],
        "parameter_rows": len(rows),
        "order": {
            "min": min(row["n"] for row in rows),
            "max": max(row["n"] for row in rows),
            "distinct": len({row["n"] for row in rows}),
            "greater_than_23": sum(row["n"] > 23 for row in rows),
            "histogram_sha256": sha256_json(replayed_exact_fields["n_hist"]),
        },
        "d_leaf_hist": histogram(row["max_d_leaf"] for row in rows),
        "unimodal_count": sum(row["unimodal"] for row in rows),
        "non_log_concave_count": sum(not row["log_concave"] for row in rows),
        "mode_le_ceil_mean_count": sum(
            row["mode_first"] <= row["ceil_mean"] for row in rows
        ),
        "mode_le_ceil_mean_failures": [
            {"params": row["params"], "mode": row["mode_first"], "ceil_mean": row["ceil_mean"]}
            for row in rows
            if row["mode_first"] > row["ceil_mean"]
        ],
        "ceil_mean_minus_mode_hist": histogram(
            row["ceil_mean_minus_mode"] for row in rows
        ),
        "mean": {
            "minimum": str(min(exact_means)),
            "maximum": str(max(exact_means)),
        },
        "archived_exact_fields_match": not mismatches,
        "rows_sha256": sha256_json(rows),
    }
    return certificate, mismatches


def replay_n28_lc_witnesses(path: Path) -> tuple[dict[str, Any], list[str]]:
    archived = json.loads(path.read_text(encoding="utf-8"))
    retained = archived["top_lc_failures"]
    mismatches: list[str] = []
    check_equal(mismatches, "lc_failure_count", len(retained), archived["lc_failure_count"])
    check_equal(mismatches, "retained witness count", len(retained), 19)

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for rank, item in enumerate(retained, start=1):
        encoded = item["graph6"]
        if encoded in seen:
            mismatches.append(f"duplicate retained graph6 witness at rank {rank}")
        seen.add(encoded)
        n, adj = parse_graph6(encoded.encode("ascii"))
        edge_count = sum(map(len, adj)) // 2
        poly = independence_poly(n, adj)
        defects = exact_lc_defects(poly)
        ratio, ratio_position = exact_lc_ratio_position(poly)
        mode_first, mode_last = mode_interval(poly)
        mean = exact_mean(poly)
        structure = structural_metrics(adj)

        if n != item["n"] or n != 28:
            mismatches.append(f"witness {rank}: decoded order {n}, archived {item['n']}")
        if edge_count != n - 1:
            mismatches.append(f"witness {rank}: graph has {edge_count} edges, expected {n - 1}")
        if not is_connected(adj):
            mismatches.append(f"witness {rank}: decoded graph is disconnected")
        if ratio_position != item["lc_pos"]:
            mismatches.append(
                f"witness {rank}: exact LC-ratio position {ratio_position}, "
                f"archived {item['lc_pos']}"
            )
        if not math.isclose(float(ratio), item["lc_ratio"], rel_tol=0.0, abs_tol=1e-15):
            mismatches.append(
                f"witness {rank}: exact LC ratio {ratio}, archived float {item['lc_ratio']}"
            )
        if not defects:
            mismatches.append(f"witness {rank}: recomputed polynomial is log-concave")

        rows.append(
            {
                "rank": rank,
                "graph6": encoded,
                "n": n,
                "edges": edge_count,
                "alpha": len(poly) - 1,
                "mode_first": mode_first,
                "mode_last": mode_last,
                "mean_numerator": mean.numerator,
                "mean_denominator": mean.denominator,
                "ceil_mean": math.ceil(mean),
                "unimodal": is_unimodal(poly),
                "lc_defects": defects,
                "worst_lc_ratio_numerator": ratio.numerator,
                "worst_lc_ratio_denominator": ratio.denominator,
                "worst_lc_ratio_position": ratio_position,
                "max_d_leaf": structure["max_d_leaf"],
                "poly_sha256": sha256_json(poly),
            }
        )

    if rows:
        exact_worst = max(
            rows,
            key=lambda row: Fraction(
                row["worst_lc_ratio_numerator"],
                row["worst_lc_ratio_denominator"],
            ),
        )
        archived_worst = archived["worst_lc_item"]
        check_equal(
            mismatches,
            "worst_lc_item.graph6",
            exact_worst["graph6"],
            archived_worst["graph6"],
        )
        check_equal(
            mismatches,
            "worst_lc_item.lc_pos",
            exact_worst["worst_lc_ratio_position"],
            archived_worst["lc_pos"],
        )

    certificate = {
        "source_artifact": str(path.relative_to(REPO_ROOT)),
        "source_sha256": sha256_file(path),
        "retained_witnesses": len(rows),
        "unique_graph6": len(seen),
        "order_hist": histogram(row["n"] for row in rows),
        "alpha_hist": histogram(row["alpha"] for row in rows),
        "d_leaf_hist": histogram(row["max_d_leaf"] for row in rows),
        "lc_defect_position_hist": histogram(
            defect["k"] for row in rows for defect in row["lc_defects"]
        ),
        "unimodal_count": sum(row["unimodal"] for row in rows),
        "mode_le_ceil_mean_count": sum(
            row["mode_first"] <= row["ceil_mean"] for row in rows
        ),
        "mode_le_ceil_mean_failures": [
            row["graph6"] for row in rows if row["mode_first"] > row["ceil_mean"]
        ],
        "archived_fields_match": not mismatches,
        "rows_sha256": sha256_json(sorted(rows, key=lambda row: row["graph6"])),
    }
    return certificate, mismatches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=REPO_ROOT / "results",
        help="Directory containing the archived JSON summaries",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "results" / DEFAULT_OUT,
        help="Deterministic JSON certificate to write",
    )
    args = parser.parse_args()

    mismatches: list[str] = []
    n28, n28_mismatches = replay_n28_lc_witnesses(args.results_dir / N28_ARTIFACT)
    mismatches.extend(f"{N28_ARTIFACT}: {message}" for message in n28_mismatches)

    families: list[dict[str, Any]] = []
    for filename in FAMILY_ARTIFACTS:
        certificate, family_mismatches = replay_family(args.results_dir / filename)
        families.append(certificate)
        mismatches.extend(f"{filename}: {message}" for message in family_mismatches)

    hard_family_rows = sum(item["parameter_rows"] for item in families)
    hard_family_rows_greater_than_23 = sum(
        item["order"]["greater_than_23"] for item in families
    )
    certificate = {
        "schema": "erdos993-mode-mean-artifact-replay-v1",
        "checker": "scripts/verify_mode_mean_artifacts.py",
        "n28_lc_witnesses": n28,
        "hard_families": families,
        "totals": {
            "n28_lc_witnesses": n28["retained_witnesses"],
            "hard_family_rows": hard_family_rows,
            "hard_family_rows_greater_than_23": hard_family_rows_greater_than_23,
            "all_replayed_rows": n28["retained_witnesses"] + hard_family_rows,
            "hard_family_rows_expected": 5823,
            "mode_le_ceil_mean_count": n28["mode_le_ceil_mean_count"]
            + sum(item["mode_le_ceil_mean_count"] for item in families),
            "unimodal_count": n28["unimodal_count"]
            + sum(item["unimodal_count"] for item in families),
        },
        "all_checks_passed": not mismatches and hard_family_rows == 5823,
    }

    if hard_family_rows != 5823:
        mismatches.append(f"hard-family row total: replayed {hard_family_rows}, expected 5823")
        certificate["all_checks_passed"] = False
    if mismatches:
        for message in mismatches:
            print(f"ERROR: {message}", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"replayed={certificate['totals']['all_replayed_rows']} "
        f"hard_family_rows={hard_family_rows} n28_lc_witnesses={n28['retained_witnesses']} "
        f"saved={args.out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
