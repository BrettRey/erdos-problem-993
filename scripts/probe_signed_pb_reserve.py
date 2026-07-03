#!/usr/bin/env python3
"""Probe reserve for signed low-probability PB laws.

After splitting Bernoulli variables with p_i > 1/2, an arbitrary
Poisson-binomial law can be represented as a deterministic shift plus

    X - Y,

where X and Y are independent sums of Bernoulli variables with parameters at
most 1/2.  This script stress-tests the first post-descent reserve for that
signed law.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - environment guard
    raise SystemExit("This exploratory probe requires numpy.") from exc


Block = tuple[int, float]


def first_descent(seq: np.ndarray) -> int | None:
    for k in range(1, len(seq)):
        if seq[k] < seq[k - 1]:
            return k
    return None


def binomial_shape(n: int, p: float) -> np.ndarray:
    if n == 0:
        return np.array([1.0], dtype=np.float64)
    if p <= 0.0:
        out = np.zeros(n + 1, dtype=np.float64)
        out[0] = 1.0
        return out
    if p >= 1.0:
        out = np.zeros(n + 1, dtype=np.float64)
        out[n] = 1.0
        return out

    mode = min(max(int((n + 1) * p), 0), n)
    pmf = np.zeros(n + 1, dtype=np.float64)
    pmf[mode] = 1.0
    odds = p / (1.0 - p)
    for k in range(mode, n):
        pmf[k + 1] = pmf[k] * (n - k) / float(k + 1) * odds
    for k in range(mode, 0, -1):
        pmf[k - 1] = pmf[k] * k / float(n - k + 1) / odds
    return pmf / float(pmf.sum())


def grouped_pmf(blocks: list[Block]) -> np.ndarray:
    pmf = np.array([1.0], dtype=np.float64)
    for count, p in blocks:
        if count <= 0:
            continue
        pmf = np.convolve(pmf, binomial_shape(count, p))
        pmf /= float(pmf.sum())
    return pmf


def signed_metric(
    *,
    x_blocks: list[Block],
    y_blocks: list[Block],
    kind: str,
    details: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    x_pmf = grouped_pmf(x_blocks)
    y_pmf = grouped_pmf(y_blocks)
    signed = np.convolve(x_pmf, y_pmf[::-1])
    support_start = -(len(y_pmf) - 1)
    descent_index = first_descent(signed)
    if descent_index is None or descent_index >= len(signed) - 1 or signed[descent_index] <= 0:
        return None

    x_mean = sum(count * p for count, p in x_blocks)
    y_mean = sum(count * p for count, p in y_blocks)
    x_var = sum(count * p * (1.0 - p) for count, p in x_blocks)
    y_var = sum(count * p * (1.0 - p) for count, p in y_blocks)
    variance = x_var + y_var
    if variance <= 0:
        return None

    pressure = float(signed[descent_index + 1] / signed[descent_index])
    reserve = 1.0 - pressure
    row = {
        "kind": kind,
        "x_n": sum(count for count, _ in x_blocks),
        "y_n": sum(count for count, _ in y_blocks),
        "x_blocks": [[count, p] for count, p in x_blocks],
        "y_blocks": [[count, p] for count, p in y_blocks],
        "signed_mean": x_mean - y_mean,
        "variance": variance,
        "first_descent_index": descent_index,
        "first_descent_value": support_start + descent_index,
        "pressure": pressure,
        "reserve": reserve,
        "variance_times_reserve": variance * reserve,
        "quarter_ok": variance < 1.0 or reserve >= 1.0 / (4.0 * variance),
        "fifth_ok": variance < 1.0 or reserve >= 1.0 / (5.0 * variance),
    }
    if details:
        row.update(details)
    return row


def compact(row: dict[str, Any]) -> dict[str, Any]:
    keys = [
        "kind",
        "x_n",
        "y_n",
        "signed_mean",
        "variance",
        "first_descent_value",
        "pressure",
        "reserve",
        "variance_times_reserve",
        "lambda",
        "eta",
        "p",
        "q",
        "x_blocks",
        "y_blocks",
    ]
    return {key: row[key] for key in keys if key in row}


def random_grouped_blocks(rng: random.Random, n: int, max_groups: int) -> list[Block]:
    groups = rng.randint(1, min(max_groups, n))
    cuts = sorted(rng.sample(range(1, n), groups - 1)) if groups > 1 else []
    counts = [*cuts, n]
    counts = [counts[0], *[counts[i] - counts[i - 1] for i in range(1, len(counts))]]
    family = rng.choice(["uniform", "sparse", "middle", "edge"])
    probs: list[float] = []
    for _ in counts:
        if family == "uniform":
            probs.append(0.5 * rng.random())
        elif family == "sparse":
            probs.append(0.5 * rng.betavariate(0.35, 4.0))
        elif family == "middle":
            probs.append(0.5 * rng.betavariate(2.0, 2.0))
        else:
            probs.append(0.5 * rng.betavariate(0.4, 0.4))
    return sorted(zip(counts, probs), key=lambda item: item[1])


def add_random_rows(
    rows: list[dict[str, Any]],
    *,
    rng: random.Random,
    samples: int,
    max_groups: int,
) -> None:
    n_values = [5, 10, 20, 50, 100, 200, 300]
    for _ in range(samples):
        x_n = rng.choice(n_values)
        y_n = rng.choice(n_values)
        row = signed_metric(
            x_blocks=random_grouped_blocks(rng, x_n, max_groups),
            y_blocks=random_grouped_blocks(rng, y_n, max_groups),
            kind="random_grouped",
        )
        if row is not None:
            rows.append(row)


def add_two_binomial_rows(rows: list[dict[str, Any]], *, grid_size: int) -> None:
    n_pairs = [
        (10, 10),
        (20, 20),
        (50, 50),
        (100, 100),
        (200, 200),
        (300, 300),
        (50, 200),
        (200, 50),
        (100, 300),
        (300, 100),
    ]
    p_values = np.linspace(0.005, 0.5, grid_size)
    for x_n, y_n in n_pairs:
        for p_raw in p_values:
            for q_raw in p_values:
                p = float(p_raw)
                q = float(q_raw)
                row = signed_metric(
                    x_blocks=[(x_n, p)],
                    y_blocks=[(y_n, q)],
                    kind="two_binomial",
                    details={"p": p, "q": q},
                )
                if row is not None:
                    rows.append(row)


def add_finite_skellam_rows(rows: list[dict[str, Any]]) -> None:
    totals = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    approximants = [50, 100, 200, 500, 1000]
    split_values = [0.0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5]
    for total in totals:
        for split in split_values:
            lam = total * (1.0 - split)
            eta = total * split
            for base in approximants:
                x_n = max(1, int(math.ceil(base * max(lam, 0.05))))
                y_n = max(1, int(math.ceil(base * max(eta, 0.05))))
                p = min(0.5, lam / x_n) if lam > 0 else 0.0
                q = min(0.5, eta / y_n) if eta > 0 else 0.0
                row = signed_metric(
                    x_blocks=[(x_n, p)],
                    y_blocks=[(y_n, q)],
                    kind="finite_skellam",
                    details={"lambda": lam, "eta": eta, "p": p, "q": q},
                )
                if row is not None:
                    rows.append(row)


def best_by_cutoff(rows: list[dict[str, Any]], cutoffs: list[float]) -> list[dict[str, Any]]:
    out = []
    for cutoff in cutoffs:
        candidates = [row for row in rows if row["variance"] >= cutoff]
        best = min(candidates, key=lambda row: row["variance_times_reserve"], default=None)
        out.append({"variance_cutoff": cutoff, "best": compact(best) if best else None})
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=993)
    parser.add_argument("--random-samples", type=int, default=2000)
    parser.add_argument("--max-groups", type=int, default=6)
    parser.add_argument("--binomial-grid-size", type=int, default=31)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/signed_pb_reserve_probe_2026-07-03.json"),
    )
    args = parser.parse_args()
    rng = random.Random(args.seed)

    rows: list[dict[str, Any]] = []
    add_random_rows(
        rows,
        rng=rng,
        samples=args.random_samples,
        max_groups=args.max_groups,
    )
    add_two_binomial_rows(rows, grid_size=args.binomial_grid_size)
    add_finite_skellam_rows(rows)

    v_ge_1 = [row for row in rows if row["variance"] >= 1.0]
    quarter_failures = [row for row in rows if not row["quarter_ok"]]
    fifth_failures = [row for row in rows if not row["fifth_ok"]]
    summary = {
        "source": {
            "kind": "signed_low_probability_pb_reserve_probe",
            "seed": args.seed,
            "random_samples": args.random_samples,
            "max_groups": args.max_groups,
            "binomial_grid_size": args.binomial_grid_size,
        },
        "processed": len(rows),
        "counts": {
            "variance_ge_1": len(v_ge_1),
            "quarter_failures": len(quarter_failures),
            "fifth_failures": len(fifth_failures),
        },
        "best_by_variance_cutoff": best_by_cutoff(rows, [1, 2, 5, 10, 20, 50]),
        "top_smallest_v_reserve": [
            compact(row)
            for row in sorted(v_ge_1, key=lambda row: row["variance_times_reserve"])[
                : args.top
            ]
        ],
        "quarter_failures": [compact(row) for row in quarter_failures[: args.top]],
        "fifth_failures": [compact(row) for row in fifth_failures[: args.top]],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(
        f"Wrote {args.out}; processed={len(rows)}, "
        f"quarter_failures={len(quarter_failures)}, fifth_failures={len(fifth_failures)}",
        flush=True,
    )
    return 1 if quarter_failures or fifth_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
