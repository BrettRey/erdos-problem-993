"""Route-2 endpoint scan for spider lanes.

This targets the all-leaf extremal pattern found on 2026-05-20:

    S(1, 2^7, 4)

Rather than enumerating all spider partitions, it scans the structured lane
S(1^j, 2^a, r) using the spider product formula:

    I_T(x) = prod_L P_L(x) + x prod_L P_{L-1}(x),

where P_L is the independence polynomial of a path on L vertices.

Run from the repository root:

    python3 gpt_attack/route2_spider_lane_scan.py --a-max 200 --r-max 40
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from math import inf
from pathlib import Path
from typing import Iterable

_POLY_POW_CACHE: dict[tuple[int, int], list[int]] = {}
_SPIDER_POLY_CACHE: dict[tuple[tuple[int, int], ...], list[int]] = {}


def poly_add(a: list[int], b: list[int]) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def poly_pow(p: list[int], e: int) -> list[int]:
    out = [1]
    base = p[:]
    while e:
        if e & 1:
            out = poly_mul(out, base)
        e >>= 1
        if e:
            base = poly_mul(base, base)
    return out


def trim(p: list[int]) -> list[int]:
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def path_polys(max_len: int) -> list[list[int]]:
    """P_L = independence polynomial of a path on L vertices."""
    polys = [[1], [1, 1]]
    for _ in range(2, max_len + 1):
        polys.append(poly_add(polys[-1], [0] + polys[-2]))
    return polys[: max_len + 1]


def spider_poly_from_counts(counts: dict[int, int], paths: list[list[int]]) -> list[int]:
    key = tuple(sorted((length, count) for length, count in counts.items() if count > 0))
    cached = _SPIDER_POLY_CACHE.get(key)
    if cached is not None:
        return cached

    hub_off = [1]
    hub_on = [1]
    for length, count in sorted(counts.items()):
        if count <= 0:
            continue
        off_key = (length, count)
        on_key = (length - 1, count)
        off_pow = _POLY_POW_CACHE.get(off_key)
        if off_pow is None:
            off_pow = poly_pow(paths[length], count)
            _POLY_POW_CACHE[off_key] = off_pow
        on_pow = _POLY_POW_CACHE.get(on_key)
        if on_pow is None:
            on_pow = poly_pow(paths[length - 1], count)
            _POLY_POW_CACHE[on_key] = on_pow
        hub_off = poly_mul(hub_off, off_pow)
        hub_on = poly_mul(hub_on, on_pow)
    out = poly_add(hub_off, [0] + hub_on)
    _SPIDER_POLY_CACHE[key] = out
    return out


def mode_left(poly: list[int]) -> int:
    mx = max(poly)
    return next(i for i, a in enumerate(poly) if a == mx)


def mean_var(poly: list[int], lam: Fraction) -> tuple[Fraction, Fraction]:
    weights = [Fraction(c) * lam**k for k, c in enumerate(poly)]
    z = sum(weights, Fraction(0))
    mu = sum(Fraction(k) * w for k, w in enumerate(weights)) / z
    second = sum(Fraction(k * k) * w for k, w in enumerate(weights)) / z
    return mu, second - mu * mu


def mean_var_float(poly: list[int], lam: float) -> tuple[float, float]:
    if lam <= 0.0:
        return 0.0, 0.0
    log_lam = math.log(lam)
    logs = [math.log(c) + k * log_lam for k, c in enumerate(poly) if c > 0]
    if not logs:
        return 0.0, 0.0
    shift = max(logs)
    z = 0.0
    mu = 0.0
    second = 0.0
    for k, c in enumerate(poly):
        if c <= 0:
            continue
        w = math.exp(math.log(c) + k * log_lam - shift)
        z += w
        mu += k * w
        second += k * k * w
    mu /= z
    return mu, second / z - mu * mu


def route2_for_removed_arm(
    counts: dict[int, int],
    remove_length: int,
    paths: list[list[int]],
) -> dict | None:
    if counts.get(remove_length, 0) <= 0 or remove_length < 2:
        return None

    poly_t = spider_poly_from_counts(counts, paths)
    m = mode_left(poly_t)
    if m < 2 or poly_t[m - 1] <= 0 or poly_t[m] <= 0:
        return None
    lam = Fraction(poly_t[m - 1], poly_t[m])

    b_counts = dict(counts)
    b_counts[remove_length] -= 1
    if b_counts[remove_length] == 0:
        del b_counts[remove_length]
    if remove_length > 2:
        b_counts[remove_length - 2] = b_counts.get(remove_length - 2, 0) + 1

    poly_b = spider_poly_from_counts(b_counts, paths)
    if m - 1 >= len(poly_b) or poly_b[m - 2] <= 0 or poly_b[m - 1] <= 0:
        return None

    tau = Fraction(poly_b[m - 2], poly_b[m - 1])
    mu_tau, var_tau = mean_var(poly_b, tau)
    mu_lam, var_lam = mean_var(poly_b, lam)
    threshold_half = Fraction(2 * m - 3, 2)
    threshold_exact = Fraction(m - 1) - lam / (1 + lam)
    deficit = threshold_half - mu_tau
    gap = lam - tau
    gain = mu_lam - mu_tau
    route2_slack = mu_lam - threshold_half
    exact_slack = mu_lam - threshold_exact

    lam_margin = None
    tau_margin = None
    if deficit > 0 and gap > 0:
        tau_margin = var_tau / tau * gap - deficit
        lam_margin = var_lam / lam * gap - deficit

    return {
        "counts": dict(sorted(counts.items())),
        "b_counts": dict(sorted(b_counts.items())),
        "remove_length": remove_length,
        "n": 1 + sum(length * count for length, count in counts.items()),
        "m": m,
        "lambda": lam,
        "tau": tau,
        "gap": gap,
        "mu_tau": mu_tau,
        "mu_lam": mu_lam,
        "var_tau": var_tau,
        "var_lam": var_lam,
        "deficit_tau": deficit,
        "gain": gain,
        "route2_slack": route2_slack,
        "exact_slack": exact_slack,
        "tau_margin": tau_margin,
        "lam_margin": lam_margin,
    }


def route2_float_for_removed_arm(
    counts: dict[int, int],
    remove_length: int,
    paths: list[list[int]],
) -> dict | None:
    if counts.get(remove_length, 0) <= 0 or remove_length < 2:
        return None

    poly_t = spider_poly_from_counts(counts, paths)
    m = mode_left(poly_t)
    if m < 2 or poly_t[m - 1] <= 0 or poly_t[m] <= 0:
        return None
    lam = poly_t[m - 1] / poly_t[m]

    b_counts = dict(counts)
    b_counts[remove_length] -= 1
    if b_counts[remove_length] == 0:
        del b_counts[remove_length]
    if remove_length > 2:
        b_counts[remove_length - 2] = b_counts.get(remove_length - 2, 0) + 1

    poly_b = spider_poly_from_counts(b_counts, paths)
    if m - 1 >= len(poly_b) or poly_b[m - 2] <= 0 or poly_b[m - 1] <= 0:
        return None

    tau = poly_b[m - 2] / poly_b[m - 1]
    mu_tau, var_tau = mean_var_float(poly_b, tau)
    mu_lam, var_lam = mean_var_float(poly_b, lam)
    threshold_half = m - 1.5
    threshold_exact = m - 1 - lam / (1 + lam)
    deficit = threshold_half - mu_tau
    gap = lam - tau
    gain = mu_lam - mu_tau
    route2_slack = mu_lam - threshold_half
    exact_slack = mu_lam - threshold_exact

    lam_margin = None
    tau_margin = None
    if deficit > 0.0 and gap > 0.0:
        tau_margin = var_tau / tau * gap - deficit
        lam_margin = var_lam / lam * gap - deficit

    return {
        "counts": dict(sorted(counts.items())),
        "b_counts": dict(sorted(b_counts.items())),
        "remove_length": remove_length,
        "n": 1 + sum(length * count for length, count in counts.items()),
        "m": m,
        "lambda": lam,
        "tau": tau,
        "gap": gap,
        "mu_tau": mu_tau,
        "mu_lam": mu_lam,
        "var_tau": var_tau,
        "var_lam": var_lam,
        "deficit_tau": deficit,
        "gain": gain,
        "route2_slack": route2_slack,
        "exact_slack": exact_slack,
        "tau_margin": tau_margin,
        "lam_margin": lam_margin,
    }


def as_jsonable(rec: dict) -> dict:
    out = {}
    for k, v in rec.items():
        if isinstance(v, Fraction):
            out[k] = {"num": v.numerator, "den": v.denominator, "float": float(v)}
        elif isinstance(v, dict):
            out[k] = {str(kk): vv for kk, vv in v.items()}
        else:
            out[k] = v
    return out


def fmt_counts(counts: dict[int, int]) -> str:
    parts = []
    for length, count in sorted(counts.items()):
        if count == 1:
            parts.append(str(length))
        else:
            parts.append(f"{length}^{count}")
    return "S(" + ",".join(parts) + ")"


def fmt_frac(v: Fraction | None) -> str:
    if v is None:
        return "None"
    return f"{float(v):.12g}"


def fmt_num(v) -> str:
    if v is None:
        return "None"
    return f"{float(v):.12g}"


def update_min(best: dict[str, dict | None], key: str, rec: dict) -> None:
    val = rec.get(key)
    if val is None:
        return
    if best.get(key) is None or val < best[key][key]:
        best[key] = rec


def scan_lane(j_values: Iterable[int], a_max: int, r_max: int, exact: bool) -> list[dict]:
    paths = path_polys(max(2, r_max))
    rows: list[dict] = []
    best: dict[str, dict | None] = {
        "lam_margin": None,
        "route2_slack": None,
        "exact_slack": None,
        "deficit_tau": None,
        "gap": None,
    }

    for j in j_values:
        for a in range(0, a_max + 1):
            for r in range(2, r_max + 1):
                counts: dict[int, int] = {}
                if j:
                    counts[1] = j
                if a:
                    counts[2] = a
                counts[r] = counts.get(r, 0) + 1
                if sum(counts.values()) < 3:
                    continue

                for remove_length in sorted(k for k in counts if k >= 2):
                    if exact:
                        rec = route2_for_removed_arm(counts, remove_length, paths)
                    else:
                        rec = route2_float_for_removed_arm(counts, remove_length, paths)
                    if rec is None:
                        continue
                    rows.append(rec)
                    for key in best:
                        if key == "deficit_tau":
                            val = rec[key]
                            if best[key] is None or val > best[key][key]:
                                best[key] = rec
                        else:
                            update_min(best, key, rec)

    print(f"checked {len(rows):,} spider lane records")
    for key, rec in best.items():
        if rec is None:
            continue
        val = rec[key]
        print(
            f"{key:>14}: {fmt_frac(val)}  "
            f"{fmt_counts(rec['counts'])} remove={rec['remove_length']} "
            f"m={rec['m']} n={rec['n']} route2={fmt_frac(rec['route2_slack'])}"
        )
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description="Route-2 scan for S(1^j,2^a,r) spider lanes.")
    ap.add_argument("--a-max", type=int, default=100)
    ap.add_argument("--r-max", type=int, default=40)
    ap.add_argument("--j-max", type=int, default=1)
    ap.add_argument("--float", action="store_true", help="Use float evaluation for wider prescans.")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    rows = scan_lane(range(args.j_max + 1), args.a_max, args.r_max, exact=not args.float)
    if args.out:
        out = {
            "params": {
                "a_max": args.a_max,
                "r_max": args.r_max,
                "j_max": args.j_max,
            },
            "rows": [as_jsonable(r) for r in rows],
        }
        path = Path(args.out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(out, indent=2))
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
