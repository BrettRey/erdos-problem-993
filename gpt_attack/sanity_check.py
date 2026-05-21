"""Small exact sanity checks for the GPT attack packet.

Run from the repository root:

    python3 gpt_attack/sanity_check.py

This is not a proof checker. It guards against two cheap failure modes:
false global mean/mode claims, and over-interpreting sparse LC failures as
good density targets for DRC-style arguments.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from graph6 import parse_graph6  # noqa: E402
from indpoly import independence_poly  # noqa: E402


def star_poly(d: int) -> list[int]:
    """Independence polynomial coefficients for K_{1,d}."""
    return [1, d + 1] + [comb(d, k) for k in range(2, d + 1)]


def modes(poly: list[int]) -> list[int]:
    m = max(poly)
    return [i for i, a in enumerate(poly) if a == m]


def mean(poly: list[int]) -> Fraction:
    return Fraction(sum(k * a for k, a in enumerate(poly)), sum(poly))


def check_star_guard() -> None:
    print("Star guard: global mode <= floor(n/3)+1 and global mu < n/3 are false.")
    for d in (3, 10, 20):
        poly = star_poly(d)
        n = d + 1
        mu = mean(poly)
        print(
            f"  K_1,{d}: n={n}, modes={modes(poly)}, "
            f"floor(n/3)+1={n // 3 + 1}, mu={mu}, n/3={Fraction(n, 3)}"
        )


def poly_from_item(item: dict) -> list[int]:
    if item.get("poly") is not None:
        return item["poly"]
    n, adj = parse_graph6((item["graph6"] + "\n").encode("ascii"))
    if n != item["n"]:
        raise ValueError(f"graph6 n mismatch: {n} != {item['n']}")
    return independence_poly(n, adj)


def density_summary(path: Path) -> None:
    if not path.exists():
        return

    data = json.loads(path.read_text())
    print(f"\nDensity probe from {path}:")

    lc_densities = []
    for item in data.get("top_lc_failures", []):
        poly = poly_from_item(item)
        k = item["lc_pos"]
        lc_densities.append(Fraction(poly[k], comb(item["n"], k)))
    if lc_densities:
        print(
            "  LC failures at failure index: "
            f"min={min(lc_densities)} ({float(min(lc_densities)):.6g}), "
            f"max={max(lc_densities)} ({float(max(lc_densities)):.6g})"
        )

    nm_densities = []
    for item in data.get("top_near_misses", [])[:50]:
        poly = poly_from_item(item)
        k = item["nm_pos"]
        nm_densities.append(Fraction(poly[k], comb(item["n"], k)))
    if nm_densities:
        print(
            "  Top near misses at near-miss index: "
            f"min={min(nm_densities)} ({float(min(nm_densities)):.6g}), "
            f"max={max(nm_densities)} ({float(max(nm_densities)):.6g})"
        )


def main() -> None:
    check_star_guard()
    density_summary(ROOT / "results" / "analysis_n28_modal_lc_nm.json")


if __name__ == "__main__":
    main()
