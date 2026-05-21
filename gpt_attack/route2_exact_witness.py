"""Extract exact rational Route-2 data for selected witnesses.

Run from the repository root:

    python3 gpt_attack/route2_exact_witness.py
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import prod
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from graph6 import parse_graph6  # noqa: E402
from indpoly import independence_poly  # noqa: E402


WITNESSES = {
    "all_deg2_min_lam_endpoint_margin": {
        "g6": "S???????C?G?G?C?@??G??_?@?C@?B~_?",
        "leaf": 11,
        "support": 1,
    },
    "min_lam_endpoint_margin": {
        "g6": "S???????C?G?G?C?@??G??_?@??@?F~_?",
        "leaf": 10,
        "support": 0,
    },
    "min_exact_slack": {
        "g6": "T???????C?G?G?C?@??G??_?@??@???_B~o?",
        "leaf": 10,
        "support": 0,
    },
    "min_gap_deficit": {
        "g6": "V?????????O?O?G?A??O?@??A??A??@???O??A?F~o??",
        "leaf": 11,
        "support": 0,
    },
}


def remove_vertices(adj: list[list[int]], remove: set[int]) -> list[list[int]]:
    keep = [v for v in range(len(adj)) if v not in remove]
    idx = {v: i for i, v in enumerate(keep)}
    out = [[] for _ in keep]
    for v in keep:
        vv = idx[v]
        for u in adj[v]:
            if u in idx:
                out[vv].append(idx[u])
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


def fmt(fr: Fraction) -> str:
    return f"{fr} ~= {float(fr):.12g}"


def spider_signature(adj: list[list[int]]) -> str:
    deg = [len(nb) for nb in adj]
    hubs = [v for v, d in enumerate(deg) if d >= 3]
    if len(hubs) != 1:
        return "not a one-hub spider"
    h = hubs[0]
    lengths = []
    for nbr in adj[h]:
        prev = h
        cur = nbr
        length = 1
        while len(adj[cur]) == 2:
            nxt = adj[cur][0] if adj[cur][1] == prev else adj[cur][1]
            prev, cur = cur, nxt
            length += 1
        lengths.append(length)
    return "S(" + ",".join(map(str, sorted(lengths))) + ")"


def analyze(name: str, g6: str, leaf: int, support: int) -> None:
    n, adj = parse_graph6((g6 + "\n").encode("ascii"))
    poly_t = independence_poly(n, adj)
    m = mode_left(poly_t)
    lam = Fraction(poly_t[m - 1], poly_t[m])

    b_adj = remove_vertices(adj, {leaf, support})
    poly_b = independence_poly(len(b_adj), b_adj)
    tau = Fraction(poly_b[m - 2], poly_b[m - 1])

    mu_tau, var_tau = mean_var(poly_b, tau)
    mu_lam, var_lam = mean_var(poly_b, lam)
    threshold_half = Fraction(2 * m - 3, 2)
    threshold_exact = Fraction(m - 1) - lam / (1 + lam)
    deficit_tau = threshold_half - mu_tau
    gain = mu_lam - mu_tau
    gap = lam - tau
    tau_margin = var_tau / tau * gap - deficit_tau
    lam_margin = var_lam / lam * gap - deficit_tau

    print(f"\n{name}")
    print("-" * len(name))
    print(f"n = {n}")
    print(f"degree signature = {dict(sorted(Counter(map(len, adj)).items()))}")
    print(f"spider signature = {spider_signature(adj)}")
    print(f"mode m = {m}")
    print(f"lambda = {fmt(lam)}")
    print(f"tau = {fmt(tau)}")
    print(f"gap = {fmt(gap)}")
    print(f"mu_B(tau) = {fmt(mu_tau)}")
    print(f"mu_B(lambda) = {fmt(mu_lam)}")
    print(f"Var_B(tau) = {fmt(var_tau)}")
    print(f"Var_B(lambda) = {fmt(var_lam)}")
    print(f"deficit_tau = {fmt(deficit_tau)}")
    print(f"gain = {fmt(gain)}")
    print(f"route2_slack = {fmt(mu_lam - threshold_half)}")
    print(f"exact_slack = {fmt(mu_lam - threshold_exact)}")
    print(f"tau_endpoint_margin = {fmt(tau_margin)}")
    print(f"lambda_endpoint_margin = {fmt(lam_margin)}")
    print(f"I(T) = {poly_t}")
    print(f"I(B) = {poly_b}")


def main() -> None:
    for name, data in WITNESSES.items():
        analyze(name, data["g6"], data["leaf"], data["support"])


if __name__ == "__main__":
    main()
