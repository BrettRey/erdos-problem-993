#!/usr/bin/env python3
"""Exact counterexample to the independent-majority auxiliary conjecture.

This does *not* give a tree counterexample to Erdős 993.  It refutes the
proposed general lemma saying that every graph with a maximum independent set
of size m and complementary vertex-cover size p, where m >= 2p + 2, has a
unimodal independence polynomial.

Start with the Bhattacharyya--Kahn bipartite graph G(a,b):

* |V1| = b-a, |V2| = |V3| = a;
* V1--V2 is complete bipartite;
* V2--V3 is a perfect matching.

Its independence polynomial is

    (1+2x)^a - (1+x)^a + (1+x)^b.

Add s isolated vertices.  The resulting polynomial is

    (1+x)^s [(1+2x)^a - (1+x)^a + (1+x)^b].

For (a,b,s)=(105,167,45), the maximum independent set has size
b+s=212=2a+2, while the complementary minimum vertex cover has size a=105.
The exact coefficient sequence has a strict fall at rank 101 followed by a
strict rise at rank 102.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


A = 105
B = 167
S = 45


def convolve(left: list[int], right: list[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] += x * y
    return out


def binomial(n: int, scale: int = 1) -> list[int]:
    return [math.comb(n, k) * scale**k for k in range(n + 1)]


def polynomial(a: int, b: int, s: int) -> list[int]:
    first = convolve(binomial(s), binomial(a, 2))
    second = binomial(a + s)
    third = binomial(b + s)
    return [
        (first[k] if k < len(first) else 0)
        - (second[k] if k < len(second) else 0)
        + third[k]
        for k in range(len(third))
    ]


def first_rebound(poly: list[int]) -> tuple[int, int] | None:
    falling = False
    valley = -1
    for k in range(1, len(poly)):
        if poly[k] < poly[k - 1]:
            falling = True
            valley = k
        elif falling and poly[k] > poly[k - 1]:
            return valley, k
    return None


def main() -> None:
    poly = polynomial(A, B, S)
    rebound = first_rebound(poly)
    assert B > A
    alpha = B + S
    vertex_cover = A
    assert alpha == 2 * vertex_cover + 2
    assert rebound == (101, 102)
    assert poly[100] > poly[101] < poly[102]

    result = {
        "claim_scope": "exact counterexample to an auxiliary graph lemma; not a tree",
        "construction": "Bhattacharyya--Kahn G(a,b) plus s isolated vertices",
        "a": A,
        "b": B,
        "s": S,
        "order": A + B + S,
        "alpha": alpha,
        "minimum_vertex_cover": vertex_cover,
        "threshold_identity": f"{alpha} = 2*{vertex_cover} + 2",
        "first_rebound": {"valley_rank": rebound[0], "rise_rank": rebound[1]},
        "valley_window": {
            str(k): poly[k] for k in range(rebound[0] - 2, rebound[1] + 3)
        },
        "polynomial": poly,
    }
    out = Path("results/independent_majority_counterexample_20260807.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "polynomial"}, indent=2))


if __name__ == "__main__":
    main()
