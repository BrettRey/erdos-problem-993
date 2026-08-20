#!/usr/bin/env python3
"""Exact audits for notes/two_hub_B_logconcavity_proof_2026-08-20.md."""

from itertools import combinations_with_replacement
from math import comb
from pathlib import Path
import runpy


def conv(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def centrally_unimodal(a: list[int]) -> bool:
    if a != a[::-1]:
        return False
    return all(a[i] <= a[i + 1] for i in range((len(a) - 1) // 2))


def local_block(r: int, c: int) -> list[int]:
    a = [c * comb(r, j) for j in range(r + 1)]
    a[0] += 1
    a[-1] += 1
    return a


def even_connector_block(r: int, s: int, c: int, d: int) -> list[int]:
    """Rank coefficients of equation (2) in the connector note."""
    n = r + s + 1
    out = [c * d * comb(n, j) for j in range(n + 1)]
    for j in range(s + 1):
        out[j] += d * comb(s, j)
        out[r + 1 + j] += d * comb(s, j)
    for j in range(r + 1):
        out[j] += c * comb(r, j)
        out[s + 1 + j] += c * comb(r, j)
    out[0] += 1
    out[-1] += 1
    return out


def main() -> None:
    local_checks = quartet_checks = even_checks = 0
    for r in range(1, 61):
        for c in range(1, 7):
            ar = local_block(r, c)
            assert centrally_unimodal(ar)
            local_checks += 1

            # One active side, unmatched present side with q=0.
            assert centrally_unimodal(local_block(r + 1, c))

            # q=1: a positive scalar binomial row plus endpoint units.
            for scalar in range(1, 5):
                q1 = [scalar * c * comb(r, j) for j in range(r + 1)]
                q1[0] += 1
                q1[-1] += 1
                assert centrally_unimodal(q1)

            for s in range(1, 61):
                for d in range(1, 7):
                    four = conv(ar, local_block(s, d))
                    four[0] -= 1
                    four[-1] -= 1
                    assert min(four) >= 0
                    assert centrally_unimodal(four)
                    quartet_checks += 1

    for r in range(1, 31):
        for s in range(1, 31):
            for c in range(1, 5):
                for d in range(1, 5):
                    assert centrally_unimodal(even_connector_block(r, s, c, d))
                    even_checks += 1

    # Independent bounded check of the actual B coefficient sequence.
    packet = runpy.run_path(str(Path(__file__).parent /
        "gpt_attack/law_v_packet_2026-08-12/law_v_packet.py"))
    BC = packet["BC"]

    arms: list[tuple[int, ...]] = [()]
    for count in range(1, 4):
        for a in combinations_with_replacement(range(1, 11), count):
            if sum(a) <= 14:
                arms.append(a)

    pair_checks = 0
    for i, a in enumerate(arms):
        for b in arms[i:]:
            if not (a or b) or sum(a) + sum(b) > 18:
                continue
            B, _ = BC(a, b)
            assert all(B[k] * B[k] >= B[k - 1] * B[k + 1]
                       for k in range(1, len(B) - 1))
            pair_checks += 1

    print(f"PASS: {local_checks} local blocks and {quartet_checks} four-map cores")
    print(f"PASS: {even_checks} even-connector binomial cores (finite evidence)")
    print(f"PASS: B log-concave on {pair_checks} exact bounded arm pairs")


if __name__ == "__main__":
    main()
