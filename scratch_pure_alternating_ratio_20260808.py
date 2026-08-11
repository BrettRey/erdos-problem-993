#!/usr/bin/env python3
"""Explore exact ratio-order invariants in pure alternating rooted trees.

The grammar is

    N := leaf | N(P, ..., P),
    P := P(N, ..., N), with at least three N children.

For every rooted state v,

    E_v = product(F_c), J_v = product(E_c), F_v = E_v + x J_v.

This is a finite exact catalogue/kill-test for candidate induction invariants.
It is computational evidence, not a proof.
"""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from indpoly import _polyadd, _polymul, is_log_concave


Poly = tuple[int, ...]


def trim(poly: tuple[int, ...] | list[int]) -> Poly:
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def add(left: Poly, right: Poly) -> Poly:
    return trim(_polyadd(list(left), list(right)))


def mul(left: Poly, right: Poly) -> Poly:
    return trim(_polymul(list(left), list(right)))


def shift(poly: Poly, amount: int = 1) -> Poly:
    return (0,) * amount + poly


def partial_ratio_leq(left: Poly, right: Poly) -> bool:
    left = trim(left)
    right = trim(right)
    degree = max(len(left), len(right)) - 1
    a = left + (0,) * (degree + 1 - len(left))
    b = right + (0,) * (degree + 1 - len(right))
    delta_left = next(k for k, value in enumerate(a) if value)
    degree_right = max(k for k, value in enumerate(b) if value)
    if delta_left > degree_right + 1:
        return False
    return all(
        a[k] * b[k - 1] <= b[k] * a[k - 1]
        for k in range(1, degree + 1)
    )


def ratio_dominated(left: Poly, right: Poly) -> bool:
    return partial_ratio_leq(left, right) and partial_ratio_leq(
        right, shift(left),
    )


def relation(left: Poly, right: Poly, symbol: str) -> bool:
    if symbol == "pre":
        return partial_ratio_leq(left, right)
    if symbol == "rd":
        return ratio_dominated(left, right)
    if symbol == "coeff":
        degree = max(len(left), len(right))
        a = left + (0,) * (degree - len(left))
        b = right + (0,) * (degree - len(right))
        return all(x <= y for x, y in zip(a, b))
    raise ValueError(symbol)


@dataclass(frozen=True)
class State:
    kind: str
    order: int
    children: tuple["State", ...]
    e: Poly
    j: Poly
    f: Poly

    @property
    def signature(self) -> tuple[str, Poly, Poly, Poly]:
        return self.kind, self.e, self.j, self.f

    def short(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "order": self.order,
            "child_orders": [child.order for child in self.children],
            "E": self.e,
            "J": self.j,
            "F": self.f,
        }


def make_state(kind: str, children: tuple[State, ...]) -> State:
    e = (1,)
    j = (1,)
    for child in children:
        e = mul(e, child.f)
        j = mul(j, child.e)
    return State(
        kind=kind,
        order=1 + sum(child.order for child in children),
        children=children,
        e=e,
        j=j,
        f=add(e, shift(j)),
    )


def integer_partitions(total: int, minimum: int, count: int | None = None):
    """Yield nondecreasing positive partitions of total."""
    def rec(remaining: int, lower: int, values: tuple[int, ...]):
        if remaining == 0:
            if count is None or len(values) == count:
                yield values
            return
        if count is not None and len(values) >= count:
            return
        max_part = remaining
        if count is not None:
            slots = count - len(values)
            max_part = remaining - lower * (slots - 1)
        for part in range(lower, max_part + 1):
            yield from rec(remaining - part, part, values + (part,))
    yield from rec(total, minimum, ())


def build_catalogue(max_order: int, max_children: int) -> list[State]:
    by_kind_order: dict[tuple[str, int], list[State]] = defaultdict(list)
    leaf = make_state("N", ())
    by_kind_order[("N", 1)].append(leaf)
    all_states = [leaf]

    for order in range(2, max_order + 1):
        for kind, child_kind, min_children in (("N", "P", 1), ("P", "N", 3)):
            seen: set[tuple[str, Poly, Poly, Poly]] = set()
            child_total = order - 1
            for child_count in range(min_children, min(max_children, child_total) + 1):
                for partition in integer_partitions(child_total, 1, child_count):
                    pools = [by_kind_order[(child_kind, size)] for size in partition]
                    if any(not pool for pool in pools):
                        continue
                    for children in itertools.product(*pools):
                        # Canonicalize repeated equal-order child choices.
                        keys = [child.signature for child in children]
                        if any(
                            partition[i] == partition[i + 1] and keys[i] > keys[i + 1]
                            for i in range(len(children) - 1)
                        ):
                            continue
                        state = make_state(kind, children)
                        if state.signature in seen:
                            continue
                        seen.add(state.signature)
                        by_kind_order[(kind, order)].append(state)
                        all_states.append(state)
    return all_states


def candidate_polys(state: State) -> dict[str, Poly]:
    l = (1, 1)
    return {
        "J": state.j,
        "E": state.e,
        "F": state.f,
        "LJ": mul(l, state.j),
        "LE": mul(l, state.e),
        "LF": mul(l, state.f),
        "xJ": shift(state.j),
        "xE": shift(state.e),
        "xF": shift(state.f),
        "xLJ": shift(mul(l, state.j)),
        "xLE": shift(mul(l, state.e)),
        "xLF": shift(mul(l, state.f)),
    }


def audit(states: list[State]) -> dict[str, object]:
    counts: Counter[str] = Counter()
    first_failures: dict[str, dict[str, object]] = {}
    names = tuple(candidate_polys(states[0]))
    for state in states:
        counts[f"states_{state.kind}"] += 1
        counts[f"states_{state.kind}_order_{state.order}"] += 1
        polynomials = candidate_polys(state)
        for name, poly in polynomials.items():
            key = f"{state.kind}:{name}:lc"
            holds = is_log_concave(list(poly))
            counts[f"{key}:{'pass' if holds else 'fail'}"] += 1
            if not holds and key not in first_failures:
                first_failures[key] = state.short()
        for left, right in itertools.permutations(names, 2):
            for symbol in ("pre", "rd", "coeff"):
                key = f"{state.kind}:{left}:{symbol}:{right}"
                holds = relation(polynomials[left], polynomials[right], symbol)
                counts[f"{key}:{'pass' if holds else 'fail'}"] += 1
                if not holds and key not in first_failures:
                    first_failures[key] = state.short()

    universal = []
    for key in sorted(counts):
        if not key.endswith(":pass"):
            continue
        stem = key[:-5]
        if f"{stem}:fail" not in counts:
            universal.append({"relation": stem, "count": counts[key]})
    return {
        "state_count": len(states),
        "counts": dict(sorted(counts.items())),
        "universal_relations": universal,
        "first_failures": first_failures,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=25)
    parser.add_argument("--max-children", type=int, default=6)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    states = build_catalogue(args.max_order, args.max_children)
    result = {
        "scope": "exact finite grammar catalogue; evidence only",
        "max_order": args.max_order,
        "max_children": args.max_children,
        **audit(states),
    }
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
    concise = {
        "max_order": args.max_order,
        "max_children": args.max_children,
        "state_count": result["state_count"],
        "universal_relations": result["universal_relations"],
        "out": str(args.out) if args.out else None,
    }
    print(json.dumps(concise, indent=2))


if __name__ == "__main__":
    main()
