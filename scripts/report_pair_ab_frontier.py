#!/usr/bin/env python3
"""Report nearest-miss frontier from pair-AB attachment search artifact."""

import argparse
import json
from fractions import Fraction


def F(nd):
    return Fraction(int(nd["num"]), int(nd["den"]))


def fmt_frac(nd):
    return f"{nd['num']}/{nd['den']}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--top", type=int, default=5)
    args = ap.parse_args()

    d = json.load(open(args.inp))
    t = d["totals"]

    print("checked_total", t["checked_total"])
    print("checked_A", t["checked_A"], "checked_B", t["checked_B"])
    print("kept_A", t["kept_A"], "kept_B", t["kept_B"])
    print("shared_keys", t["shared_keys"])
    print("witness_found", t["witness_found"])

    if d.get("witness") is not None:
        w = d["witness"]
        A = w["A"]; B = w["B"]
        print("\nWITNESS")
        print("A.g6", A["g6"], "A.N", A["N"], "A.m", A["m"], "A.lambda", fmt_frac(A["lambda"]), "A.rho", fmt_frac(A["rho"]))
        print("B.g6", B["g6"], "B.N", B["N"], "B.m", B["m"], "B.lambda", fmt_frac(B["lambda"]), "B.rho", fmt_frac(B["rho"]))
        # exact assertions
        assert A["m"] == B["m"], "m mismatch"
        assert F(A["lambda"]) == F(B["lambda"]), "lambda mismatch"
        assert F(A["rho"]) == F(B["rho"]), "rho mismatch"
        assert A["N"] != B["N"], "N not different"
        print("assertions_ok", True)

    frontier = d.get("nearest_miss_frontier", [])
    print("\nTOP_NEAREST_MISSES")
    for i, e in enumerate(frontier[: args.top]):
        print(
            i + 1,
            "m", e["m"],
            "lambda", fmt_frac(e["lambda"]),
            "NA", e["NA"],
            "NB", e["NB"],
            "deltaN", e["deltaN"],
            "rhoA", fmt_frac(e["rhoA"]),
            "rhoB", fmt_frac(e["rhoB"]),
            "|drho|", fmt_frac(e["abs_diff_rho"]),
        )


if __name__ == "__main__":
    main()
