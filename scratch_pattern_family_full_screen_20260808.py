#!/usr/bin/env python3
"""Numerical full-coefficient screen of Theorem 4.11 pattern trees.

For U_m=(T^v_{1,l}:S^w_{2,n})_k^(m), Theorem 4.11 gives

    I(U_m)=S_{2,l} T_{k,n}^m + x(1+2x)^l S_{2,n}^{km}.

The exact tail scanner only sees a fixed number of leading coefficients.
This program keeps every coefficient, scales after every FFT convolution,
and reports only valleys whose three coefficients exceed a relative floor.
Any reported candidate must be replayed in exact integer arithmetic.
No absence result from this floating-point screen is a proof.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass

import numpy as np
from scipy.signal import fftconvolve


@dataclass
class ScaledPoly:
    coeffs: np.ndarray
    log_scale: float = 0.0


def normalize(coeffs: np.ndarray, log_scale: float = 0.0) -> ScaledPoly:
    coeffs = np.asarray(coeffs, dtype=np.float64)
    # FFT roundoff can create tiny negative values at the support boundary.
    peak = float(np.max(np.abs(coeffs))) if coeffs.size else 0.0
    if peak == 0.0:
        return ScaledPoly(coeffs, -math.inf)
    coeffs[np.abs(coeffs) < peak * 1e-13] = 0.0
    coeffs[coeffs < 0.0] = 0.0
    peak = float(np.max(coeffs))
    return ScaledPoly(coeffs / peak, log_scale + math.log(peak))


def add(left: ScaledPoly, right: ScaledPoly) -> ScaledPoly:
    scale = max(left.log_scale, right.log_scale)
    size = max(left.coeffs.size, right.coeffs.size)
    out = np.zeros(size)
    out[: left.coeffs.size] += left.coeffs * math.exp(left.log_scale - scale)
    out[: right.coeffs.size] += right.coeffs * math.exp(right.log_scale - scale)
    return normalize(out, scale)


def mul(left: ScaledPoly, right: ScaledPoly) -> ScaledPoly:
    return normalize(
        fftconvolve(left.coeffs, right.coeffs),
        left.log_scale + right.log_scale,
    )


def power(base: ScaledPoly, exponent: int) -> ScaledPoly:
    result = ScaledPoly(np.array([1.0]))
    factor = base
    while exponent:
        if exponent & 1:
            result = mul(result, factor)
        exponent >>= 1
        if exponent:
            factor = mul(factor, factor)
    return result


def binomial(linear: tuple[float, float], exponent: int) -> ScaledPoly:
    return power(ScaledPoly(np.array(linear, dtype=np.float64)), exponent)


def shifted(poly: ScaledPoly, amount: int = 1) -> ScaledPoly:
    return ScaledPoly(np.pad(poly.coeffs, (amount, 0)), poly.log_scale)


def s2(t: int) -> ScaledPoly:
    return add(binomial((1.0, 2.0), t), shifted(binomial((1.0, 1.0), t)))


def tree_t(k: int, n: int, sn: ScaledPoly | None = None) -> ScaledPoly:
    sn = sn or s2(n)
    return add(power(sn, k), shifted(binomial((1.0, 2.0), k * n)))


def robust_valley(poly: ScaledPoly, relative_floor: float) -> tuple[int, float] | None:
    a = poly.coeffs
    floor = float(np.max(a)) * relative_floor
    decreasing = False
    best = None
    for index in range(a.size - 1):
        left = float(a[index])
        right = float(a[index + 1])
        if min(left, right) < floor:
            continue
        if right < left * (1.0 - 1e-10):
            decreasing = True
        elif decreasing and right > left * (1.0 + 1e-10):
            ratio = right / left
            if best is None or ratio > best[1]:
                best = (index, ratio)
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k-max", type=int, default=10)
    parser.add_argument("--ell-max", type=int, default=40)
    parser.add_argument("--n-max", type=int, default=40)
    parser.add_argument("--m-max", type=int, default=500)
    parser.add_argument("--degree-max", type=int, default=30_000)
    parser.add_argument("--relative-floor", type=float, default=1e-10)
    parser.add_argument("--progress-every", type=int, default=5_000)
    args = parser.parse_args()

    checked = 0
    closest = None
    for k in range(1, args.k_max + 1):
        for n in range(1, args.n_max + 1):
            stride = k * (n + 1)
            sn = s2(n)
            tk = tree_t(k, n, sn)
            sk = power(sn, k)
            tk_m = ScaledPoly(np.array([1.0]))
            sk_m = ScaledPoly(np.array([1.0]))
            max_m = min(args.m_max, (args.degree_max - 1) // stride)
            for m in range(1, max_m + 1):
                tk_m = mul(tk_m, tk)
                sk_m = mul(sk_m, sk)
                for ell in range(args.ell_max + 1):
                    degree = stride * m + ell + 1
                    if degree > args.degree_max:
                        break
                    first = mul(s2(ell), tk_m)
                    second = mul(shifted(binomial((1.0, 2.0), ell)), sk_m)
                    full = add(first, second)
                    checked += 1
                    hit = robust_valley(full, args.relative_floor)
                    if hit is not None:
                        index, ratio = hit
                        row = {
                            "k": k, "n": n, "ell": ell, "m": m,
                            "degree": degree, "index": index, "ratio": ratio,
                            "relative_floor": args.relative_floor,
                        }
                        print(json.dumps({"event": "candidate", **row}), flush=True)
                        return
                    if checked % args.progress_every == 0:
                        print(json.dumps({
                            "event": "progress", "checked": checked,
                            "parameters": [k, n, ell, m], "closest": closest,
                        }), flush=True)
    print(json.dumps({"event": "done", "checked": checked, "candidate": None}), flush=True)


if __name__ == "__main__":
    main()
