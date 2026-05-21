"""Fixed-r hub-off certificate using a C++ reserve checker.

The SymPy certificate is convenient for low-degree margin and fugacity checks,
but the shifted reserve polynomial becomes dense and high-degree for larger
``r``.  This driver keeps the symbolic setup in Python and sends the reserve
coefficient arithmetic to ``fixed_r_huboff_reserve_cpp.cpp``.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from fractions import Fraction
from math import ceil
from pathlib import Path

import sympy as sp

from fixed_r_huboff_certificate import (
    normalized_coeff_ratio,
    poly_pair_common_integer_coeffs,
    shifted_linear_int_poly,
    shifted_poly,
    shifted_positive,
    stabilized_shifts,
    t,
)
from route2_spider_lane_scan import path_polys


def compile_helper() -> Path:
    source = Path(__file__).with_name("fixed_r_huboff_reserve_cpp.cpp")
    binary = Path(tempfile.gettempdir()) / "fixed_r_huboff_reserve_cpp"
    if not binary.exists() or source.stat().st_mtime > binary.stat().st_mtime:
        include_flags = []
        lib_flags = []
        for prefix in [Path("/opt/homebrew"), Path("/usr/local")]:
            if (prefix / "include" / "gmpxx.h").exists():
                include_flags = [f"-I{prefix / 'include'}"]
                lib_flags = [f"-L{prefix / 'lib'}"]
                break
        subprocess.run(
            [
                "clang++",
                "-O2",
                "-std=c++17",
                *include_flags,
                str(source),
                "-o",
                str(binary),
                *lib_flags,
                "-lgmpxx",
                "-lgmp",
            ],
            check=True,
        )
    return binary


def serialize_poly(poly: list[int]) -> str:
    return f"{len(poly)} " + " ".join(str(coeff) for coeff in poly)


def reserve_check_cpp(
    binary: Path,
    r: int,
    reserve_denom: int,
    l_coeffs: list[int],
    z_coeffs: list[int],
    a_poly: list[int],
    m_poly: list[int],
    reserve_threads: int,
) -> str:
    payload = "\n".join(
        [
            f"{r} {reserve_denom}",
            serialize_poly(l_coeffs),
            serialize_poly(z_coeffs),
            serialize_poly(a_poly),
            serialize_poly(m_poly),
            "",
        ]
    )
    env = os.environ.copy()
    env["FIXED_R_RESERVE_THREADS"] = str(reserve_threads)
    proc = subprocess.run(
        [str(binary)],
        input=payload,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    if proc.returncode not in (0, 1):
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip())
    return proc.stdout.strip()


def check_r(
    r: int,
    threshold: int,
    margin_denom: int,
    reserve_denom: int,
    reserve_threads: int,
) -> bool:
    binary = compile_helper()
    paths = path_polys(max(2, r))
    pr = paths[r]
    shifts = stabilized_shifts(r, threshold, paths)
    print(f"r={r}, threshold a>={threshold}")
    print(f"stabilized shifts D_q=3m-2a: {shifts}")

    all_ok = True
    for residue, shift in shifts.items():
        threshold_t = ceil((threshold - residue) / 3)
        a = 3 * t + residue
        m = sp.simplify((2 * a + shift) / 3)
        f_minus = normalized_coeff_ratio(pr, a, m, -1)
        f_zero = normalized_coeff_ratio(pr, a, m, 0)
        f_plus = normalized_coeff_ratio(pr, a, m, 1)
        left_margin = sp.cancel((f_zero - f_minus) / f_zero)
        right_margin = sp.cancel((f_zero - f_plus) / f_zero)
        lam0 = sp.cancel(f_minus / f_zero)
        checks = {
            f"left_margin_ge_1_over_{margin_denom}a": left_margin
            - Fraction(1, margin_denom) / a,
            f"right_margin_ge_1_over_{margin_denom}a": right_margin
            - Fraction(1, margin_denom) / a,
            "lambda0_ge_3_over_4": lam0 - sp.Rational(3, 4),
            "lambda0_ge_1_over_2": lam0 - sp.Rational(1, 2),
            "lambda0_le_2": sp.Rational(2) - lam0,
        }
        print(f"  residue q={residue}, t>={threshold_t}, m={m}")
        for name, expr in checks.items():
            ok = shifted_positive(expr, threshold_t)
            all_ok = all_ok and ok
            print(f"    {name}: {ok}")

        lam_num, lam_den = sp.fraction(sp.cancel(lam0))
        l_poly = shifted_poly(sp.Poly(lam_num, t), threshold_t)
        z_poly = shifted_poly(sp.Poly(lam_den, t), threshold_t)
        l_coeffs, z_coeffs = poly_pair_common_integer_coeffs(l_poly, z_poly)
        a_poly = shifted_linear_int_poly(a.subs(t, t + threshold_t))
        m_poly = shifted_linear_int_poly(m.subs(t, t + threshold_t))
        reserve_output = reserve_check_cpp(
            binary,
            r,
            reserve_denom,
            l_coeffs,
            z_coeffs,
            a_poly,
            m_poly,
            reserve_threads,
        )
        reserve_ok = "reserve_ok=true" in reserve_output
        all_ok = all_ok and reserve_ok
        print(f"    reserve_ge_1_over_{reserve_denom}a: {reserve_ok} ({reserve_output})")

    return all_ok


def main() -> None:
    ap = argparse.ArgumentParser(description="Fixed-r hub-off certificate with C++ reserve.")
    ap.add_argument("--r", type=int, required=True)
    ap.add_argument("--threshold", type=int, default=200)
    ap.add_argument("--margin-denom", type=int, default=1000)
    ap.add_argument("--reserve-denom", type=int, default=1000)
    ap.add_argument("--reserve-threads", type=int, default=1)
    args = ap.parse_args()

    ok = check_r(
        args.r,
        args.threshold,
        args.margin_denom,
        args.reserve_denom,
        args.reserve_threads,
    )
    assert ok


if __name__ == "__main__":
    main()
