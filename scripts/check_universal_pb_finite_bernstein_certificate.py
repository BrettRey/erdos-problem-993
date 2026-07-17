#!/usr/bin/env python3
"""Independently check the full finite Bernstein certificate.

This checker intentionally uses only the Python standard library.  It neither
imports the SymPy generator nor evaluates formula strings from the certificate.
Instead it implements the two window formulas directly with exact
``fractions.Fraction`` polynomial arithmetic, checks their certified
numerators by cross multiplication, and reconstructs every Bernstein
identity.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, gcd
from pathlib import Path
import re
import sys
from typing import Any, Iterable


DEFAULT_CERTIFICATE = Path(
    "results/universal_pb_finite_bernstein_full_certificate_2026-07-16.json"
)

LEGACY_CERTIFICATE_PATH = (
    "results/universal_pb_finite_bernstein_certificate_2026-07-10.json"
)
LEGACY_CERTIFICATE_SHA256 = (
    "6b91554d9ab1f43151e36c94c5c8c427c7bb057130b7f39d233b14c7ab3860c6"
)

SOURCE_FORMULA_IDS = {
    "asymmetric_window_v1": (
        "u_r=((H-1)/H)^r*prod_{j=1}^{r-1}(1-j/(H+1)); "
        "v_r=prod_{j=1}^{r}(1-j/H); "
        "A=(1+sum(u_r+v_r))*sum(r^2*(u_r+v_r))"
        "-(sum(r*(u_r-v_r)))^2"
    ),
    "symmetric_window_v1": (
        "b_r=prod_{j=1}^{r}(1-j/H); "
        "A=(1+2*sum(b_r))*(2*sum(r^2*b_r))"
    ),
    "scalar_target_v1": "Q=(H+1)*(3*H+4)/4",
}

RATIONAL_PATTERN = re.compile(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?\Z")

Polynomial = tuple[Fraction, ...]


class CertificateError(Exception):
    """Raised when a certificate check fails."""


def fail(message: str) -> None:
    raise CertificateError(message)


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def json_digest(value: Any) -> str:
    return sha256(canonical_json(value).encode("ascii")).hexdigest()


def exact_string(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def parse_rational(raw: Any, context: str) -> Fraction:
    if not isinstance(raw, str) or RATIONAL_PATTERN.fullmatch(raw) is None:
        fail(f"{context}: expected a canonical rational string")
    if "/" in raw:
        numerator_text, denominator_text = raw.split("/", 1)
        numerator = int(numerator_text)
        denominator = int(denominator_text)
        if gcd(abs(numerator), denominator) != 1:
            fail(f"{context}: rational is not reduced")
        value = Fraction(numerator, denominator)
    else:
        value = Fraction(int(raw), 1)
    if exact_string(value) != raw:
        fail(f"{context}: rational is not canonically encoded")
    return value


def require_exact_keys(value: Any, keys: set[str], context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        fail(f"{context}: expected an object")
    actual = set(value)
    if actual != keys:
        missing = sorted(keys - actual)
        extra = sorted(actual - keys)
        fail(f"{context}: key mismatch; missing={missing}, extra={extra}")
    return value


def require_int(value: Any, context: str) -> int:
    if type(value) is not int:
        fail(f"{context}: expected an integer")
    return value


def require_string(value: Any, context: str) -> str:
    if not isinstance(value, str):
        fail(f"{context}: expected a string")
    return value


def p_trim(coefficients: Iterable[Fraction]) -> Polynomial:
    result = list(coefficients)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return tuple(result or [Fraction(0)])


def p_const(value: int | Fraction) -> Polynomial:
    return (Fraction(value),)


def p_add(left: Polynomial, right: Polynomial) -> Polynomial:
    size = max(len(left), len(right))
    return p_trim(
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    )


def p_neg(value: Polynomial) -> Polynomial:
    return tuple(-coefficient for coefficient in value)


def p_sub(left: Polynomial, right: Polynomial) -> Polynomial:
    return p_add(left, p_neg(right))


def p_scale(value: Polynomial, scalar: int | Fraction) -> Polynomial:
    return p_trim(Fraction(scalar) * coefficient for coefficient in value)


def p_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] += left_value * right_value
    return p_trim(result)


def p_pow(value: Polynomial, exponent: int) -> Polynomial:
    if exponent < 0:
        fail("internal error: negative polynomial exponent")
    result = p_const(1)
    base = value
    power = exponent
    while power:
        if power & 1:
            result = p_mul(result, base)
        base = p_mul(base, base)
        power //= 2
    return result


def p_product(values: Iterable[Polynomial]) -> Polynomial:
    result = p_const(1)
    for value in values:
        result = p_mul(result, value)
    return result


def h_minus(h: Polynomial, constant: int) -> Polynomial:
    return p_sub(h, p_const(constant))


def expected_symmetric_numerator(left_endpoint: int) -> Polynomial:
    """Return 4 H^(2m) (S_m T_m - Q(H)) after H=m+t."""
    m = left_endpoint
    h: Polynomial = (Fraction(m), Fraction(1))
    h_to_m = p_pow(h, m)
    total_mass_numerator = h_to_m
    second_moment_numerator = p_const(0)
    b_numerator = p_const(1)
    for index in range(1, m + 1):
        b_numerator = p_mul(b_numerator, h_minus(h, index))
        lifted = p_mul(b_numerator, p_pow(h, m - index))
        total_mass_numerator = p_add(
            total_mass_numerator, p_scale(lifted, 2)
        )
        second_moment_numerator = p_add(
            second_moment_numerator, p_scale(lifted, 2 * index * index)
        )
    window_numerator = p_scale(
        p_mul(total_mass_numerator, second_moment_numerator), 4
    )
    target_numerator = p_mul(
        p_pow(h, 2 * m),
        p_mul(p_add(h, p_const(1)), p_add(p_scale(h, 3), p_const(4))),
    )
    return p_sub(window_numerator, target_numerator)


def asymmetric_cross_products(
    certified_numerator: Polynomial,
) -> tuple[Polynomial, Polynomial]:
    """Return the two sides of P/E = A-Q for the [3,4] cell."""
    h: Polynomial = (Fraction(3), Fraction(1))
    h_plus_one = p_add(h, p_const(1))
    common_denominator = p_mul(p_pow(h, 3), p_pow(h_plus_one, 2))
    total_mass_numerator = common_denominator
    first_moment_numerator = p_const(0)
    second_moment_numerator = p_const(0)

    for index in range(1, 4):
        u_numerator = p_mul(
            p_pow(h_minus(h, 1), index),
            p_product(h_minus(h, j - 1) for j in range(1, index)),
        )
        u_lifted = p_mul(
            u_numerator,
            p_mul(p_pow(h, 3 - index), p_pow(h_plus_one, 3 - index)),
        )

        v_numerator = p_product(h_minus(h, j) for j in range(1, index + 1))
        v_lifted = p_mul(
            v_numerator,
            p_mul(p_pow(h, 3 - index), p_pow(h_plus_one, 2)),
        )

        total_mass_numerator = p_add(
            total_mass_numerator, p_add(u_lifted, v_lifted)
        )
        first_moment_numerator = p_add(
            first_moment_numerator, p_scale(p_sub(u_lifted, v_lifted), index)
        )
        second_moment_numerator = p_add(
            second_moment_numerator,
            p_scale(p_add(u_lifted, v_lifted), index * index),
        )

    window_numerator = p_sub(
        p_mul(total_mass_numerator, second_moment_numerator),
        p_pow(first_moment_numerator, 2),
    )
    raw_numerator = p_sub(
        p_scale(window_numerator, 4),
        p_mul(
            p_pow(common_denominator, 2),
            p_mul(h_plus_one, p_add(p_scale(h, 3), p_const(4))),
        ),
    )
    raw_denominator = p_scale(p_pow(common_denominator, 2), 4)
    certified_denominator = p_scale(
        p_mul(p_pow(h, 5), p_pow(h_plus_one, 3)), 4
    )
    return (
        p_mul(certified_numerator, raw_denominator),
        p_mul(raw_numerator, certified_denominator),
    )


def bernstein_reconstruction(coefficients: list[Fraction]) -> Polynomial:
    degree = len(coefficients) - 1
    t: Polynomial = (Fraction(0), Fraction(1))
    one_minus_t: Polynomial = (Fraction(1), Fraction(-1))
    result = p_const(0)
    for index, coefficient in enumerate(coefficients):
        basis = p_mul(
            p_pow(t, index), p_pow(one_minus_t, degree - index)
        )
        result = p_add(result, p_scale(basis, coefficient * comb(degree, index)))
    return result


def expected_cell_roster() -> list[dict[str, Any]]:
    cells = [
        {
            "cell_id": "asymmetric_3_4",
            "formula_id": "asymmetric_window_v1",
            "window_count": 3,
            "left_endpoint": 3,
            "degree": 10,
            "expected_denominator": {
                "constant": 4,
                "factors": [["H", 5], ["H+1", 3]],
            },
        }
    ]
    cells.extend(
        {
            "cell_id": f"symmetric_{left}_{left + 1}",
            "formula_id": "symmetric_window_v1",
            "window_count": left,
            "left_endpoint": left,
            "degree": 2 * left + 2,
            "expected_denominator": {
                "constant": 4,
                "factors": [["H", 2 * left]],
            },
        }
        for left in range(4, 16)
    )
    return cells


CELL_KEYS = {
    "cell_id",
    "formula_id",
    "window_count",
    "interval",
    "left_endpoint",
    "degree",
    "expected_denominator",
    "numerator_power_coefficients_low_to_high",
    "bernstein_coefficients_low_to_high",
    "coefficient_count",
    "bernstein_vector_sha256",
    "cell_payload_sha256",
}


def validate_cell(record: Any, expected: dict[str, Any], index: int) -> int:
    context = f"payload.cells[{index}]"
    cell = require_exact_keys(record, CELL_KEYS, context)
    for field in ("cell_id", "formula_id"):
        if require_string(cell[field], f"{context}.{field}") != expected[field]:
            fail(f"{context}.{field}: unexpected value")
    for field in ("window_count", "left_endpoint", "degree"):
        if require_int(cell[field], f"{context}.{field}") != expected[field]:
            fail(f"{context}.{field}: unexpected value")
    left_endpoint = expected["left_endpoint"]
    if cell["interval"] != [left_endpoint, left_endpoint + 1]:
        fail(f"{context}.interval: unexpected value")
    if cell["expected_denominator"] != expected["expected_denominator"]:
        fail(f"{context}.expected_denominator: unexpected value")

    degree = expected["degree"]
    power_raw = cell["numerator_power_coefficients_low_to_high"]
    bernstein_raw = cell["bernstein_coefficients_low_to_high"]
    if not isinstance(power_raw, list) or not isinstance(bernstein_raw, list):
        fail(f"{context}: coefficient vectors must be arrays")
    if len(power_raw) != degree + 1 or len(bernstein_raw) != degree + 1:
        fail(f"{context}: coefficient-vector length does not match degree")
    coefficient_count = require_int(
        cell["coefficient_count"], f"{context}.coefficient_count"
    )
    if coefficient_count != degree + 1:
        fail(f"{context}.coefficient_count: unexpected value")

    power = [
        parse_rational(value, f"{context}.numerator_power[{position}]")
        for position, value in enumerate(power_raw)
    ]
    bernstein = [
        parse_rational(value, f"{context}.bernstein[{position}]")
        for position, value in enumerate(bernstein_raw)
    ]
    certified_numerator = p_trim(power)
    if len(certified_numerator) != degree + 1:
        fail(f"{context}: certified numerator does not have the stated degree")
    if any(coefficient <= 0 for coefficient in bernstein):
        fail(f"{context}: Bernstein coefficient is not strictly positive")

    expected_vector_digest = sha256(
        ",".join(bernstein_raw).encode("ascii")
    ).hexdigest()
    if cell["bernstein_vector_sha256"] != expected_vector_digest:
        fail(f"{context}: Bernstein-vector digest mismatch")
    cell_without_digest = dict(cell)
    claimed_cell_digest = cell_without_digest.pop("cell_payload_sha256")
    if claimed_cell_digest != json_digest(cell_without_digest):
        fail(f"{context}: cell-payload digest mismatch")

    reconstructed = bernstein_reconstruction(bernstein)
    if reconstructed != certified_numerator:
        fail(f"{context}: Bernstein reconstruction does not match numerator")

    if expected["formula_id"] == "asymmetric_window_v1":
        left_cross_product, right_cross_product = asymmetric_cross_products(
            certified_numerator
        )
        if left_cross_product != right_cross_product:
            fail(f"{context}: numerator does not match asymmetric source formula")
    else:
        source_numerator = expected_symmetric_numerator(left_endpoint)
        if certified_numerator != source_numerator:
            fail(f"{context}: numerator does not match symmetric source formula")

    return coefficient_count


TOP_LEVEL_KEYS = {
    "schema",
    "certificate_date",
    "claim",
    "rational_encoding",
    "variable",
    "domain",
    "basis",
    "coefficient_order",
    "substitution",
    "source_formula_ids",
    "expands_summary_artifact",
    "generator",
    "payload",
    "summary",
    "payload_sha256",
}


def validate_certificate(value: Any) -> dict[str, Any]:
    certificate = require_exact_keys(value, TOP_LEVEL_KEYS, "certificate")
    fixed_metadata = {
        "schema": "universal-pb-finite-bernstein/v2",
        "certificate_date": "2026-07-16",
        "claim": "positivity of the thirteen compact scalar-cell numerators",
        "rational_encoding": "reduced decimal p or p/q with q>0",
        "variable": "t",
        "domain": ["0", "1"],
        "basis": "binomial(d,i)*t^i*(1-t)^(d-i)",
        "coefficient_order": "ascending index 0 through degree",
        "substitution": "H=left_endpoint+t",
        "source_formula_ids": SOURCE_FORMULA_IDS,
        "expands_summary_artifact": {
            "path": LEGACY_CERTIFICATE_PATH,
            "sha256": LEGACY_CERTIFICATE_SHA256,
        },
        "generator": {
            "path": "scripts/verify_universal_pb_finite_bernstein.py",
            "sympy_version": "1.14.0",
        },
    }
    for field, expected in fixed_metadata.items():
        if certificate[field] != expected:
            fail(f"certificate.{field}: unexpected value")

    payload = require_exact_keys(certificate["payload"], {"cells"}, "payload")
    cells = payload["cells"]
    if not isinstance(cells, list):
        fail("payload.cells: expected an array")
    roster = expected_cell_roster()
    if len(cells) != len(roster):
        fail("payload.cells: unexpected number of cells")

    total_coefficients = sum(
        validate_cell(cell, expected, index)
        for index, (cell, expected) in enumerate(zip(cells, roster, strict=True))
    )
    expected_summary = {
        "total_cells": 13,
        "total_coefficients": 275,
        "all_coefficients_strictly_positive": True,
        "all_source_identities_exact": True,
        "all_bernstein_identities_exact": True,
    }
    if total_coefficients != expected_summary["total_coefficients"]:
        fail("payload.cells: unexpected total coefficient count")
    if certificate["summary"] != expected_summary:
        fail("certificate.summary: unexpected value")
    expected_payload_digest = json_digest(payload)
    if certificate["payload_sha256"] != expected_payload_digest:
        fail("certificate.payload_sha256: digest mismatch")
    return {
        "status": "passed",
        "cells": len(cells),
        "coefficients": total_coefficients,
        "payload_sha256": expected_payload_digest,
    }


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_float(raw: str) -> Any:
    fail(f"JSON floating-point number is not permitted: {raw}")


def load_certificate(path: str) -> Any:
    if path == "-":
        text = sys.stdin.read()
    else:
        text = Path(path).read_text(encoding="utf-8")
    return json.loads(
        text,
        object_pairs_hook=reject_duplicate_keys,
        parse_float=reject_float,
        parse_constant=reject_float,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "certificate",
        nargs="?",
        default=str(DEFAULT_CERTIFICATE),
        help="Full v2 certificate path, or '-' to read JSON from standard input.",
    )
    args = parser.parse_args()
    try:
        certificate = load_certificate(args.certificate)
        result = validate_certificate(certificate)
    except (CertificateError, OSError, json.JSONDecodeError) as error:
        print(f"certificate check failed: {error}", file=sys.stderr)
        return 1
    result["certificate"] = args.certificate
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
