# Lean-Friendly Bernstein Replacement for the Universal PBD Finite Cells

Date: 2026-07-10
Status: exact certificate complete; generated Lean stubs compile under the
repository's Lean 4.28.0 toolchain

## Result

The thirteen finite scalar cells in the universal Poisson-binomial
effective-drop proof do not require Sturm theory. After translating each cell
by

\[
H=m+t,\qquad 0\le t\le1,
\]

the numerator has a full-degree Bernstein expansion

\[
P_m(m+t)=\sum_{i=0}^{d_m}c_{m,i}
  \binom{d_m}{i}t^i(1-t)^{d_m-i}
\]

in which every coefficient \(c_{m,i}\) is strictly positive. This covers:

- the exact asymmetric numerator on \(3\le H\le4\), with degree 10 and
  11 positive coefficients;
- all twelve symmetric numerators on \(m\le H\le m+1\),
  \(m=4,\ldots,15\), with degrees \(2m+2\) and 264 positive coefficients.

There are 275 coefficients in total. Since every Bernstein basis function is
nonnegative on \([0,1]\), and at least one endpoint basis function is positive
at every point, these expansions prove strict positivity throughout all
thirteen closed cells. The denominators are already known to be positive:

\[
4H^5(H+1)^3
\]

for the asymmetric cell and \(4H^{2m}\) for compact cell \(m\).

This is logically interchangeable with the finite Sturm portion of the
universal proof, but much easier to replay in Lean: each cell reduces to one
`ring` identity and positivity of rational coefficients and Bernstein basis
monomials.

## Exact Check

Run:

```bash
python3 scripts/verify_universal_pb_finite_bernstein.py
```

The script independently reconstructs the asymmetric and symmetric window
forms, extracts the rational numerator on every cell, converts power
coefficients to the Bernstein basis, and then checks exactly that:

1. all 275 coefficients are strictly positive;
2. every Bernstein expansion reconstructs its numerator identically;
3. the expected thirteen cells and degrees are present.

The compact machine-readable result is:

```text
results/universal_pb_finite_bernstein_certificate_2026-07-10.json
```

It records per-cell counts, minimum coefficients, size statistics, and SHA-256
digests. The combined cell digest is

```text
c920cc3bc11eb1564047645ef6b8dd4221efa834b439115bcbad6fb8fdfe4330
```

## Lean Emission

The same script can emit exact Lean-friendly coefficient lists, positivity
proofs, and symbolic identity stubs:

```bash
python3 scripts/verify_universal_pb_finite_bernstein.py \
  --emit-lean /tmp/universal_pb_finite_bernstein_stubs.lean
lake env lean /tmp/universal_pb_finite_bernstein_stubs.lean
```

The generated file has 13 `List ℚ` coefficient tables, 13 `norm_num`
positivity proofs, and 13 real-polynomial identities discharged by `ring`.
The command above compiled successfully under Lean 4.28.0 on 2026-07-10.

## Formalization Consequence

For mechanical formalization, the universal theorem is now preferable to the
redundant direct three-reflected-factor proof. The universal proof artifacts
are complete, and its only formerly Sturm-specific finite component has a
direct positivity certificate. The best first formalization packet is still
the endpoint-aware curvature propagation and mass-window layer; the generated
Bernstein identities can then be imported as a separate finite certificate.
