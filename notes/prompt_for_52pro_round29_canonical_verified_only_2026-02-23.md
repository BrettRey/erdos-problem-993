# Round 29: Canonical Verified-Only Gate (no fabricated witnesses)

Context: Erdos #993 closure project, canonical degree-2 bridge setup.

You must use the repository definitions exactly (from `attack4_common.py` and `scripts/verify_route1_moment_class_scan.py`).

## Hard validity gates (must satisfy all)
A claimed witness/counterexample is valid only if all are true:
1. Parsed from a valid g6 string as a tree (`|E| = n-1`).
2. Passes `is_dleaf_le_1(n, adj) == True`.
3. Passes canonical decomposition: `bridge_decomposition(n, adj, g6, require_dleaf=True) is not None`.

Any witness failing any gate is invalid and must be discarded.

## Definitions (exact)
Use these exact quantities from `scripts/verify_route1_moment_class_scan.py`:
- `Threshold = (m-2) - (p_{m-1}*lambda^(m-1) + 2*p_m*lambda^m)/P(lambda)`.
- `B_max` = subset-DP optimum over child-factor orders of summed one-step LP-dual lower bounds in class `C_{deg,mu,mu2,mu3,mu4}`.
- One-step LP constraints:
  `q(i)=c0+c1*i+c2*i(i-1)+c3*i(i-1)(i-2)+c4*i(i-1)(i-2)(i-3) <= beta_i(F,m,lambda)`.

Canonical key (for pair matching) is:
`K4 := (deg(P), m, lambda, mu1, mu2, mu3, mu4)`.

## Required task (binary)
Return exactly one:

A) `SUCCESS`:
- One constructive canonical symbolic lemma that implies `B_max >= Threshold` (or directly `R1_tail2`) for all canonical instances.
- Must provide full chain to `E_route1`.
- No geometric-series tails, no truncation assumptions.

B) `BLOCKED`:
- Explicit canonical pair `(T1,T2)` such that:
  - both pass all three validity gates,
  - same `K4` (and any extra invariants you claim the class uses),
  - different target side (`Threshold` differs or sign of `B_max-Threshold` differs).

If neither is achieved, output `INSUFFICIENT`.

## Mandatory computation evidence
You must include a reproducible Python snippet and its output (or equivalent script output) that:
1. Verifies all three validity gates for every witness you cite.
2. Prints `(m, lambda, Threshold, B_max, B_max-Threshold)`.
3. If claiming a pair, prints both keys to show equality.

## Mandatory output format
1) `Canonical class inputs used`
2) `Lemma or explicit canonical pair`
3) `Derivation / verification` (must include validity-gate outputs)
4) `Binary verdict: SUCCESS or BLOCKED or INSUFFICIENT`
5) `If BLOCKED/INSUFFICIENT: minimal next canonical invariant or class`
