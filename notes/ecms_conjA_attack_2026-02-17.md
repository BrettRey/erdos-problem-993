# ECMS + Conjecture A Attack (2026-02-17)

Branch: `attack-ecms-conjA`

This note records two new computational attack directions:
- an ECMS bridge inequality via endpoint deletions,
- a weighted heavy-neighborhood compensation inequality for Conjecture A.

## 1) Conjecture A attack: weighted heavy-neighborhood compensation

Let `P(v)` be the hard-core occupation probability at fugacity `lambda = 1`.
Define
- `H = {v : P(v) > 1/3}`,
- `N(H) = (union_{h in H} N(h)) \ H`.

Candidate inequality (WHNC):

`sum_{h in H} (P(h) - 1/3) <= sum_{u in N(H)} (1/3 - P(u))`.

Interpretation: total heavy excess is compensated by deficit on neighbors of the heavy set.

If WHNC holds, then `mu = sum_v P(v) <= n/3`, because vertices outside `H` all have
`P(v) <= 1/3` and therefore contribute nonnegative deficit.

### Exact verification (Fraction arithmetic)

Computation was run in-session with the same formulas now packaged in
`attack_conjA_weighted_compensation.py`.
Reproduction command:

```bash
python3 attack_conjA_weighted_compensation.py --min-n 3 --max-n 23
```

Observed output summary:
- Checked `931,596` trees with `d_leaf <= 1` (through `n = 23`).
- Failures: `0`
- Equalities: `0`
- Global minimum margin `deficit - excess`: `0.020952965577536575` (strictly positive).

This gives a concrete strengthening candidate for the mean-bound route to Conjecture A.

### 1A) Local overlap control test (single-vertex compensation)

To test the strongest form of overlap control, define for each `u in N(H)`:
- `deficit(u) = 1/3 - P(u)`,
- `local_excess(u) = sum_{h in N(u) cap H} (P(h) - 1/3)`.

The local domination condition is:

`local_excess(u) <= deficit(u)` for all `u in N(H)`.

This was checked by `conjecture_a_overlap_control.py`:

```bash
python3 conjecture_a_overlap_control.py --min-n 3 --max-n 23
```

Result:
- **False** in general (first failures at `n = 14`).
- Through `n <= 23`: `23,277` local failures across `7,749,703` tested vertices `u`.
- Worst ratio `local_excess/deficit = 1.839520398220` at `n=23`.
- Worst witness has `|N(u) cap H| = 8` heavy neighbors.

So strict per-vertex compensation is too strong; WHNC is genuinely global.

### 1C) Extremal margin scan beyond n=23

A parallel scanner (`conjecture_a_whnc_extremal.py`) was run for larger `n`:

```bash
python3 conjecture_a_whnc_extremal.py --min-n 21 --max-n 25 --workers 8
```

Key minima (all with zero WHNC failures):
- `n=21`: min margin `0.028409879091549715`, minimizer spider with arms `[2 x 10]` (`S(2^10)`).
- `n=22`: min margin `0.09726535541752991`, minimizer is a path.
- `n=23`: min margin `0.020952965577536575`, minimizer `S(2^11)`.
- `n=24`: min margin `0.097265355773`, minimizer is a path.
- `n=25`: min margin `0.015296795552874954`, minimizer `S(2^12)`.

This falsifies the working guess that the WHNC minimizer is `S(2^k,1)`. For odd `n`,
the extremal appears to be the balanced spider `S(2^k)`; for even `n`, a path.

### Closed-form WHNC margin for balanced spider S(2^k)

Let `A = (2/3)^(k-1)` and `n = 2k + 1`.
In `S(2^k)`, heavy vertices are exactly the leaves, and each leaf is paired with one
middle vertex (no overlap in `N(H)`).

Using cavity formulas:
- `P(mid) = 1 / (3 + 2A)`,
- `P(leaf) = (1 + A) / (3 + 2A)`.

Therefore:

`margin(S(2^k)) = k * (2/3 - P(mid) - P(leaf)) = k*A / (3*(3 + 2A)).`

Asymptotically:

`margin(S(2^k)) ~ (k/9) * (2/3)^(k-1) -> 0.`

So WHNC is strictly positive in tested range but not bounded away from zero along this
extremal family.

## 2) ECMS attack: contraction mode vs endpoint-deletion modes

For edge `e = uv` in `T`, define:
- `m_T = mode(I(T))`,
- `m_C = mode(I(T/e))`,
- `m_u = mode(I(T-u))`,
- `m_v = mode(I(T-v))`,
- `m_* = max(m_u, m_v)`.

Candidate bridge inequality:

`|m_C - m_*| <= 1`.

### Exhaustive verification through n = 18

Computation was run in-session with the same formulas now packaged in
`attack_ecms_deletion_bridge.py`.
Reproduction command:

```bash
python3 attack_ecms_deletion_bridge.py --min-n 3 --max-n 18
```

Observed output summary:
- Checked `3,348,674` edges (all non-isomorphic trees through `n = 18`).
- Bridge failures (`|m_C - m_*| > 1`): `0`
- One-sided ECMS failures (`m_C < m_T - 1`): `0`
- Full ECMS failures (`|m_C - m_T| > 1`): `0`

Distributions:
- `m_C - m_*`: `-1` in `1,252,868` edges, `0` in `2,095,224`, `+1` in `582`.
- `m_C - m_T`: `-1` in `1,341,769` edges, `0` in `2,004,755`, `+1` in `2,150`.

The bridge statistic is much tighter than raw ECMS shift and may offer a better handle
for a proof attempt (localizing contraction mode against endpoint deletions).
