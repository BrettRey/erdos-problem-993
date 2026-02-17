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
