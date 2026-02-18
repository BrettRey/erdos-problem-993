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

### 1B) Allocation-based attacks (new)

#### Equal-share allocation (fails)

Tested the proposed split

`P(h) - 1/3 <= sum_{u~h} (1/3 - P(u)) / deg_H(u)`, where `deg_H(u)=|N(u) cap H|`.

Repro script:

```bash
python3 conjecture_a_equal_share_allocation.py --min-n 23 --max-n 23 \
  --out results/whnc_equal_share_n23.json
```

Result through `n<=23` (`931,596` d_leaf<=1 trees; full run done in-session):
- Fails heavily (`277,921` heavy-vertex failures).
- First failures at `n=8`.
- Worst ratio `(P(h)-1/3) / RHS = 4.0` at `n=23`.

So equal-share is too rigid.

#### Edge-constrained transport feasibility (holds)

Defined bipartite transport from heavy excess to neighbor deficits:
- demands on `H`: `d(h)=P(h)-1/3`
- supplies on `U=N(H)`: `s(u)=1/3-P(u)`
- flow allowed only on edges `(h,u)` with `u~h`.

Max-flow feasibility holds for all tested trees:

```bash
python3 conjecture_a_flow_allocation.py --min-n 3 --max-n 23 --workers 8
```

Result:
- Checked all `931,596` d_leaf<=1 trees through `n=23`.
- Feasibility failures: `0`.
- Maximum unmet-demand gap: `0`.

Per-`n=23` certificate saved as:
- `results/whnc_flow_alloc_n23.json`

This is strictly stronger evidence than WHNC itself: every instance admits a
local edge-respecting compensation plan.

#### Weighted Hall subset scan (strong structure)

For non-empty `S subseteq H`, define Hall slack

`slack(S) = supply(N(S)) - demand(S)`.

Scanned all non-empty subsets of `H` (bitmask DP) and found:

```bash
python3 conjecture_a_hall_subset_scan.py --min-n 23 --max-n 23 \
  --out results/whnc_hall_subset_n23.json
```

And separately for `n<=20,21,22`:
- `n<=20`: all `77,141` trees
- `n=21`: all `98,581` trees
- `n=22`: all `227,678` trees
- `n=23`: all `528,196` trees

In every tested tree (through `n=23`):
- Minimum non-empty Hall slack is attained by a singleton subset `S={h}`.
- Argmin subset size distribution is **100% at size 1**.
- Hall failures: `0`.

This suggests a sharpened conjecture:

`min_{nonempty S subseteq H} slack(S) = min_{h in H} slack({h}).`

If proved, WHNC and full transport feasibility follow immediately from local
singleton bounds.

### Additional diagnostics from this pass

#### Private-neighbor strategy is false in this form

Two candidate structural statements were tested and both fail (even with `d_leaf<=1`):

1. **Strong form (false):** for every nonempty `S subseteq H`, every `h in S` has a
   private neighbor (a `u` with `N(u) cap S = {h}`).
2. **Weak form (false):** for every nonempty `S subseteq H`, at least one `h in S`
   has a private neighbor.

Repro script:

```bash
python3 conjecture_a_private_neighbor_scan.py --min-n 3 --max-n 23 \
  --stop-on-first --out results/whnc_private_neighbor_first.json
```

First counterexamples already at `n=8`:
- `g6 = G?B@dO`
- `H = [0,1,2,3,4]`
- Strong-form failure witness: `S = {0,1}` (`0` private, `1` not private).
- For `S = {0,1,2,3}`, every neighbor of each `h in S` is shared.

So a proof through universal private-neighbor existence cannot work directly.

#### Stronger singleton-dominance law holds (through n=23)

Beyond global singleton argmin, we checked the stronger statement:

`slack(S) >= min_{h in S} slack({h})` for every non-empty `S subseteq H`.

Repro script:

```bash
python3 conjecture_a_singleton_dominance_scan.py --min-n 3 --max-n 23 \
  --out results/whnc_singleton_dominance_n23.json
```

Result through all `931,596` d_leaf<=1 trees (`n<=23`):
- Global singleton-argmin failures: `0`.
- Strong local singleton-dominance failures: `0`.
- Argmin subset size distribution: `{1: 931596}`.

This is stronger than the earlier Hall-scan summary and further supports a proof
strategy based on singleton lower bounds.

#### Degree-only fractional allocation attempt fails

Tested the purely combinatorial weighting

`x_{h,u} = 1/deg(h)` on heavy-side edges (which would prove WHNC if
`sum_{h~u} 1/deg(h) <= 1` for every `u in N(H)`).

Repro script:

```bash
python3 conjecture_a_degree_allocation_scan.py --min-n 3 --max-n 23 \
  --out results/whnc_degree_allocation_n23.json
```

This fails frequently for d_leaf<=1 trees (`n>=8`):
- Fail trees through `n<=23`: `650,763 / 931,596`.
- Worst observed `sum_{h~u} 1/deg(h) = 4.5` at `n=23`.

So a proof cannot ignore the probability weights `P(v)` and rely only on
heavy-neighborhood incidence degrees.

#### Local-overlap failures are almost always leaf-driven (diagnostic)

Repro script:

```bash
python3 conjecture_a_overlap_leaf_split.py --min-n 3 --max-n 23 \
  --out results/whnc_overlap_leaf_split_n23.json
```

This scans local-overlap failures (`local_excess(u) > deficit(u)`) through
`n<=23` and splits by whether `u` has at least one heavy leaf neighbor.

Counts:
- Total local-overlap failures: `23,277`.
- With a heavy leaf adjacent to `u`: `23,245`.
- Without any heavy leaf adjacent to `u`: `32`.

So while local-overlap control is false globally, the obstruction is highly
concentrated in configurations where an overloaded `u` sees heavy leaves.
This suggests a proof decomposition by handling leaf-heavy overlaps first.

#### Submodularity/marginal-removal scan (new)

For

`F(S) = supply(N(S)) - demand(S)`,

define removal marginal

`M(h,S) = F(S) - F(S\\{h}) = supply(N(h)\\N(S\\{h})) - demand(h)`.

A natural singleton-argmin proof attempt is a peeling argument:
for every non-singleton `S`, find `h in S` with `M(h,S) >= 0`, then delete `h`.
If always possible, this reduces to a singleton and proves
`F(S) >= min_{h in S} F({h})`.

Repro script:

```bash
python3 conjecture_a_marginal_scan.py --min-n 3 --max-n 22 \
  --out results/whnc_marginal_scan_n22.json
```

Status (through `n<=22`):
- Peeling criterion is **false**: many subsets have `M(h,S) < 0` for all `h in S`.
- Quantitatively: among `403,400` d_leaf<=1 trees (`n<=22`), `193,208` contain at
  least one all-negative subset, with `665,029` all-negative subsets total.
- But these all-negative subsets are **always leaf-involving**:
  - `allneg_subsets_without_leaf = 0` through `n<=22`.
  - Every all-negative subset contains at least one heavy leaf.
- Even for all-negative subsets, singleton dominance still has positive slack:
  - minimum observed `F(S) - min_{h in S} F({h})` is `0.022372600469548365` (at `n=21`).

Artifacts:
- `results/whnc_marginal_scan_n20.json`
- `results/whnc_marginal_scan_n21.json`
- `results/whnc_marginal_scan_n22.json`

This supports a two-case proof strategy for singleton argmin:
1) non-leaf-heavy subsets (where all-negative marginals may be impossible),
2) leaf-heavy subsets (the only observed obstruction class).

#### Non-leaf heavy case appears structurally solved (computationally)

Define `H_nl = {h in H : deg(h) >= 2}`.
For non-empty `S subseteq H_nl` with `|S| >= 2`, check whether some `h in S`
has a private neighbor relative to `S`.

Repro script:

```bash
python3 conjecture_a_nonleaf_private_scan.py --min-n 3 --max-n 23 \
  --out results/whnc_nonleaf_private_n23.json
```

Result through all `931,596` d_leaf<=1 trees (`n<=23`):
- Failures: `0` (no non-leaf subset without a private-neighbor vertex).

Consequence with edge surplus:
if `h` has a private neighbor `u`, then
`M(h,S) = supply(N(h)\\N(S\\{h})) - demand(h) >= supply(u) - demand(h) > 0`.
So every non-leaf subset has at least one strictly positive removal marginal.

This isolates the singleton-argmin obstruction to subsets containing heavy leaves.

#### Leaf-block certificate prototype (rewrite + exact arithmetic)

A dedicated prototype now emits exact (Fraction) certificates for the hard
leaf-heavy all-negative subsets.

Repro scripts:

```bash
python3 conjecture_a_leaf_block_certificate.py --min-n 8 --max-n 20 \
  --max-certs 40 --out results/whnc_leaf_block_certs_n20.json

python3 conjecture_a_leaf_block_certificate.py --min-n 21 --max-n 21 \
  --max-certs 40 --out results/whnc_leaf_block_certs_n21.json
```

For each certificate subset `S`, it records exact values of:
- `F(S)`,
- `m = min_{h in S} F({h})`,
- `gap = F(S) - m`,
- singleton credits `credit(h) = F({h}) - m`,
- overlap costs `overlap_cost(u) = (deg_S(u)-1) * supply(u)`.

and verifies the canonical identity:

`F(S) - m = (|S|-1)m + sum_h credit(h) - sum_u overlap_cost(u)`.

Coverage:
- `n<=20`: `108,447` leaf-heavy all-negative subsets scanned.
- `n=21`: `159,204` leaf-heavy all-negative subsets scanned.
- For all kept hardest certificates, identity check passes exactly.

Hardest recorded certificates:
- `n=20` hardest kept gap: `15/514 ≈ 0.029182879377`.
- `n=21` hardest kept gap: `162/7241 ≈ 0.022372600470`.

This does not yet prove the leaf-heavy case, but it gives a formal algebraic
normal form to target and concrete minimal witnesses for inequality design.

#### General-tree check: singleton-Hall is not universal

Running `conjecture_a_hall_subset_scan.py --all-trees` on general trees
(`n=10..16` tested) shows:
- many Hall failures (`17,625`/`32,413` trees in that range),
- singleton argmin only ~`9.65%` overall,
- worst slack strongly negative.

Hence the singleton-min phenomenon appears specific to the `d_leaf<=1` regime,
which supports using that structural hypothesis in the proof.

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

## 3) Packaged submodularity scans + analytic updates (2026-02-18)

Today’s in-session ad hoc scans were packaged into reproducible scripts and rerun.

### Packaged scripts

- `conjecture_a_local_overlap_profile.py`
  - profiles local-overlap failure by key `(deg_H(u), heavy_leaf_count(u))`.
- `conjecture_a_leaf_signature_scan.py`
  - scans all-negative leaf-heavy subsets and aggregates structural signatures.
- `conjecture_a_leaf_overload_exceptions.py`
  - isolates overload-key behavior and extracts rare non-`leaf_S=1` exceptions.
- `conjecture_a_signature321_slack.py`
  - verifies identity/positivity for the dominant signature `(3,2,1,(2,2))`.

### New artifacts

- `results/whnc_local_overlap_profile_n23.json`
- `results/whnc_leaf_signature_allneg_n21.json`
- `results/whnc_leaf_overload_exceptions_n21.json`
- `results/whnc_signature321_slack_n21.json`

### Key numerical outcomes

From `whnc_local_overlap_profile_n23.json` (full `n<=23` frontier):

- totals: `seen=23,942,357`, `considered=931,596`, `u_count=7,749,703`,
  `fail_count=23,277`.
- leaf-overlap concentration is explicit:
  - `(deg_H,leaf_H)=(3,1)`: `9,837 / 236,006` fails,
  - `(4,1)`: `10,583 / 38,205`,
  - `(5,1)`: `2,580 / 3,461`,
  - `(6,1)`: `230 / 233`,
  - `(8,1)`: `1 / 1`.
- non-leaf keys remain tiny:
  - `(4,0)`: `2 / 5,136`,
  - `(5,0)`: `11 / 544`,
  - `(6,0)`: `17 / 37`,
  - `(7,0)`: `2 / 2`.

From `whnc_leaf_signature_allneg_n21.json`:

- all-negative leaf-heavy subsets (`n<=21`): `267,651`.
- dominant signature: `(k,l,r,degs)=(3,2,1,(2,2))`, count `104,976`.
- global minimum gap remains `0.022372600469548365` at this signature.

From `whnc_leaf_overload_exceptions_n21.json`:

- subset-with-overload count: `10,815`.
- subset with at least one overloaded `u` and `leaf_S(u)=1`: `10,808`.
- subset with any non-`leaf_S=1` overload: `7` total.
- best gap among those 7 exceptions: `0.28533719979180583` (far from extremal).

From `whnc_signature321_slack_n21.json`:

- signature instances checked: `104,976`.
- identity failures: `0`,
- positivity failures: `0`,
- structure failures: `0`.

### Formal lemmas/proof note

See `notes/whnc_submodularity_progress_2026-02-17.md`:

- Lemma A: non-leaf private-neighbor peeling terminates at singleton and yields
  singleton dominance on `H_nonleaf`.
- Lemma B: dominant signature `(3,2,1,(2,2))` has explicit algebraic slack
  formula and strict positivity from Theorem 6 edge bounds.

## 4) New analogy trial: statistical-physics decimation (2026-02-18)

To try a genuinely different lane, we decimated leaves exactly and rewrote the
problem as an inhomogeneous hard-core model on the leaf-stripped core.

### Exact decimation identities

For `d_leaf<=1` tree `T`, let:

- `L` = leaves,
- `A` = supports adjacent to leaves,
- `C = V(T) \ L` (core graph).

On `C`, define vertex activities:

- `lambda_v = 1/2` for `v in A`,
- `lambda_v = 1` for `v in C\A`.

Then at `lambda=1`:

- `Z_T = 2^{|L|} * Z_C^(lambda_v)`,
- `P_T(v) = P_C(v)` for `v in C`,
- `P_T(leaf l with support s) = (1 - P_C(s))/2`.

Hence:

`mu(T) = |A|/2 + sum_{v in C\A} P_C(v) + (1/2) sum_{v in A} P_C(v).`

And equivalently:

`n/3 - mu(T) = (|C|/3 - |A|/6) - [sum_{v in C\A} P_C(v) + (1/2) sum_{v in A} P_C(v)].`

Verified computationally (`conjecture_a_decimation_core_model.py`):

```bash
python3 conjecture_a_decimation_core_model.py --min-n 3 --max-n 21 \
  --out results/whnc_decimation_core_n21.json
```

Result:

- `seen=3,490,527`, `considered=175,722`,
- `z_fail=0`, `core_prob_fail=0`, `leaf_prob_fail=0`, `gap_identity_fail=0`.

### Decimated weighted-WHNC (new structural inequality)

On `C`, define heavy set:

`H = {v in C\A : P(v) > 1/3}`.

Define weighted supply:

- `s(u)=1/3-P(u)` for `u in C\A`,
- `s(u)= (1/2)*(2/3-P(u))` for `u in A`.

Candidate inequality:

`sum_{h in H}(P(h)-1/3) <= sum_{u in N(H)} s(u).`

Scanned by `conjecture_a_decimated_whnc.py`.

Through full `n<=23` frontier:

```bash
python3 conjecture_a_decimated_whnc.py --min-n 3 --max-n 23 \
  --out results/whnc_decimated_whnc_n23.json
```

Result:

- `seen=23,942,357`, `considered=931,596`,
- **global failures: `0`**,
- minimum global margin: `0.27615847319174214`.

Also tested local weighted overlap at each `u in N(H)`:

`sum_{h~u, h in H}(P(h)-1/3) <= s(u)`.

This stronger local form is almost true but not universal:

- local failures: `54` through `n<=23`,
- worst local margin: `-0.10158506746109458`.

Interpretation: decimation substantially regularizes the overlap pathology
(global weighted inequality is empirically robust), while preserving a small
residual nonlocal overlap effect.

### Decimated marginal scan (new strongest empirical signal)

In the decimated model, define for non-empty `S subseteq H`:

`M(h,S)=s(N_priv(h,S))-(P(h)-1/3)`.

A peeling proof would follow if every `S` has some `h` with `M(h,S) >= 0`.
So we scanned for all-negative subsets (`M(h,S)<0` for all `h in S`).

```bash
python3 conjecture_a_decimated_marginal_scan.py --min-n 3 --max-n 21 \
  --out results/whnc_decimated_marginal_scan_n21.json

python3 conjecture_a_decimated_marginal_scan.py --min-n 22 --max-n 23 \
  --out results/whnc_decimated_marginal_scan_n23_tail.json
```

Merged artifact:

- `results/whnc_decimated_marginal_scan_n23.json`

Result through full `n<=23`:

- `seen=23,942,357`, `considered=931,596`,
- `trees_with_h=674,393`,
- `trees_with_allneg=0`,
- `allneg_subsets=0`.

So the all-negative marginal obstruction observed in the original (non-decimated)
WHNC scan vanishes completely after exact leaf decimation on the full frontier.

### Decimated peeling theorem attempt (new)

See `notes/decimated_peeling_proof_attempt_2026-02-18.md`.

Main point:

- In the decimated model, every non-empty `S subseteq H` has a private-neighbor
  heavy vertex (forest degree argument + `deg_C(h)>=2` for `h in H`).
- Combined with Theorem 6 edge surplus, this gives a strict positive removal
  marginal at each peeling step.
- Therefore decimated Hall feasibility and singleton-argmin follow analytically
  (no brute-force subset optimization needed).

The remaining gap is the final bridge from decimated Hall to `n/3-mu(T)`:
an extra support penalty term `-|A cap N(H)|/6` appears in the exact identity.

### Gap-formula Steiner scan (exact full frontier)

To avoid the support-penalty bridge term, we also scanned the gap-formula
weights directly:

- `s_gap(u)=1/3-P(u)` for `u in C\\A`,
- `s_gap(u)=(1/2)(1/3-P(u))` for `u in A`.

Script:

```bash
python3 conjecture_a_steiner_gap_scan.py --min-n 3 --max-n 23 \
  --out results/whnc_steiner_gap_scan_n23.json
```

Result through full `n<=23` (`931,596` d_leaf<=1 trees):

- `F_gap(S) >= 0` for all `2,881,985` non-empty subsets `S subseteq H_core`
  (`fgap_fail=0`),
- for all `1,580,936` non-singleton subsets, some Steiner leaf has
  `M_gap(h,S) >= 0` (`steiner_fail=0`),
- strict margins on this frontier:
  - `minimum F_gap = 0.1380801129345332`,
  - `minimum best-Steiner marginal = 0.0031308091184815007`.

This is currently the strongest computational support for the Steiner-proof lane.

### Mode-vs-mean bridge: tiepoint scan (new)

To push the remaining gap (`mode` from `mu`) without assuming full LC, we scanned
adjacent tie fugacities:

`lambda_k = i_{k-1}/i_k`, with `mu(lambda_k)` from the tilted IS distribution.

Script:

```bash
python3 conjecture_a_tie_mean_scan.py --min-n 3 --max-n 23 \
  --out results/whnc_tie_mean_scan_n23_modefocus.json
```

Results through full `n<=23` (`931,596` d_leaf<=1 trees):

- all adjacent ties checked: `11,172,104`,
- global lower bound `mu(lambda_k) >= k-1` is **false** (2 failures; worst margin
  `-0.04470192405278617` at `n=22`, `k=12`, `lambda=33`),
- global upper bound `mu(lambda_k) <= k` has 0 failures.

At `lambda=1`, still:

- `mode <= ceil(mu)`: 0 failures,
- Darroch set-membership (`mode in {floor(mu), ceil(mu)}`): 0 failures.

Most promising refinement:

- For `m = mode(I_T)` (at `lambda=1`), focused tie condition
  `mu(lambda_m) >= m-1` had **0 failures** over all `931,596` trees,
  with minimum margin `0.36225903819956606`.

This suggests a narrower proof target: prove the single mode-centered tie
inequality, then lift via monotonicity of `mu(lambda)`.

Dedicated packaged run:

```bash
python3 conjecture_a_mode_tie_focused_scan.py \
  --min-n 3 --max-n 23 \
  --out results/whnc_mode_tie_focused_dleaf_n23.json
```

- `d_leaf<=1` summary (`n<=23`): `checked=931,596`, `fail=0`,
  `minimum_margin=0.36225903819956606`.

All-tree extension:

```bash
python3 conjecture_a_mode_tie_focused_scan.py --all-trees \
  --min-n 3 --max-n 22 \
  --out results/whnc_mode_tie_focused_alltrees_n22.json
```

- all-tree summary (`n<=22`): `checked=9,114,283`, `fail=0`,
  `minimum_margin=0.36225903819956606`.
- global minimum witness in both runs is the balanced length-2 spider at `n=21`
  (`g6=T???????C?G?G?C?@??G??_?@??@???_B~o?`, degree signature `{1:10,2:10,10:1}`).
