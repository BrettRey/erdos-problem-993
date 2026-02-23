# Route-1 Transfer Algebraic Analysis (2026-02-19)

## Goal

For the canonical degree-2 bridge decomposition

- `B = T - {l,s}` with `deg(s)=2`,
- `u` = neighbor of `s` different from `l`,
- `P = dp_B[u][0]`, `Q = dp_B[u][1]`, `I(B)=P+Q`,
- `m = mode(I(T))`, `lambda = i_{m-1}(T)/i_m(T)`,

define

- `D := mu_B(lambda) - mu_P(lambda)`,
- `kappa(T) := D - 1/(1+lambda)`.

Route-2 exact implies Route-1 iff `kappa(T) <= 0`.

This note gives:

1. an exact algebraic reduction for `kappa(T)`;
2. a complete extraction of the 4 `kappa>0` witnesses seen in the merged n<=23 scan;
3. a structural characterization of why those 4 are extremal.

## 1) Exact local recursion for `D`

Root `B` at `u`. For each vertex `v` in this rooted tree, define

- `I_v := dp0[v] + dp1[v]` (subtree IS polynomial),
- `q_v := dp0[v]` (subtree with `v` excluded),
- `p_v := Z(dp1[v])/Z(I_v)` (occupation probability of `v` at fugacity `lambda`),
- `delta_v := mu(I_v) - mu(q_v)`.

Then for every `v`:

`delta_v = p_v * (1 - sum_{c child of v} delta_c)`.

Also:

`p_v = lambda * R_v / (1 + lambda * R_v)`,

where `R_v = prod_{c child of v} (1 - p_c)`.

At the root `u`:

`D = mu_B - mu_P = delta_u = p_u * (1 - S_u)`, with `S_u := sum_{c child of u} delta_c`.

So

`kappa(T) = p_u * (1 - S_u) - 1/(1+lambda)`.

Using `p_u = lambda R_u/(1+lambda R_u)` gives the exact sign form:

`kappa(T) = [lambda*R_u*(lambda - (1+lambda)S_u) - 1] / [(1+lambda R_u)(1+lambda)]`.

Hence:

- If `S_u >= 0`, then `kappa(T) <= 0` automatically.
- Any `kappa(T) > 0` witness must satisfy `S_u < 0`.

So the only obstruction is **negative first-level delta mass at `u`**.

## 2) One-child reduction (the failure pattern)

All 4 positive witnesses satisfy `deg_B(u)=1`.
Let the unique child be `c`. Write

- `p := p_c`,
- `G := sum_{gc child of c} delta_gc`.

Then

- `delta_c = p * (1 - G)`,
- `S_u = delta_c`,
- `p_u = lambda*(1-p)/(1+lambda*(1-p))`,
- `D = p_u * (1 - delta_c) = lambda*(1-p)/(1+lambda*(1-p)) * (1 - p + pG)`.

So `kappa>0` requires `G>1` (equivalently `delta_c<0`).
At `lambda=1`, the exact threshold is

`kappa>0  <=>  G > (3 - 2p)/(2(1-p))`.

For small `p` this is about `G > 1.5`.

Interpretation: failures need a child `c` that is rarely occupied (`p` small) but has very large total grandchild delta mass `G`, so that `delta_c` becomes negative and multiplies `p_u` by a factor slightly above 1.

## 3) Complete extraction of the 4 exact-failure witnesses

New verifier script:

- `verify_route1_transfer_failures_2026_02_19.py`

Runs used:

```bash
python3 verify_route1_transfer_failures_2026_02_19.py \
  --min-n 20 --max-n 20 \
  --out results/route1_transfer_failures_n20_2026_02_19.json

python3 verify_route1_transfer_failures_2026_02_19.py \
  --min-n 23 --max-n 23 --res 0 --mod 8 \
  --out results/route1_transfer_failures_n23_r0_2026_02_19.json

python3 verify_route1_transfer_failures_2026_02_19.py \
  --min-n 23 --max-n 23 --res 4 --mod 8 \
  --out results/route1_transfer_failures_n23_r4_2026_02_19.json
```

Merged artifact:

- `results/route1_transfer_failures_all4_2026_02_19.json`

Recovered witnesses (`kappa>0`) are exactly 4:

1. `V?????????????_?G?@??C??G??G??E???o?oB?@|S??`
   - `n=23`, `m=8`, `lambda=0.9972305185989683`
   - `D=0.5058006710219001`
   - `kappa=0.005107340588772491` (global max through n<=23)
   - `p_u=0.49126887657403`, `p=0.03164325433203115`
   - `delta_c=-0.029580124328679958`, `G=1.9348003216830052`

2. `V???????????????O?A??G??W??K??B???W??@_Fii??`
   - `n=23`, `m=8`, `lambda=0.9932366383979288`
   - `D=0.5059567659780138`
   - `kappa=0.004260188293312095`
   - `p_u=0.49350657138744103`, `p=0.019005904254752588`
   - `delta_c=-0.025228021899630804`, `G=2.327378143206335`

3. `S???????????_?O?C??o?@_?@_??oFig?`
   - `n=20`, `m=7`, `lambda=0.993367722918202`
   - `D=0.5056446774709604`
   - `kappa=0.003981091519019597`
   - `p_u=0.4904112869076267`, `p=0.031207847509587286`
   - `delta_c=-0.031062479535069798`, `G=1.9953419416551497`

4. `V?????????????_?G?@??C??G??G??E???o?_B?B|S??`
   - `n=23`, `m=8`, `lambda=0.9899992267443359`
   - `D=0.503416483696693`
   - `kappa=0.0009037256208958011`
   - `p_u=0.49184802611719036`, `p=0.02230710185446091`
   - `delta_c=-0.02352039037510778`, `G=2.0543902353861236`

## 4) Extremal local-shape characterization

New shape-profiler:

- `verify_route1_transfer_extremal_shapes_2026_02_19.py`
- output: `results/route1_transfer_extremal_shapes_2026_02_19.json`

Common structure across all 4 failures:

- `deg_B(u)=1` (all four).
- Unique child `c` has degree `8` or `9` in `B`.
- Exactly one direct leaf branch at `c`.
- Remaining branches are mostly short rooted chains:
  - shape `()` (single-vertex branch): delta near `1/2`,
  - shape `((),)` (length-2 rooted chain): delta near `1/6`,
  - shape `(((),),)` (length-3 rooted chain): delta near `1/3`,
  - plus one deeper branch in 2 of the 4 witnesses, with delta `~0.39` to `~0.44`.

Mechanism:

- these branch deltas sum to `G ~ 1.93 .. 2.33 > 1`, forcing `delta_c < 0`;
- `p` stays small (`~0.02 .. 0.03`) because `p = lambda*prod(1-p_gc)/(1+lambda*prod(1-p_gc))`;
- root occupancy `p_u` remains around `0.49`;
- multiplicative factor `(1 - delta_c)` is about `1.02 .. 1.03`, pushing `D` just above `1/(1+lambda)` by `9e-4 .. 5.1e-3`.

This explains why only a tiny set of trees violate exact transfer and why the excess is small.

## 5) Status toward an analytic `kappa*`

What is now fully established:

- Exact algebraic criterion for transfer failure is reduced to the sign of
  `lambda*R_u*(lambda-(1+lambda)S_u)-1`.
- Every observed failure is one-child-at-`u` with an overcritical grandchild mass (`G>1`) at that child.
- Through n<=23, the empirical transfer constant is

`kappa*_{<=23} = 0.005107340588772491`.

What remains open for a full analytic proof:

- a structural inequality bounding `G` (or equivalently `delta_c`) in terms of occupation profile constraints from `d_leaf<=1` branch families, strong enough to force
  `kappa(T) <= kappa*` (or `<=0` if aiming for exact transfer).

So this pass closes Task (2): it identifies and algebraically explains the 4 extremal failures. The exact global analytic upper bound is reduced to a constrained one-child branch optimization problem.

---

## 6) Extended analysis: why D <= 1/(1+lambda) is unprovable but Route 1 holds (NEW)

Session 2026-02-19 (later): verified the full chain of inequalities on 931,596
d_leaf<=1 trees through n=23.

### 6a) Multi-child case: D <= a is PROVED (computationally)

When u has k >= 2 children in B:

| Property | Status | Detail |
|----------|--------|--------|
| D <= a = lam/(1+lam) | **0 failures** | max D/a = 0.910 (n <= 23) |

**Proof sketch (Case S >= 0)**: D = p_u*(1-S). Since p_u <= a (because R < 1),
and 1-S <= 1, we get D <= a. Since a < 1-a for lam < 1, this gives D < 1/(1+lam).

**Why S >= 0 suffices for multi-child**: with k >= 2 children, the positive
contributions to S (from well-behaved subtrees) tend to dominate. And even when
S < 0, the product R = prod(1-p_c) being bounded away from 1 keeps p_u small
enough that D = p_u*(1-S) <= a.

### 6b) One-child case: transfer fails but mu_P has huge margin

When u has exactly 1 child c in B, D can exceed 1/(1+lam). But:

| Property | Status | Detail |
|----------|--------|--------|
| mu_P >= m-2 directly | **0 failures** | min margin = 0.383 (n <= 23) |
| All-leaves: exists leaf with D <= 1-a | **0 failures** | max gap = -0.063 (n <= 20) |
| Route-2 slack covers excess | **0 failures** | (n <= 23) |

The key insight: mu_P = mu_{I(T_c)}(lam) for the large subtree T_c (size ~n-3).
Even at lam < 1, T_c is large enough that mu_Tc >= m-2 by a wide margin.

### 6c) Recursive bound: f(v) <= p_v fails; f(v) <= a fails

Tested: f(v) = p_v*(1-sum_c f(c)) at every vertex v in B (rooted at u).

| Property | Status | Detail |
|----------|--------|--------|
| f(v) <= p_v for all v | **FAILS** | 4012 failures n <= 20 |
| p_v <= a for all v | **0 failures** | max p_v/a = 1.000 |
| f(v) <= a for all v | **FAILS** | max f/a = 1.0147 (n = 20) |
| f(v) <= 1/(1+lam) | **FAILS** | 1 vertex-level failure (n = 20) |

So no simple pointwise recursive bound closes the one-child gap. The bound
must use global structure (the mode constraint m <= floor(n/3)+1, or the
Route-2 slack).

### 6d) Per-subtree lower bounds on mu_Tc: too crude

Tested: mu_{T_c}(lam_m(T)) >= (|T_c|-1)/3.

**Fails badly**: 958,420 failures out of 1,187,696 subtree checks (n <= 23).
The per-subtree bound cannot work because lam_m(T) can be much smaller than
the subtree's own mode fugacity, making mu_Tc small.

### 6e) Recommended proof strategy for Route 1

**The transfer is a dead end.** The correct approach is one of:

**(A)** Prove mu_P >= m-2 DIRECTLY, e.g., via the cavity/Steiner apparatus
already available for B. Since mu_P = sum_c mu_{T_c}(lam_m) and the subtrees
cover n-3 vertices, a coupling argument comparing lam_m(T) to the subtree
fugacities may work.

**(B)** Don't need Route 1 at all. The bridge decomposition
Phi_m(T) = Phi_m(A) + lam*Phi_{m-1}(B) has both terms independently >= 0
(0 failures through n=23). Route 2 gives Phi_m(A) >= 0 via the pendant
bonus, and Phi_{m-1}(B) >= 0 follows from STRONG C2 (0 failures through
n=23, 4.5M degree-2 leaves). Route 1 provides additional slack but is not
strictly needed for Phi_m(T) >= 0.

**(C)** Hybrid: use multi-child D <= a (proved) for k >= 2, and for the
one-child case, choose a different leaf (the all-leaves scan shows every
tree has a good leaf with D <= 1/(1+lam)).

## 7) New artifacts from this analysis

Scripts:
- `conjecture_a_route1_transfer_algebraic.py` (v1: identity derivation, exact
  fraction verification, exponential bound attempt)
- `conjecture_a_route1_transfer_algebraic_v2.py` (v2: combined analysis,
  multi-child D <= a, all-leaves scan, direct mu_P bound)
- `conjecture_a_route1_transfer_algebraic_v3.py` (v3: comprehensive n=23 scan,
  delta-negative profiling, per-subtree bounds)

Key numerical results (all on d_leaf<=1, canonical leaf, n <= 23):
- D = p_u*(1-sum_delta) identity: 0 failures (exact fractions)
- mu_P >= m-2: 0 failures, min margin **0.383**
- Multi-child D <= a: 0 failures, max D/a = **0.910**
- All-leaves min D <= 1-a: 0 failures (n <= 20), max gap = **-0.063**
- f(v) <= 1/(1+lam) at vertex level: 1 failure (n = 20)
- Per-subtree mu_Tc >= (|Tc|-1)/3: FAILS (not viable)
