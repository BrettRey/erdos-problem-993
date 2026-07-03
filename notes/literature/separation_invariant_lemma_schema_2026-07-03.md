# Separation Invariant Lemma Schema
Date: 2026-07-03

## Purpose

This is a proof-facing draft, not a claimed theorem. Its job is to turn the empirical pattern from the Ramos--Sun, Galvin, Bautista--Ramos, Li, `n=28`, forest-product, product-power, and beam-search audits into a small set of precise obligations.

The key shift is:

> Do not try to prove log-concavity, or even "almost log-concavity." Prove that any log-concavity-breaking ratio bump after the first descent has its right endpoint still below `1`.

That is exactly the kind of thing the current stress data supports. LC defects can be large as LC ratios, and they can be multiple, but their absolute right-hand ratios are deep-tail quantities, not valley-making quantities.

## Definitions

Let

```text
P(x) = sum_{k=0}^alpha a_k x^k
```

be an independence polynomial with strictly positive coefficients on `0 <= k <= alpha`. Define consecutive coefficient ratios

```text
r_k(P) = a_k / a_{k-1},   1 <= k <= alpha.
```

Let the first strict descent be

```text
D(P) = min { k : r_k(P) < 1 },
```

with `D(P) = +infinity` if no such `k` exists. For tree independence polynomials this only matters after the leftmost mode; plateaux are harmless.

An LC defect at index `k` is the inequality

```text
a_{k-1} a_{k+1} > a_k^2.
```

In ratio language this is exactly

```text
r_{k+1}(P) > r_k(P).
```

Define the late-bump pressure

```text
B(P) = max { r_{k+1}(P) :
             k >= D(P),
             r_{k+1}(P) > r_k(P) },
```

with `B(P) = 0` if the set is empty. Also define the post-descent tail pressure

```text
R(P) = max { r_j(P) : j >= D(P) }.
```

The computations show that `R(P)` can be close to `1`, but `B(P)` remains far below `1` in the current hard sources.

## Lemma 1: Ratio-Crossing Certificate

**Lemma.** Let `P` have positive coefficients. If `B(P) < 1`, then the coefficient sequence of `P` is unimodal.

Equivalently, if `P` is not unimodal, then there is an LC defect after the first descent whose right-hand ratio exceeds `1`.

**Proof sketch.** If `P` is not unimodal, then after the first strict descent there is a later strict rise. In ratio terms, there are indices `s < t` with `r_s < 1` and `r_t > 1`, where all relevant indices are after the first descent. Moving from a ratio below `1` to a ratio above `1` requires some adjacent increase `r_{k+1} > r_k` whose right endpoint is already above `1`. That index `k` is an LC defect and contributes `r_{k+1} > 1` to `B(P)`. Contraposition gives the result.

This lemma is elementary but useful: it says a non-unimodal tree cannot hide in arbitrary LC failure. It must produce a very specific kind of LC failure, namely a post-descent bump whose right endpoint crosses above `1`.

## Lemma 2: Tail-Lock Reduction

Let `L(P)` be any tail-lock index such that

```text
r_j(P) <= 1 for all j >= L(P).
```

For tree independence polynomials, Levit--Mandrescu-type tail monotonicity gives such an `L(P)` in the last part of the sequence.

**Lemma.** To prove `B(P) < 1`, it is enough to control LC bumps with

```text
D(P) <= k < L(P).
```

**Proof sketch.** If `k >= L(P)`, then `r_{k+1}(P) <= 1` by the tail-lock hypothesis. So only the bridge between first descent and tail lock can contain a dangerous bump.

This is the first place where the proof should use known tree-specific tail monotonicity. The separation problem is not the whole tail; it is the bridge between the mode region and the certified decreasing tail.

## Target Lemma 3: Tree Late-Bump Reserve

**Target.** For every tree `T`, with `P = I(T)`, every bridge LC bump satisfies

```text
D(P) <= k < L(P) and r_{k+1}(P) > r_k(P)
    implies
r_{k+1}(P) < 1.
```

Equivalently,

```text
B(I(T)) < 1.
```

Together with Lemma 1, this proves unimodality directly.

This target is logically close to unimodality, so it should not be sold as a major simplification by itself. Its value is narrower: it identifies the exact local obstruction that any proof must rule out. The empirical audits show that known LC failures are not close to this obstruction.

## Stronger, More Useful Version

The proof should aim for a quantitative reserve:

```text
r_{k+1}(I(T)) <= 1 - Delta(T,k)
```

where `Delta(T,k) > 0` is obtained from structural data rather than from coefficient inspection.

Candidate sources for `Delta`:

1. **Mean reserve.** In the `d_leaf <= 1` lane, the manuscript's mean bound gives `mu(T) < n/3`. If a late bump creates `r_{k+1} >= 1` with `k >= floor(n/3)+1`, then a tie-fugacity argument should force the mean at some `lambda <= 1` to be at least `k`, contradicting monotonicity of `mu(lambda)` and the mean bound.

2. **Tie-fugacity margin.** The existing manuscript already uses the condition `mu(lambda_m) >= m - 1` at the mode tie. The separation variant would need the analogous statement for any post-descent upward crossing:

   ```text
   if r_k < 1 <= r_{k+1}, and lambda = 1 / r_{k+1},
   then mu(lambda) >= k.
   ```

   This statement is false for arbitrary coefficient sequences unless additional shape hypotheses are added. The proof task is to find the tree-specific hypothesis that makes it true, or to replace it with the weakest true substitute.

3. **PNP transfer surplus.** Outside the `d_leaf <= 1` lane, multi-leaf hubs are already routed through Hub Exclusion and Transfer. The separation lemma should not require mean reserve for those trees. It should either reduce to residual `d_leaf <= 1` components or use the Case B hub surplus to keep any post-descent crossing below `1`.

4. **Subdivision identity.** For subdivided trees, use

   ```text
   I(T_e) = I(T) + x I(T/e)
   ```

   plus ECMS/combined-tail machinery. The separation version should say that in the gap case, the two ambiguous positions cannot create a post-descent LC bump with right endpoint above `1`.

## Lane Split

### Lane A: `d_leaf <= 1`

Desired theorem:

> If `T` is a tree with `d_leaf(v) <= 1` for every vertex, then every post-descent LC bump of `I(T)` has right-hand ratio below `1`.

Likely proof route:

1. Use the manuscript's `mu(T) < n/3` bound.
2. Use mode-mean or the tie-fugacity condition to put the first descent at or before the low-mode threshold.
3. Use a direct ratio-reserve lemma to rule out `r_{k+1} >= 1` after first descent. A tie-fugacity crossing argument may still be useful for actual upward crossings, but it should not be applied to every harmless LC bump.
4. Use Levit--Mandrescu tail monotonicity to avoid proving anything separately in the final tail.

What must be falsified next:

```text
For every d_leaf <= 1 stress row and every post-descent LC bump,
does the candidate tie-fugacity lower bound hold?
```

The existing product-power and beam-search audits already show the ratio conclusion, but not the tie-fugacity mechanism. The first tie-fugacity audit below shows that the broad LC-bump version of the mechanism is false.

### Lane B: Multi-Leaf Hubs

Desired theorem:

> If `T` has a vertex with at least two leaf-neighbors, then any dangerous post-descent LC bump either disappears under the Hub Exclusion transfer or is bounded by the same separation statement on the residual components.

Likely proof route:

1. Choose a multi-leaf hub `v`, remove `v` and its leaf-neighbors, and pass to the residual forest `T'`.
2. Use the proved Transfer Lemma for 1-Private maximal IS as the combinatorial side condition.
3. Translate the transfer into a coefficient-ratio inequality, not merely a size inequality.
4. Apply the Lane A separation statement on residual components with `d_leaf <= 1`, recursively if needed.

This is the least developed part. The current PNP framework transfers maximal-IS size information cleanly, but it does not yet transfer coefficient-ratio information. If this translation fails, Lane B should remain PNP-based while Lane A carries the separation work.

### Lane C: Subdivision / Degree-2 Structure

Desired theorem:

> Under ECMS and the combined-tail condition, subdivision preserves not only unimodality but late-bump separation.

Likely proof route:

1. Write `I(T_e) = I(T) + x I(T/e)`.
2. If the summand modes differ by at most `1`, use the adjacent-mode sum lemma.
3. In the gap case, prove a two-index ratio inequality at `M` and `M+1`.
4. Check that later LC bumps in the sum have right endpoint below `1` because both summands are already descending.

This lane is closest to existing manuscript machinery and may be the fastest way to get a nontrivial formal lemma.

## Product/Forest Stability Lemma

The product audits suggest another useful lemma.

Let `P` and `Q` be separated polynomials, and let `PQ` be their product. Since coefficients of `PQ` are the convolution of coefficients of `P` and `Q`, unimodality itself is expected to be stable under convolution for positive unimodal sequences. But the separation invariant is stronger and needs its own statement.

Candidate:

> If `B(P) < 1` and `B(Q) < 1`, and the post-descent maxima of `P` and `Q` occur within one step of their modes, then `B(PQ) < 1`.

This is not needed to prove the tree conjecture from trees alone, but it is valuable because the conjecture is also stated for forests and because products amplify near-miss behavior. The product-power audit and beam search are falsification tests for this lemma.

Current evidence:

- Product powers: 184 rows, 52 non-LC, 0 non-unimodal, max post-descent ratio `0.9988814019`, max LC-bump right ratio `0.0166666667`.
- Beam search: 2,866 rows, 2,628 non-LC, 0 non-unimodal, max post-descent ratio `0.9971493049`, max LC-bump right ratio `0.0166666667`.

## Minimal Proof Skeleton

If the target lemmas above can be proved, the proof would read:

1. Let `T` be a minimal counterexample to unimodality.
2. Let `P = I(T)` and let `D = D(P)` be the first descent.
3. Since `T` is a counterexample, Lemma 1 gives an LC defect `k >= D` with `r_{k+1} > 1`.
4. If `k` lies in the final tail, Levit--Mandrescu tail monotonicity contradicts `r_{k+1} > 1`.
5. If `T` has `d_leaf <= 1`, Lane A gives `r_{k+1} < 1`, contradiction.
6. If `T` has a multi-leaf hub, Lane B transfers to smaller residual components until Lane A applies, or gives a direct hub-surplus contradiction.
7. If `T` is a subdivision, Lane C and minimality rule it out under the same ECMS/combined-tail assumptions already used in the manuscript.
8. No case remains.

Step 7 is conditional on the existing ECMS/combined-tail route. If the goal is an unconditional proof, the separation project should first focus on Steps 4--6 and avoid depending on subdivision.

## Immediate Next Computation

The next computation should not be another broad LC corpus audit. It should test the mechanism behind Lane A:

For every row in the current stress artifacts and every post-descent LC bump with `r_{k+1} > r_k`, compute

```text
lambda = 1 / r_{k+1}
mu(lambda)
mu(lambda) - k
```

and separately compute this for any post-descent upward crossing, not only LC defects.

If `mu(lambda) >= k` holds in the tree/family/product stress set, the tie-fugacity route is worth trying to prove. If it fails, the draft should pivot: the ratio separation is still empirically true, but the mean/tie mechanism is not the proof mechanism.

## Tie-Fugacity Audit Outcome

I ran the mechanism test as:

```bash
python3 scripts/audit_tie_fugacity_bumps.py \
  --max-power 20 \
  --mixed-max-power 4 \
  --beam-max-factors 8 \
  --beam-width 80 \
  --unsafe-threshold 0.95 \
  --top 10 \
  --out results/tie_fugacity_bump_audit_2026-07-03.json
```

The audit covered 3,050 product rows: the 184 product-power rows and the 2,866 beam-search rows. It found:

- 8,777 post-descent LC-bump events.
- 0 post-descent upward-transition events.
- 6,590 negative values of `mu(lambda) - k` at LC bumps.
- Minimum `mu(lambda) - k = -3.8776208619`, attained by `TG_{8,8}` at `k = 214`.
- Maximum LC-bump right ratio still only `0.0166666667`.

Interpretation: the broad tie-fugacity mechanism is false for ordinary LC bumps. Tying the two coefficients at a harmless deep-tail LC defect often leaves the mean below the left tied index. This does **not** threaten the separation invariant, because none of these bumps is an upward transition after the first descent. It does mean that a proof should not try to show `mu(lambda) >= k` for all LC bumps.

The surviving proof target is narrower:

```text
if k >= D(P) and r_{k+1}(P) >= 1,
then derive a contradiction from mean reserve, tail lock, or structural transfer.
```

This crossing-only version is exactly what Lemma 1 needs. It remains untested in the current stress set because no such crossings occur; that is itself the observed separation.

## Direct Crossing-Pressure Search Outcome

I then searched for the crossing obstruction directly, using the post-first-descent ratio

```text
max { a_{j+1} / a_j : j >= D(P) }
```

as the objective rather than using LC defects as a proxy.

The evolutionary optimizer was run as:

```bash
python3 scripts/crossing_pressure_optimizer.py \
  --min-n 60 \
  --max-n 180 \
  --step-n 40 \
  --pop-size 70 \
  --generations 120 \
  --archive-size 16 \
  --seed 993 \
  --verbose-every 30 \
  --out results/crossing_pressure_optimizer_2026-07-03.json
```

It found no strict post-descent crossing and no non-unimodal tree. Its best pressures were:

| `n` | Best pressure | Reserve |
|---:|---:|---:|
| 60 | `0.9287861476` | `0.0712138524` |
| 100 | `0.9575083299` | `0.0424916701` |
| 140 | `0.9699738045` | `0.0300261955` |
| 180 | `0.9770457694` | `0.0229542306` |

The optimizer repeatedly converged to broom-like trees with one very high-degree vertex and one small branch vertex. This is useful as a sanity check, but it is not the sharpest stress source.

The deterministic family scan was run as:

```bash
python3 scripts/scan_crossing_families.py \
  --min-n 50 \
  --max-n 500 \
  --step-n 10 \
  --broom-max-path 80 \
  --arm-values 2-6 \
  --max-arms 3 \
  --top 20 \
  --out results/crossing_family_scan_2026-07-03.json
```

It scanned 5,870 broom and multi-arm-star rows, with 0 non-unimodal rows, 0 post-descent upward-transition rows, and 0 non-log-concave rows. The best crossing pressure was

```text
0.991750265410068
```

at the multi-arm star `M(488;2,3,6)` on 500 vertices. This slightly improves the old broom benchmark `broom(33,467)`, which has pressure `0.9916957223198212`.

This changes the next proof target. Brooms are not just examples; they are a ridge in the crossing-pressure landscape. But the best observed ridge has a small finite arm bouquet attached to the hub. The natural next mathematical task is an asymptotic ratio calculation for

```text
M(s;2,3,6)
```

and nearby `M(s; a,b,c)` families, compared against the broom asymptotic. If these families have an explicit reserve of order `1/s` below `1`, they may supply the model inequality for the general hub case.

## Asymptotic Hub-Ridge Outcome

For a hub with `s` pendant leaves and fixed path arms `a_1, ..., a_m`, the independence polynomial has the exact closed form

```text
(1+x)^s prod_j I(P_{a_j}) + x prod_j I(P_{a_j-1}).
```

I added `scripts/asymptotic_hub_ridge.py` to exploit this formula and scan larger hub bouquets without vertex-level DP. The main scan was:

```bash
python3 scripts/asymptotic_hub_ridge.py \
  --s-values 100,200,500,1000,2000,5000 \
  --broom-max-arm 120 \
  --arm-values 2-10 \
  --max-arms 3 \
  --top 40 \
  --self-check \
  --out results/asymptotic_hub_ridge_2026-07-03.json
```

It scanned 1,980 exact hub-bouquet rows and found 0 post-descent upward transitions. The best row in that bounded scan was `broom(118,5000)`, with crossing pressure `0.9992005303423901` and reserve `0.000799469657610`, so `s * reserve = 3.9973482880`.

Because the broom optimum was hitting the arm cap, I ran a focused `s=5000` broom scan:

```bash
python3 scripts/asymptotic_hub_ridge.py \
  --s-values 5000 \
  --broom-arm-values 1-240 \
  --arm-values '' \
  --max-arms 3 \
  --top 20 \
  --self-check \
  --out results/asymptotic_broom_ridge_2026-07-03.json
```

This scanned 240 exact broom rows, again with 0 post-descent upward transitions. The best row was `broom(194,5000)`, with crossing pressure `0.9992098719189995`, reserve `0.000790128081001`, and `s * reserve = 3.9506404050`. The second row was `broom(241,5000)`, so the long-arm broom ridge is not fully optimized yet.

The important point is not the exact maximizing arm. The important point is that the hard hub examples are now visibly governed by a binomial-tail reserve:

```text
max_{j >= D} a_{j+1}/a_j = 1 - about 4/s.
```

For comparison, `M(5000;2,3,6)` has pressure `0.999182328534265` and `s * reserve = 4.0883573287`, close to but below the best broom. Thus `M(s;2,3,6)` was a finite-range ridge, not obviously the asymptotic maximizer.

This gives a cleaner proof subproblem:

> Prove a hub-bouquet reserve theorem: for `I(T) = (1+x)^s Q(x) + x R(x)`, with fixed or controlled arm polynomial `Q`, the first post-descent ratio is at most `1 - c/s + O(s^{-2})`, with `c` bounded away from `0`.

For fixed arms, the `xR(x)` term has bounded degree and is absent near the middle of the leaf-binomial range. The problem reduces to ratios of `(1+x)^s Q(x)`, where `Q` is a fixed path-product polynomial. This looks much more tractable than the original tree recurrence: it is an explicit coefficient-ratio estimate for a binomial convolution.

## Variance Correction

The `1 - c/s` form is the right first target for fixed or controlled arms, but it is not the correct fully general asymptotic statement for arbitrary broom handles. I ran a floating-point probe of the dominant product term `(1+x)^s I(P_l)` at `s = 5000`:

```bash
python3 scripts/probe_broom_variance_scaling.py \
  --out results/broom_variance_scaling_probe_2026-07-03.json
```

This is not a certificate, because it omits the exponentially small hub-included term `x I(P_{l-1})`. It is meant only to identify the right denominator for the proof statement.

The result is decisive enough for strategy:

| Arm `l` | Pressure | `s * reserve` | `n * reserve` | `variance * reserve` |
|---:|---:|---:|---:|---:|
| 193 | `0.9992098719` | `3.9506404050` | `4.1039252527` | `1.0013777758` |
| 5,000 | `0.9992021635` | `3.9891826243` | `7.9791630851` | `1.3541778914` |
| 20,000 | `0.9995314396` | `2.3428020814` | `11.7144789672` | `1.4239332080` |
| 50,000 | `0.9997188657` | `1.4056714875` | `15.4626674969` | `1.6087164720` |

So the robust statement should be variance-based:

```text
max_{j >= D} a_{j+1}/a_j <= 1 - c / Var + lower-order terms,
```

for the coefficient distribution of `(1+x)^s Q(x)`, with the `1 - c/s` statement recovered when `Var = s/4 + O(s)` or `Var = s/4 + O(1)`.

This correction strengthens the proof program by preventing a false target. The hard local ratio is controlled by the width of the binomial-convolved path distribution, not by the number of leaves alone. For fixed arms, these coincide up to constants; for very long broom handles, they do not.

The proof-facing reduction is now separated out in

```text
notes/literature/hub_bouquet_variance_reserve_proof_schema_2026-07-03.md
```

The key new algebraic handle there is the exact conditional-ratio identity for `A(x)=(1+x)^sQ(x)`:

```text
A_{k+1}/A_k = E_k[(s-k+Y)/(k+1-Y)],
```

where `P_k(Y=y) proportional q_y binom(s,k-y)`. This turns the reserve into a conditional leaf-count excess estimate near the first descent.
