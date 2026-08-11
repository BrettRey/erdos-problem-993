# Effective unimodality for every fixed arm vector

## Scope

This note supplies the effective refinement requested in
`notes/next_move_2026-07-16.md`: for each fixed multi-arm vector, it gives a
computable threshold after which adding leaves at the hub preserves
unimodality. It does **not** prove unimodality for arbitrary trees, and it
does not handle arm vectors that grow with the number of hub leaves.

**Prior-art correction.** Every multi-arm star is a spider, and
Li--Li--Yang--Zhang already prove the stronger result that every spider has a
strongly log-concave independence polynomial. Accordingly, the theorem below
is not a new family closure. Its contribution is only an elementary explicit
threshold and finite certificate mechanism for a family already known to be
log-concave for every `s`.

## Theorem

Fix positive integers `a_1,...,a_m` and set

```text
Q(x) = product_j I(P_{a_j};x),
R(x) = product_j I(P_{a_j-1};x),
F_s(x) = (1+x)^s Q(x) + xR(x).
```

Write `e=xR=sum_{k=0}^r e_k x^k`, put `e_{r+1}=0`, and define

```text
D_k = max(0, e_k-e_{k+1})  (0 <= k <= r).
```

Let `S` be the least integer `s >= 2r+1` such that

```text
binom(s,k+1)-binom(s,k) >= D_k  for every 0 <= k <= r.   (1)
```

Then `F_s` is unimodal for every integer `s >= S`. The threshold is
effective: the left side of (1) is nondecreasing for `s >= 2r+1`, and the
finite bound

```text
S <= 2r+1 + (r+1) max_k D_k
```

guarantees termination of a binary search for `S`.

## Proof

Write

```text
d_k(s) = [x^k](1+x)^s Q(x),
f_k(s) = d_k(s)+e_k.
```

Every `Q` coefficient is nonnegative and `q_0=1`. For `0 <= k <= r` and
`s >= 2r+1`, expand the difference as

```text
d_{k+1}-d_k
 = sum_j q_j (binom(s,k+1-j)-binom(s,k-j)),
```

with out-of-range binomial coefficients interpreted as zero. For `j<=k`,
the binomial index `k-j` is at most `r`, so each displayed difference is
nonnegative. The possible `j=k+1` term is also nonnegative. Retaining only
the `j=0` contribution gives

```text
d_{k+1}-d_k >= binom(s,k+1)-binom(s,k).
```

Condition (1) therefore yields

```text
f_{k+1}-f_k
 >= D_k + e_{k+1}-e_k
 >= 0
```

for every `k<=r`. Thus `F_s` is nondecreasing through coefficient `r+1`.

For `k>=r+1`, the correction `xR` vanishes and `F_s` agrees with
`(1+x)^sQ`. Each path independence polynomial is real-rooted with negative
zeros; hence so are `Q` and `(1+x)^sQ`. Newton's inequalities make the
coefficient sequence of `(1+x)^sQ` log-concave and therefore unimodal. At
the join, `f_{r+1}-f_r>=0` also implies `d_{r+1}-d_r>=0`, because
`e_r>=0` and `e_{r+1}=0`. Consequently the unimodal tail of `d` cannot
create a descent-then-rise after the nondecreasing prefix of `F_s`.
Therefore `F_s` is unimodal.

Finally, if `D=max D_k` and
`s=2r+1+(r+1)D`, then for every `k<=r`,

```text
binom(s,k+1)-binom(s,k)
 = binom(s,k)(s-2k-1)/(k+1) >= D >= D_k.
```

This proves the stated finite upper bound and computability. QED.

## Verification boundary

The proof above is symbolic and applies to every fixed positive arm vector,
but the prior spider theorem already implies unimodality without a threshold.
The accompanying script performs three separate exact checks:

1. it computes `S` and replays every inequality defining it;
2. it compares the closed formula with the generic tree DP at small `s`;
3. for selected vectors, it checks every `s<S` for a finite all-`s`
   corollary.

Run:

```text
python3 scripts/fixed_arm_unimodality_threshold.py
python3 -m unittest test_fixed_arm_threshold.py
```

Output:

```text
results/fixed_arm_unimodality_threshold_20260808.json
```

The default exact replay gives:

| arm vector | certified `S` | finite range checked | all `s` certified? |
|---|---:|---:|---|
| `(5,5,4,2)` | 18 | `0..17` | yes |
| `(2,3,6)` | 14 | `0..13` | yes |
| `(6,6,1)` | 16 | `0..15` | yes |
| `(8)` | 12 | `0..11` | yes |

The unit-test kill sweep also covered every nondecreasing arm vector of
arity one through three with arm lengths at most five. It replayed the
threshold inequality and checked exact unimodality through `S+2`; no failure
occurred. This sweep tests the implementation and examples, not the general
proof.

The exact finite scan is evidence only for the named vectors. The eventual
statement for every fixed vector rests on the proof, not on that scan.

## Consequence and limitation

Internally, the threshold proof turns an all-`s` statement into a finite
computation. Strategically, however, this does not advance beyond the known
spider theorem. The relevant next class has multiple interacting hubs, where
the path-product correction is no longer a single spider decomposition. See
`notes/two_hub_ratio_order_killtest_2026-08-08.md` for exact failures of the
most direct ratio-order extension.
