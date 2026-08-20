# Bounded pendant cores with arbitrary connectors — 2026-08-08

> **2026-08-20 update.** The bounded theorem and partial-synchronization
> route below are superseded by a direct clan-cancellation proof draft for
> arbitrary pendant weight and every connector length. See
> `notes/two_hub_B_logconcavity_proof_2026-08-20.md` and
> `notes/c2_connector_clan_reduction_2026-08-20.md`. The original finite
> theorem and its verifier remain valid independent evidence.

## Computer-assisted theorem

Let `T(a,b,t)` be the tree obtained from two hubs as follows. Attach pendant
paths with positive lengths given by the multisets

```text
a=(a_1,...,a_r),   b=(b_1,...,b_s),
```

to the two hubs, and join the hubs by a path having `t >= 0` internal
vertices. If

$$
\sum_i a_i+\sum_j b_j\leq24,
$$

then the independence polynomial of `T(a,b,t)` is log-concave, and hence
unimodal, for every `t >= 0`.

The connector length is unbounded, so this is an infinite-family theorem,
not an order-24 tree scan. The bound 24 applies only to vertices on pendant
paths, excluding the hubs and the internal connector vertices.

This does **not** prove the corresponding statement for pendant weight above
24, the full class of trees with at most two branch vertices, or Erdős
Problem 993.

## Algebraic transfer for a fixed pendant core

Let `P_n=I(P_n;x)`, with `P_{-2}=0` and `P_{-1}=P_0=1`, and put

$$
Q_a=\prod_iP_{a_i},\qquad R_a=\prod_iP_{a_i-1},
$$

with analogous definitions for `b`. Conditioning on the two hubs gives

$$
\begin{aligned}
F_t={}&Q_aQ_bP_t
+x(R_aQ_b+Q_aR_b)P_{t-1}
+x^2R_aR_bP_{t-2}. \tag{1}
\end{aligned}
$$

Define the adjacent-hub and contracted-hub polynomials

$$
\begin{aligned}
B&=Q_aQ_b+x(R_aQ_b+Q_aR_b),\\
C&=Q_aQ_b+xR_aR_b.
\end{aligned} \tag{2}
$$

Then

$$
F_0=B,\qquad F_1=B+xC,
\qquad F_t=F_{t-1}+xF_{t-2}\quad(t\geq2). \tag{3}
$$

For coefficient sequences, write `A ~_p D` for partial synchronization. The
closure theorem of Hu--Wang--Zhao--Zhao says that partial synchronization is
preserved by common convolution and by nonnegative sums of pairwise
partially synchronized log-concave sequences; such sums are log-concave.

Consequently, the following four base facts suffice:

$$
B,C\text{ are log-concave},\qquad C\sim_pB,qquad xC\sim_pB. \tag{4}
$$

Indeed, `F_1=B+xC`. The first relation needed by the recurrence,
`F_1 ~_p B`, follows by adding `B ~_p B` and `xC ~_p B`. For the second,
common convolution turns `C ~_p B` into `xC ~_p xB`, and adding this to
`B ~_p xB` gives `F_1 ~_p xB`. The same two-relation induction used for the
unit-arm connector theorem now proves that every `F_t` is log-concave.

Thus all dependence on the unbounded connector has been removed: for a fixed
pair `(a,b)`, only the finite base conditions (4) need verification.

## Exact finite certificate

The verifier enumerates every arm multiset as an integer partition. Through
total pendant weight 24 there are:

- 7,338 rooted arm multisets; and
- 163,523 unordered pairs `(a,b)` with total weight at most 24.

For every pair, using integer coefficients only, it verifies all four facts
in (4). The canonical stream consisting of `(a,b,B,C)` for every checked pair
has SHA-256 digest

```text
25c37f1435984d3bbdb275424ea4f6f5766c33a7c6e792a67700a90dae91efb1
```

An independent tree-DP audit checks (1) on 1,115 trees with total pendant
weight at most 8 and `0 <= t <= 4`, and checks the recurrence in 669 of those
cases.

Replay:

```text
python3 verify_c2_bounded_pendant_core_20260808.py
python3 -m unittest test_c2_bounded_pendant_core.py
```

Certificate:

```text
results/c2_bounded_pendant_core_20260808.json
```

The finite enumeration is the proof of the bounded base statement; it is not
evidence for unenumerated pendant weights. The extension from those exact
base cases to all connector lengths is algebraic.

## New universal target

The same computation suggests the following theorem-strength target without
the bound in (4): for arbitrary pendant-arm multisets, `B` and `C` are
log-concave and satisfy `C ~_p B` and `xC ~_p B`. Proving it would establish
log-concavity for every tree with at most two branch vertices.

Ordinary synchronization is too strong: it already fails for arm vectors
`()`, `(1,3)`, even though both required partial-synchronization relations
hold. Any universal proof must therefore work at the partial-synchronization
or aggregate Turán-gap level.

Weak synchronization is not a sufficient replacement. Hu--Wang--Zhao--Zhao
prove that partial synchronization is strictly stronger than weak
synchronization, and their support-reduction criterion still requires the
full two-index partial-synchronization inequalities on the overlapping
support. Thus checking only the diagonal (`m=n`) cannot establish the
universal base invariant.

## Asymptotically tight family

An exact profile through total pendant weight 20 evaluates 9,105,016
two-index inequalities. The tightest strict inequality occurs for arm vectors
`(1,1)` and `(9,9)`, in the diagonal of `xC ~_p B` at the top degree. Its
required sides are `64 >= 58`, a margin of 6.

This is the `h=4` member of an algebraically tractable family. Put the arm
vectors equal to `(1,1)` and `(2h+1,2h+1)`, and let `N=2h+4` be the common
top degree of `B` and `C`. Directly from
`[x^k]P_n=binom(n+1-k,k)` one obtains

$$
B_N=C_N=1,
$$

and

$$
B_{N-1}=2h^2+5h+6,
\qquad
C_{N-1}=h^2+3h+4.
$$

The top diagonal partial-synchronization inequality for `xC` and `B` is
therefore

$$
2C_{N-1}B_N\geq C_NB_{N-1},
$$

with exact margin

$$
2C_{N-1}-B_{N-1}=h+2>0. \tag{5}
$$

Its relative margin tends to zero. Thus this family explains the observed
bottleneck and supplies an asymptotically sharp target for any universal
proof; it does not refute the invariant.

`verify_c2_bottleneck_family_20260808.py` checks the coefficient identities
in (5) and, as separate finite evidence, verifies log-concavity and both full
base partial-synchronization relations for `1 <= m <= 500`. The latter finite
check is not an all-`m` proof. Its report is
`results/c2_bottleneck_family_20260808.json`.

No prior-art novelty claim is made for the bounded theorem without a fuller
authenticated search.
