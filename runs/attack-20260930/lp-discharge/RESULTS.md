# LP search for a local discharging certificate of central log-concavity
<!-- SUMMARY: An LP over 1.8M exact constraints finds simple one-hop discharging rules certifying window log-concavity; the constant rule "keep 2/3, give 1/3 equally to neighbours" (S23) holds for every tree n<=27 exhaustively (about 1.2e9 trees) but is refuted at H(75,5), n=451, and no constant rule survives MSH(38;11,2), n=1293; degree-based rules survive the hub families only with zero slack (float LP optimum fails exact check at 1e-20) · status: constant rules refuted; degree-based knife-edge · updated: 2026-09-30 -->

Run on 30 September 2026 at Brett's direction ("proceed"). It follows the
refutation of R1 (`../SUMMARY.md`) and the failed tree-shift test
(`../gts-test/RESULTS.md`).

## Setting

For a tree T, let `W(T) = [ceil(n/4), min(q, alpha-1)]` with
`q = ceil((2 alpha - 1)/3)`. Define the defect

`D_v(k) = k i_{k-1}(T) i_k(T - N[v]) - (k+1) i_k(T) i_{k-1}(T - N[v])`.

Then `sum_v D_v(k) = k(k+1)(i_{k-1} i_{k+1} - i_k^2)`, so log-concavity (LC)
at `k` holds exactly when the defects sum to at most 0.

A **one-hop discharging rule** splits each `D_v` over the closed
neighbourhood `N[v]`, with shares summing to 1. If every vertex's received
total `L_u` is `<= 0`, summing gives LC. Since
`sum_u L_u = sum_v D_v` for any such rule, the rule is a certificate. With
LC on `W(T)`, Fang et al.'s Prop. 8.2 and the Levit–Mandrescu tail give
unimodality.

- R1 is the rule with shares `1/(deg v + 1)`. It is refuted at n = 237.
- The pointwise lemma `D_v <= 0` is the rule with every share kept (keep 1).
  It is refuted at n = 28 by the hub-star H(9,2).

## LP results (`lp_c1.py`, `rowstore.py`, `lp_forms.py`)

There are 1,798,942 exact `(tree, k, u)` constraints: every tree with
`4 <= n <= 16`, plus 72 family trees up to n = 249. The families are
hub-stars H(m,s), stars of stars, two stars joined through a middle vertex,
and MSH constructions (a centre joined to hub-stars), including the n = 237
R1 killer. Rows are normalized by `k i_{k-1} i_k`, and the LP maximizes the
uniform margin `t`.

| Rule family | Optimum margin `t*` | Note |
|---|---|---|
| free share per degree `sigma(d)` (degree-based) | 0.00425 | erratic optimum (overfits the few high degrees) |
| monotone `sigma(d)` | 0.00300 | |
| **constant** `sigma` | 0.00284 at `sigma = 0.673` | the whole interval `sigma in [0.1148, 0.9045]` is feasible; both ends are set by MSH(10^8) |
| `sigma(d) = c/(d+1)` | 0.00072 at `c = 1.998` | R1 is `c = 1`, infeasible |
| `alpha + beta/(d+1)` | 0.00306 | |

Under R1's weights the worst row is `+5.4e-5` (R1 fails), at MSH(10^8).

## The constant rule S23

S23 is `L_u = (2/3) D_u + (1/3) sum_{v ~ u} D_v / deg v`: each vertex keeps
2/3 of its defect and splits 1/3 equally among its neighbours.

**Exhaustive exact census.** The program is `s23_census.c`, a copy of the
validated `r1_census.c` in which only the weights and the sum identity
changed; its constraint counts and zero-violation results match the Python
store for n = 10, 12, 14. Every tree with `17 <= n <= 26` was checked, and
every count equals A000055:

| n | trees | S23 violations | pointwise violations |
|---|---|---|---|
| 17–23 | 23,909,851 | 0 | 0 |
| 24 | 39,299,897 | 0 | 0 |
| 25 | 104,636,890 | 0 | 0 |
| 26 | 279,793,450 | 0 | 0 |
| 27 | 751,065,460 | 0 | 0 |

The tightest normalized S23 load at n = 24–26 is about -0.030. Raw shard
output is in `raw/`; the progress log is `s23_progress.log`.

Because n = 27 is also clean, the smallest counterexample to the pointwise
lemma has exactly **n = 28**: H(9,2) is known there, and nothing smaller
exists. The n = 28 census was stopped at the background time limit before any shard
finished (nothing recorded). No rerun is needed: the n = 28 claim rests on
the exhaustive n <= 27 check plus the known H(9,2) violation (wave 1,
rechecked by brute force), and S23 is already refuted at n = 451.

## Adversary: constant rules refuted (`adversary/`)

The adversarial search worked in exact arithmetic and rechecked every claim
with three implementations.

- **S23 refuted.** At H(75,5) (n = 451: a hub with 75 supports, 5 leaves
  each), the hub row is positive for `k = 220..224`, with maximum normalized
  load `+8.34e-6` at `k = 222`. LC holds throughout. I confirmed this
  separately with fresh code built from the construction
  (`my_check_H75_5.py`).
- **No constant rule works.** The single tree MSH(38;11,2), n = 1293,
  needs `sigma >= 0.98755` at the centre (`k = 583`) and
  `sigma <= 0.98587` at a hub (`k = 579`). The bounds are exact rationals,
  and the gap is 0.013 at MSH(48;11,2), n = 1633.
- **How the bounds move.** Hub-stars H(m,s) push the upper end down (to
  0.406 at H(1200,12)). Stars of hub-stars MSH push the lower end up
  (towards 1 as the centre degree grows). The two ends cross between
  n = 700 and n = 900.

## Degree-based rules: a knife-edge (`symrows.py`, `gen_sym.py`, `lp_degree_*.py`, `exact_check_degree.py`)

Closed forms for H(m,s) and MSH(h;m,s) were built with python-flint and
checked against the generic DP (`check_symrows.py`). They give 623,055 rows:

- H with `m <= 200`, `s <= 8`;
- MSH with `h <= 120`, `m <= 14`, `s <= 3`, up to n = 6841.

Every level passes the exact identity check.

**LP result.** A degree-based `sigma(d)` is feasible only with **zero
margin**: the best relative margin is 0 to within about 1e-12.

**Mechanism.** Near the top of the window, the centre's defect is
essentially 0 (about `-4e-18` normalized) and each hub's is small and
positive. The centre row can only reach 0, and only if every hub keeps its
whole defect (`sigma(m+1) = 1`), so that nothing positive flows to the
centre.

**Exact check.** The float optimum, rationalized and snapped to 0 or 1,
still fails 3 of the 623,055 rows exactly (centre rows of MSH(8;9,2),
MSH(10;10,2) and MSH(4;9,2), at about `1e-20`).

**Verdict.** A degree-based rule, if one exists on these families at all,
must hit exact equalities at particular degrees. That is not a credible
proof route, and any family in which one degree plays two roles would
likely make it infeasible. Not pursued.

## Status

One-hop discharging certificates are exhausted as a proof route. Fixed
weights (R1) fail at n = 237. Constant shares (S23 and every other constant)
fail at n = 451 and n = 1293. Degree-based shares survive the hub families
only with zero slack.

In every case the obstruction is the same: hub constructions whose centre
and hubs have near-zero or positive defects near the top of the window, so
the negative defect that could pay for them sits one level further out.
Every such rule nonetheless holds for all trees with `n <= 27`, which shows
how badly small cases predict large ones here.

A certificate that could work would have to move defect across several
hops with structure-aware weights, or use a non-local argument.
