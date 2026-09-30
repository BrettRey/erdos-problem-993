# LP search for a local discharging certificate of central log-concavity
<!-- SUMMARY: An LP over 1.8M exact constraints finds simple one-hop discharging rules certifying window log-concavity; the constant rule "keep 2/3, give 1/3 equally to neighbours" (S23) holds for every tree n<=26 exhaustively (about 448M trees) and on all hub families incl. the R1 killer; n=27-28 census and an adversarial search are running · status: in progress · updated: 2026-09-30 -->

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

The tightest normalized S23 load at n = 24–26 is about -0.030. Raw shard
output is in `raw/`; the progress log is `s23_progress.log`.

**Running:** the exhaustive census for n = 27–28, where the pointwise lemma
first fails, and an adversarial search in `adversary/` trying to shrink the
feasible `sigma` interval to empty and to break S23, up to n of about 300.

## Status

S23 is a deterministic candidate lemma with no size threshold. If it holds
for every tree, #993 follows. It was fitted using hub families that had
already killed R1, and it has survived out-of-sample exhaustive testing
through n = 26. No proof mechanism is known yet.
