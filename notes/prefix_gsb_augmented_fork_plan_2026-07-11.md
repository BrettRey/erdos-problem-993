# Prefix GSB: status repair and augmented-fork prove-or-refute cycle
## Outcome
The cycle ended in the **refuted** branch. An exact last-prefix witness on a
center joined to five copies of `G(60,18)` has a negative `(AF_lambda)` gap at
`lambda=0` and negative slope, so every `lambda in [0,1]` fails. The frozen
packet is `notes/d25_augmented_fork_obstruction_packet_2026-07-11.md`. Aggregate
`qfar<=0` remains separate and parked; no broad continuation was opened.

## Objective
Complete one bounded research cycle with two deliverables:

1. make `README.md`, `STATUS.md`, `DECISIONS.md`, and the solve-today registry accurate through D27; and

2. freeze and seriously test the strongest surviving D25 proposition, returning either a proof, an exact prefix counterexample, or a precise structural no-go.


This cycle does not authorize a manuscript change or a solution claim.
## Frozen candidate
Fix a tree `T`, a rank `r` in the prefix window, and put `N=i_r`. Let `a_v` be the number of independent `r`-sets to which `v` is addable. For a center `c`, let `J_c` count, over unordered pairs of distinct neighbours of `c`, the independent `r`-sets to which both neighbours are addable. Define

```text
Q_c = 2 (N J_c - sum_{u<v in N(c)} a_u a_v),
B_c = sum_{u in N(c)} a_u^2/deg(u) + a_c sum_{u in N(c)} a_u,
S_c = sum_{u in N(c)} a_u/deg(u).
```

The candidate is that there is one universal rational `lambda in [0,1]`—with `lambda=1/2` the first frozen target—such that every center satisfies

```text
Q_c <= B_c + N (lambda a_c + (1-lambda) S_c).       (AF_lambda)
```

Summing these local budgets gives exactly the right-hand side of the global prefix-GSB inequality, but `sum_c Q_c=2q2` contains only the distance-two covariance. Thus `(AF_lambda)` implies prefix GSB only together with the separate aggregate-far statement `qfar<=0`. The summation identity and this remaining dependency will be rederived and checked before proof work; neither is being assumed from a script comment.
## Load-bearing assumptions and falsifiers
| Assumption | What would falsify it | First test |
|---|---|---|
| One fixed `lambda` works for all centers and trees in the prefix. | Exact lower and upper constraints have empty intersection, or one center is impossible for every `lambda`. | Re-run the exact exhaustive constraint intersection and preserve extremal witnesses. |
| `lambda=1/2` is a viable canonical target. | An exact negative half-gap in the prefix. | Exhaustive small trees, structured hard compositions, and exact replay of the tightest diagnostic case. |
| Summing `(AF_lambda)` supplies the distance-two part of prefix GSB, with `qfar<=0` as an independent dependency. | A missing multiplicity, degree-zero case, or residual term in the global identity. | Independent algebraic derivation plus exact coefficient checks. |
| The local claim is genuinely simpler than GSB. | Its proof requires GSB, an equivalent global cut statement, or another conjecture of the same strength. | Dependency audit of every proposed proof step. |
| Current documentation captures the whole July 11 round. | D19–D27 artifacts contain material conclusions absent from the status surfaces. | Audit every D19–D27 certificate and script header before editing prose. |

The unstated assumption to guard against is that extensive finite survival indicates a universal local inequality. It does not: computations are used to falsify and to identify extremal structure, never as proof.
## Execution stages
### 1. Repair the status surfaces
- Integrate D19–D27 into the registry and the July 11 summaries.

- Narrow the rank-four conclusion: the order-13 pair blocks exact reconstruction from the recorded radius-two moments; it does not establish that every possible rank-four proof requires radius-three data.

- Distinguish false local allocations from the still-live augmented-fork candidate and distinguish finite evidence from theorem-level claims.

- Keep `paper/main_v2.tex` unchanged.

### 2. Freeze the mathematical packet
- Derive `a_v`, `J_c`, `Q_c`, `B_c`, and `S_c` directly from independent-set counting.

- Prove the global summation identity on paper, replay it on exact trees, and state explicitly that `qfar<=0` is still required.

- State the precise prefix range, boundary cases, and whether a universal interval of `lambda` or specifically `lambda=1/2` is required.

### 3. Adversarial falsification
- Re-run exact exhaustive trees with rational arithmetic and record the intersected `lambda` interval and extremal centers.

- Stress the known 9,418-vertex hard composition and nearby parent/child types, multiplicities, and prefix ranks.

- Convert the tightest floating-point diagnostic into an exact certificate before treating its sign as evidence.

- Search specifically for incompatible lower/upper constraints, not merely a negative half-gap.


An exact prefix failure ends the proof attempt and becomes the frozen result.
### 4. Bounded proof attacks
Only the following mechanisms are in scope:

1. a branch-factorization or generating-function proof at one center;

2. a direct injection or weighted Cauchy–Schwarz argument that explains the inverse-degree terms; and

3. an extremal reduction showing why a fixed `lambda`—ideally `1/2`—controls all branch configurations.


Every intermediate lemma will be tested against the existing D25 local counterexamples. A route that invokes the desired global GSB inequality or a min-cut condition equivalent to it is recorded as circular and stopped.
### 5. Decision checkpoint and closeout
The cycle ends with exactly one of:

- **proved:** a complete universal proof of `(AF_lambda)`, together with the exact statement that it closes the distance-two budget but leaves `qfar<=0` open;

- **refuted:** an exact tree/rank/center certificate; or

- **blocked:** a frozen statement of the smallest missing lemma together with evidence that each permitted proof mechanism either fails or is circular.


No new broad approach family is opened during this cycle. Final verification will include targeted certificate replays, `python3 test_all.py`, compilation of changed Python files, and `git diff --check`.
## Intended repository result
- corrected `README.md`, `STATUS.md`, `DECISIONS.md`, and solve-today registry;

- one focused augmented-fork theorem/obstruction packet;

- only the smallest scripts or result certificates needed to replay a new mathematical conclusion;

- no changes to the manuscript and no public-status change.
