# The Vatter-order certificate for window log-concavity
<!-- SUMMARY: Follow-up to Vatter (arXiv:2608.22147): two leaves-first vertex orders certify log-concavity at every k on all 1,346,021 trees with 4<=n<=20, on every recorded killer of earlier deterministic routes (n=28 to 1293), under random relabelings, and on all 11,878 members of the 30 Sep adversary's family zoo (n<=500, zero window failures), and 25-minute annealing runs (no violation); graded PLAUSIBLE (survived every known killer; not proved; not exact above the window) · status: lead recorded, nothing queued · updated: 2026-10-03 -->

**Read this first.** On #993, small trees predict nothing. On 30 Sep, every
deterministic candidate held through n ≤ 27 and failed only on constructed
trees with hundreds or thousands of vertices: the pointwise free-count lemma
at n = 28, the neighbourhood average R1 at n = 237, and the constant
discharging rules at n = 451 and n = 1293 (`runs/attack-20260930/SUMMARY.md`).
The evidence below that counts is the killers, the adversary's family zoo
and the annealing search. The exhaustive census is a sanity check only.

## The certificate

Fix a vertex order of a tree T and classify independent sets by their least
vertex:

i_{k+1}(T) = Σ_v i_k(T[L(v)]),  where L(v) = later vertices not adjacent to v.

Hence ρ_{k+1}(T) is the i_{k−1}(L(v))-weighted average of ρ_k(L(v)), where
ρ_k = i_k / i_{k−1}. So T is log-concave at k whenever every tail with
i_{k−1}(L(v)) > 0 satisfies

ρ_k(L(v)) ≤ ρ_k(T).  (V_k)

An order satisfying (V_k) is a *Vatter order at k*. This transplants
Vatter's first-letter proof of Chase's theorem (arXiv:2608.22147), where
Claim 2 supplies (V_k) for words. For #993 only the window
[⌈n/4⌉, q], q = ⌈(2α−1)/3⌉, is needed. Zhang–Li's Lemma 3.1 would shorten
it further, to [⌈n/4⌉, ⌊(n+2α)/6⌋].

**Relation to the 30 Sep routes.**
- The pointwise free-count lemma compares ρ_k(T − N[v]) with
  ((k+1)/k)·ρ_k(T), over every v at once and without an order.
- (V_k) is stricter at the first vertex: its tail is all of T − N[v], and
  there is no (k+1)/k slack.
- After that, tails shrink. A vertex whose T − N[v] is rich (a hub or a
  centre) can be placed late, where its tail is small.

## Tools

`scripts/vatter_lib_20261003.py` (python-flint, exact) and
`scripts/probe_vatter_order_certificate_20261002.py`. Checks run on every
use:
- **Tree DP.** The library agrees with `indpoly.independence_poly` on every
  tree with n ≤ 12.
- **Closed forms.** It reproduces the 30 Sep closed forms (`symrows.H`,
  `symrows.MSH`) for I(T) and for every I(T − N[v]).
- **Identity.** `check_order` asserts i_{k+1}(T) = Σ_v i_k(L(v)) at every k
  whenever it runs with the assertion on.

## Results

**1. First-vertex filter (necessary condition).** Every Vatter order needs
some v with ρ_k(T − N[v]) ≤ ρ_k(T). This filter never binds. Almost every
vertex qualifies across the whole window on every tree tested:

| tree | qualifying vertices |
|---|---|
| H(9,2) | 27 of 28 |
| MSH(4×9,4×10) | 229 of 237 |
| H(75,5) | 450 of 451 |

The exceptions are the hubs and the centre.

**2. Order rules.** Centre-first orders fail on the hub families:
- BFS from the centre fails at 4 of 7 window k on H(9,2), at 21 of 21 on
  H(20,3), and at 26 of 48 on MSH(4×9,4×10).
- Degree-descending orders fail similarly.

Two leaves-first rules pass:
- **revbfs:** the reverse of BFS from a centre.
- **degasc:** ascending degree, ties broken by vertex index.

**3. Killers, every k (not just the window), both rules pass:**

| tree | n | α | window | time |
|---|---|---|---|---|
| H(9,2) (pointwise-lemma killer) | 28 | 19 | [7,13] | <0.1 s |
| H(20,3) | 81 | 61 | [21,41] | <0.1 s |
| MSH(4×9,4×10) (R1 killer) | 237 | 160 | [60,107] | 0.2 s |
| S(8,10,2) = MSH(8;10,2) | 249 | 168 | [63,112] | 0.1 s |
| S(9,9,2) = MSH(9;9,2) | 253 | 171 | [64,114] | 0.1 s |
| H(75,5) (2/3-rule killer) | 451 | 376 | [113,251] | 0.3 s |
| MSH(38;11,2) (excludes every constant rule) | 1293 | 874 | [324,583] | 2-5 s |

**4. Tie-breaking robustness.** Under 12 random relabelings of each of
H(9,2), H(20,3), MSH(4×9,4×10), S(8,10,2), H(75,5), a 30×3 caterpillar
and S(5^12), neither rule failed at any k. The verdicts do not depend on
which centre, adjacency order or index tie-break is used.

**5. Exhaustive census, every k.** Both rules pass on every tree with
4 ≤ n ≤ 20, 1,346,021 trees. Per n from 16 to 20 the counts were 19,320,
48,629, 123,867, 317,955 and 823,065. This is a sanity check (see the
caveat at the top).

**6. The two order-26 LC failures.** Both rules pass the whole window [7,9]
and fail only at k = 13, where log-concavity itself fails. Any certificate
must fail there.

**7. Adversary's family zoo** (`runs/attack-20260930/lp-discharge/adversary/families.py`,
all 11,878 specs with n ≤ 500: uniform and mixed MSH, mixed leaf counts,
centre decorations, subdivided edges, nesting, two centres, caterpillar
spines, single hub-stars). **All 11,878 members (n from 24 to 500) pass the whole window under both
rules: zero window failures.** By family: MSHmix 6,830, H 1,955, MSH 1,059,
twocentre 503, nest 406, MSH+cp 344, MSH+cl 340, catspine 194, MSH+sub
168, MSHsvar 79.

Outside the window, the rules fail at some k on 2,821 members (revbfs) and 2,847 (degasc), about 24%. In a
sample of 400 such failures, 385 sat at a k where log-concavity itself fails
(near the top of MSH families with s = 1). 15 sat at an LC-true k just
below the top (k ≈ 0.87α–0.9α). So the certificate is **not** exact. It can
fail where log-concavity holds, though only in the upper tail, which
unimodality does not need.

**8. Annealing adversary.** This reuses the 30 Sep mutation operator
(`sa_r1.mutate`) and seed families, with the objective swapped to the
worst normalised slack max_{v,k∈W} ρ_k(L(v))/ρ_k(T) − 1. Two runs of
25 minutes each, n ∈ [30, 260], one per rule. **No violation.**

| rule | evaluations | restarts | best slack | at | tree |
|---|---|---|---|---|---|
| revbfs | 17,761 | 45 | −0.00642 | first tail, k = 62 = ⌈n/4⌉ | n = 247 |
| degasc | 14,953 | 38 | −0.00644 | first tail, k = 62 | n = 248 |

The seeds already reached −0.00655. In both runs the annealer converged on
one hub of degree about 40.

**9. The binding case scales safely.** The tightest constraint is always
the first vertex, a leaf, with tail T − N[leaf], at the bottom of the
window. On hub-stars H(m,s), I(T) and I(T − N[leaf]) have closed forms, so
I scanned exactly up to n ≈ 27,000. The slack stays negative and shrinks
like 1/n, with n·slack converging to a negative constant.

| s | n·slack (large m) |
|---|---|
| 1 | about −4.0 |
| 2 | about −3.15 |
| 3 | about −2.28 |
| 5 | about −1.62 |
| 8 | about −1.35 |

With m fixed and s growing, n·slack tends to −4/3 for every m tested (1 to
20). That is, the margin is about 4/(3n), against the star's central LC
margin of about 4/n.

## Grade

**Plausible.** This is the 30 Sep attack's top grade short of proved.

**What it survived:**
- every tree that killed an earlier deterministic route (n = 28, 237, 249,
  253, 451, 1293), at every k and under relabeling;
- the full 11,878-member family zoo built to kill those routes;
- 32,700 annealing evaluations aimed at its own slack;
- exhaustive census to n = 20.

Its tightest constraint has a uniform scale-free margin on the extremal
family. The earlier routes had no such margin: the free-count lemma failed
by a constant factor, and the degree rules survived only with zero slack.

**What it is not:**
- not proved;
- not exact above the window (15 LC-true failures in the upper tail);
- not tested against an adversary designed with this certificate in mind
  beyond one annealing pass.

The family zoo was built against other routes, so passing it is weaker
evidence than surviving a dedicated search would be.

**What would refute it:** a tree where every leaves-first order puts a tail
over the line inside the window. The natural places to look are trees whose
first leaf's T − N[leaf] keeps a dense core, which is the opposite of a hub.

**What would promote it:** a proof that the first leaf's tail satisfies (V_k)
on the window. Later tails are smaller, and an inductive argument might
cover them. Not attempted.

## Context

Zhang–Li (Zenodo 22999166, 27 Sep) claim a full proof of #993. It is under
local review (`notes/zhang_li_2026_forest_unimodality_2026-10-03.md`). If it
stands, the Vatter certificate is no longer needed for #993. Its interest
would then be a short deterministic certificate for central log-concavity of
trees, plausibly provable by a direct argument about leaves-first tails.
That is a separate question, and nothing is queued for it.
