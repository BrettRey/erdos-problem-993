# Radius check: must positive defect travel farther as nesting deepens?
<!-- SUMMARY: Cheap pre-test for a multi-hop certificate search: in 13,444 layered trees (depth 2-6, up to n about 7,000) positive defect always sits at isolated depths adjacent to negative defect, and every positive vertex is payable by its own children and grandchildren with at least 5x surplus (all killers included); no regress, so one-hop failures come from uniform shares, not distance · status: done · updated: 2026-09-30 -->

Run on 30 September 2026 at Brett's direction ("do the cheap check"), before
any larger search for multi-hop certificates.

## Hypothesis tested

One-hop discharging rules failed on hub constructions (`../lp-discharge/`).
If each extra level of nesting pushed positive defect one step farther from
any negative defect, every fixed radius would fail and a multi-hop search
would be pointless. A radius-`r` rule, even one tailored to a single tree,
can only cancel a positive defect against negative defect within distance
`2r`.

## Method

The family is *layered* trees: every vertex at depth `j` has `b_j` children,
which covers hub-stars `[m,s]`, stars of hub-stars `[h,m,s]` and deeper
nestings. Their independence polynomials, and those of `T - N[v]` for a
vertex at each depth, have closed forms (`layered.py`, python-flint). These
match the generic forest DP exactly on six test trees, including trees with
path segments (`b_j = 1`) (`check_layered.py`). Every window level passes the
exact identity check `sum_v D_v = k(k+1)(i_{k-1}i_{k+1} - i_k^2)`.

## Results (all exact integers)

| Check | Result |
|---|---|
| Sign patterns, depths 2–6, 13,444 trees (`scan_layered.py`) | positive defect appears only at **single isolated depths**, always adjacent to a negative depth (distance 1 in every case). The widest run of positive depths is 1. |
| Private payability (`private_mass.py`, `private_mass_scan.py`) | in all 387,969 `(tree, depth, k)` positive cases, each positive vertex's **own** children and grandchildren carry at least 5 times its defect in negative defect. The tightest case is `[9,12,12,2]`, n = 4006. |
| The killers | H(9,2): 518x. H(75,5), the S23 killer: 602x. MSH(8;10,2): 128x. MSH(38;11,2), which excludes every constant rule: 1062x. MSH(96;11,2): 1147x. MSH(120;14,3): 3860x. |
| Small trees | not needed: the pointwise lemma holds for all trees `n <= 27` (`../lp-discharge/`), so there is no positive defect in the window below n = 28. |

**Float trap.** A first max-flow run (`flow_killers.py`, networkx) reported
MSH(38;11,2) and MSH(96;11,2) infeasible. The cause was float overflow: their
raw defects have 585 and 1,480 digits, far beyond about 308 for a float. The
exact private-payability check above replaces that result. Do not use
float-capacity max-flow on these trees.

## Verdict

No regress, at least on layered constructions. A one-hop certificate tailored
to each of these trees always exists, with large slack. So the failure of
uniform one-hop rules comes from the **uniformity** of their shares, not from
distance. Those rules spread a hub's share equally, including to a parent
that cannot pay, while its own descendants have hundreds of times the needed
negative defect.

What this implies for the larger search: more hops are not the issue. The
open question is whether a *uniform, local* rule can tell a vertex's paying
neighbours (typically its children in these constructions) from the ones that
cannot pay. Degree alone is ambiguous whenever degrees coincide across roles,
so shares should probably depend on each neighbour's own neighbourhood: a
radius-1 transfer decided from a radius-2 view.

**Caveat.** Only layered (spherically symmetric) trees were scanned.
Asymmetric constructions were not tested, although every killer found so far
is layered.
