# The three-fold symmetric subcubic slice: 11,892 LC failures and the unequal-arm correction
<!-- SUMMARY: The 3-fold symmetric sweep finds 11,892 LC failures through n=67, all top-corner and none non-unimodal; its gap=1 pattern is an observed condition only within that finite symmetric slice. The exhaustive n=36 census finds 15 unequal-arm failures, including 9 with gap vector (0,0,0), so the gap condition does not extend to unequal arms. A 921-tree witness-seeded sweep reproduces the exact n=36 set and finds 126 LC failures, no unimodality alarms · status: finite slices mapped; infinite-family and mechanism claims open · updated: 2026-08-26 -->

Follow-on from `matroid_intersection_number_lane_2026-08-26.md` §2.5 and §6.2. The
first subcubic log-concavity failure found by this census ($n = 34$) turned out to be a centre
with three pairwise-isomorphic 11-vertex arms and an automorphism group of order
48. This note sweeps that family.

## 1. Method

A three-fold symmetric subcubic tree is a centre joined to three copies of one
rooted arm. The centre then has degree 3, so every arm vertex has at most two
children: arms are exactly the unordered rooted trees with $\le 2$ children per
node, counted by the Wedderburn-Etherington numbers. Writing $F$ for the arm
polynomial with its root excluded and $G$ with the root included,

$$I(T) = (F+G)^3 + x\,F^3,$$

since the centre is either excluded (each arm free) or included (each arm root
excluded). Exact integer arithmetic throughout;
`scripts/symmetric_subcubic_sweep_20260826.py`.

**Reach.** All arms up to 22 vertices is $\approx 5.4$ million trees reaching
$n = 67$. The exhaustive census ladder is at $n = 36$ with $2.6 \times 10^{10}$
trees and grows about 2.3x per rung, so the symmetric sweep buys 31 vertices of
reach for roughly five orders of magnitude less work ~-- justified only because
the one known witness sits in this family.

**Positive control.** The sweep must re-find the $n=34$ witness at arm size 11,
matching its stored polynomial exactly. It does; the run aborts otherwise.

## 2. A bug worth recording

The first version of the two-fold construction multiplied the subtracted term by
$x$, and reported **three NON-UNIMODAL alarms at $n = 24$** ~-- which, taken at
face value, would have been counterexamples to Erdős Problem 993.

They were artifacts. The polynomials ended $\dots, 174, 0, 1$: an independent set
of size 14 with none of size 13, which is impossible, since deleting a vertex
from an independent set leaves one. The construction, not the trees, was wrong.

The script now validates every polynomial it computes ~-- $i_0 = 1$, $i_1 = n$, no
negative coefficients, **no internal zeros** ~-- and treats any violation as a code
error. All three checks are free, and the internal-zero one catches exactly this
class of bug instantly. With validation on, the alarm count is zero.

The lesson generalizes past this script: a search whose success condition is
"find an anomaly" will report its own bugs as discoveries, so the cheap
structural invariants of the object being computed should be asserted on every
instance, not spot-checked.

## 3. Results

**11,892 log-concavity failures**, arms 1..22 (reaching $n = 67$), all in the
three-fold family. The failure-bearing rungs are:

| $n$ | arms swept | failures | rate | break $k$ | depth (dist) |
|---|---|---|---|---|---|
| 34 | 451 | 1 | 0.22% | 17 | 5 |
| 40 | 2,179 | 6 | 0.28% | 20 | 6 |
| 46 | 10,905 | 19 | 0.17% | 23 | 7 |
| 49 | 24,631 | 34 | 0.14% | 26 | 8 |
| 52 | 56,011 | 96 | 0.17% | 26 | 8 |
| 55 | 127,912 | 198 | 0.15% | 29 | 9 |
| 58 | 293,547 | 564 | 0.19% | 29, 32 | 9, 10 |
| 61 | 676,157 | 1,185 | 0.18% | 32 | 10 |
| 64 | 1,563,372 | 3,019 | 0.19% | 35 | 11 |
| 67 | 3,626,149 | 6,770 | 0.19% | 35, 38 | 11, 12 |

Three things stand out.

1. **The observed failure rate stays near $0.2\%$** across the failure-bearing
   rungs. The finite range does not establish a limiting density, and its rate
   does not predict the rate in the full subcubic population.
2. **Every break is at $k = \alpha - 1$**, the top corner, matching the
   general census's finding across all 1,230 known non-subcubic failures. Not one
   exception in 11,892.
3. **Arm sizes 12 and 14 ($n = 37, 43$) produce nothing**, while 11, 13, and
   15 upward all do. The family is not monotone in $n$.

**The two-fold symmetric family (an edge joining two isomorphic arms) produces
zero failures at every size swept.** That finite-range asymmetry is sharp and
unexplained; it motivates studying the degree-3 junction but does not yet
identify a mechanism.

## 4. An observed condition in the symmetric slice

For an arm $A$ with root $r$, define
$$\mathrm{gap}(A) = \alpha(A) - \alpha(A \text{ with } r \text{ forbidden}).$$
$\mathrm{gap} = 1$ says **every maximum independent set of $A$ contains its
root**; $\mathrm{gap} = 0$ says some maximum independent set avoids it.

Over all arms of size 9..18:

| gap | $\#$max ind. sets $L$ | arms | failures | rate |
|---|---|---|---|---|
| 0 | 1 | 8,738 | 0 | 0% |
| 0 | 2 | 12,940 | 0 | 0% |
| 0 | $\ge 3$ | 113,773 | 0 | 0% |
| 1 | 1 | 20,252 | 273 | 1.35% |
| 1 | 2 | 20,261 | 69 | 0.34% |
| 1 | $\ge 3$ | 52,263 | 12 | 0.02% |

> **Observation in the finite symmetric sweep (NOT a theorem ~-- see §4.1).** If some
> maximum independent set of the arm avoids its root, the three-fold symmetric
> tree on that arm appears to have a log-concave independence sequence. Zero
> exceptions in **135,451** gap-0 arms, and the extremal case survives to
> $n = 358$ with a *rising* scaled margin. Instance counts alone would not
> justify believing this; §4.1 is the reason it is worth provisional credit.

And within $\mathrm{gap} = 1$ the failure rate falls by roughly an order of
magnitude per extra maximum independent set, so failures concentrate hard on arms
with a *unique* maximum independent set.

There is an exact top-degree split behind the observation. Let
$a_s = \alpha(A)$ and $a_f$ be the maximum size avoiding the root. Then
$$\alpha(T) = \max(3a_s,\; 1 + 3a_f).$$
If $\mathrm{gap} = 1$, $3a_s = 3a_f + 3 > 1 + 3a_f$, so the top of the sequence
comes from the *centre-excluded* term and $i_{\alpha}(T) = L^3$. If
$\mathrm{gap} = 0$, then $1 + 3a_f > 3a_s$ and the top comes from the
*centre-included* term instead. The two cases put a different summand at the top
corner. This explains why the gap statistic is relevant, but it does not prove
that the gap-0 case cannot break; §8 shows that the corresponding unequal-arm
extension is false. For the symmetric witness the drop is stark:
$i_{16}, i_{17}, i_{18} = 5325, 72, 1$, and $72^2 = 5184 < 5325$.

### 4.1 Adversarial test of the gap condition (Brett doubted it; correctly)

Stating this as a "clean structural law" on the strength of 135,451 exception-free
instances was the wrong move, and it is the same mistake that killed the subcubic
conjecture earlier the same day one rung past its evidence. Instance-counting is
not the test. The test is the **margin**: if gap-0 arms are only barely
log-concave and the slack is shrinking, the condition dies at larger $n$ no
matter how many clean cases are stacked up.

Measuring the normalized Turán margin
$\min_k (i_k^2 - i_{k-1}i_{k+1}) / i_k^2$ in exact rationals over every gap-0
arm at each size:

| $m$ | 8 | 10 | 12 | 14 | 16 | 18 | 20 |
|---|---|---|---|---|---|---|---|
| min margin | .245 | .210 | .154 | .147 | .134 | .108 | .105 |

**The measured margin decreases with $m$**, on a scale compatible with $c/m$.
That could be erosion toward failure, so the useful next check is to identify
the observed extremal arm rather than rely on the aggregate count.

**The observed extremal gap-0 arms form one explicit family through $m=20$.** Writing
$\mathrm{unit}(X) = ((), (X,))$ (a spine vertex carrying a pendant leaf, then
continuing), the observed worst gap-0 arm at the verified sizes $m = 3j+2$ is
$\mathrm{unit}^j(P_2)$ ~-- a caterpillar-like spine with a pendant leaf at every
third vertex. Verified as extremal at $m = 11, 14, 17, 20$ by exhaustive
comparison against all gap-0 arms of those sizes.

Pushing that family well past the sweep, to $j = 39$, i.e. arm size 119 and
$n = 358$, about ten times the reach of the exhaustive sweep:

| $j$ | 4 | 8 | 16 | 24 | 32 | 39 |
|---|---|---|---|---|---|---|
| $n$ | 43 | 79 | 151 | 223 | 295 | 358 |
| margin | .1466 | .0813 | .0430 | .0293 | .0222 | .0183 |
| $m \cdot$ margin | 2.052 | 2.115 | 2.152 | 2.166 | 2.173 | 2.177 |

The margin stays positive throughout this range, and
$m \cdot \text{margin}$ rises monotonically to 2.177. This is evidence against
an imminent crossing along this named family, not a limit calculation.

**Status after the attack.** Within the symmetric slice, the gap-0 observation
has a named extremal candidate followed ten times further than the exhaustive
data, with no crossing. It is still **not a theorem**, and one gap is explicit: the
extremal family was verified extremal only for $m \le 20$. For larger $m$ the
true worst gap-0 arm could belong to a different family with worse behaviour,
and nothing here rules that out. What the test bought is a precise target ~--
anything refuting the condition has to beat $\mathrm{unit}^j(P_2)$.

The restricted symmetric statement is definition-complete enough for a
formalization packet, but §8 removes any reason to expect the same statement for
unequal arms.

## 5. Finite negative evidence through $n=67$

**All 11,892 observed breaks lie at $k = \alpha-1$, and not one lies below
$\lceil (2\alpha-1)/3 \rceil$.** \textcite{levit2006} proved $i_k$ is strictly
*decreasing* for $k \ge \lceil (2\alpha-1)/3\rceil$. A unimodality failure is a
valley, and a valley cannot occur inside a strictly decreasing run.

Consequently, none of the 11,892 recorded failures is a counterexample to Erdős
Problem 993. Across the tested range, break depth grows from 5 at $n=34$ to 12
at $n=67$, moving away from the shallow-break frontier, which needs
$\mathrm{dist}<0$ to threaten unimodality.

This is a reason to deprioritize the symmetric slice for shallow-break searches,
not a proof about untested arm sizes. A future symmetric tree could in principle
develop a break below the decreasing zone or a separate valley even if it also
has a top-corner log-concavity break.

## 6. What is still open here

- **Prove or refute the gap condition within full three-fold symmetry.** The
  $\alpha(T)$ case split above is the start; the unequal-arm extension is already
  refuted in §8.
- **Explain the two-fold/three-fold asymmetry.** Zero failures on one side,
  0.2% on the other, with no explanation offered.
- **Why nothing at arm sizes 12 and 14.**
- What distinguishes the gap-0 unequal-arm failures from the clean gap-0
  symmetric cases.


## 7. Lead opened and killed the same evening: symmetry is not enrichment

All **21** known non-subcubic failures at $n = 26$ and $n = 28$ have nontrivial
automorphism groups, of orders 48 through 82,944. Not one is asymmetric. That
looked like it might generalize the whole approach past subcubic.

**The control refutes it.** Almost all trees are symmetric:

| $n$ | trees | with a nontrivial automorphism | rate |
|---|---|---|---|
| 12 | 551 | 522 | 94.74% |
| 14 | 3,159 | 3,020 | 95.60% |
| 16 | 19,320 | 18,653 | 96.55% |
| 18 | 123,867 | 120,623 | 97.38% |

At a base rate near 97% and rising, seeing 21 of 21 symmetric has probability
about $0.97^{21} \approx 0.53$. It is the single most likely outcome, and
carries no information at all. **The lead is dead**: "LC failures are symmetric"
is a fact about trees, not about failures.

Two things survive it, and they should not be confused with the dead claim:

- The $n=34$ witness's *specific* structure ~-- three pairwise-isomorphic arms and
  an automorphism group of order 48 ~-- is still what made the family sweep in §3
  possible, and that sweep found 11,892 failures. The construction was
  productive whether or not symmetry is statistically special.
- Whether failures have *unusually large* automorphism groups relative to the
  base distribution is a different and untested question. Not pursued.

**Method note.** The first version of this control used
`GraphMatcher.isomorphisms_iter()`, which enumerates the automorphism group; a
star on $n$ vertices has $(n-1)!$ automorphisms, so it hung indefinitely and had
to be killed. Deciding whether *any* nontrivial automorphism exists needs no
enumeration: root at the centre and compare AHU canonical forms of sibling
subtrees. That runs the full $n \le 18$ census in seconds. The distinction
between deciding a property and enumerating its witnesses is worth remembering
next time a symmetry question comes up.


## 8. The census finishes n=36 and corrects this note's premise

The exhaustive subcubic census completed its last rung after this note was
drafted: $n = 36$, **25,664,800,714 trees** (matching OEIS A000672 exactly),
**15 log-concavity failures**, no unimodality alarms. The subcubic sequence now
reads

| $n$ | 33 | 34 | 35 | 36 |
|---|---|---|---|---|
| failures | 0 | 1 | 0 | **15** |

Every one of the 15 has $\alpha = 19$, a lone break at $k = 18 = \alpha - 1$,
and $\mathrm{dist} = 5$; 6 to 9 branch vertices, 8 to 11 leaves. Independent
replay from the parent arrays confirms the stored polynomials exactly and finds
no duplicate trees:

```bash
python3 scripts/verify_subcubic_n36_failures_20260826.py
```

**They are all unicentral with a degree-3 centre and three rooted arms.** Five
have two distinct rooted-arm shapes and ten have three; none has three identical
arms. Full three-fold symmetry is arithmetically impossible at $n=36$, since it
requires $n=3m+1$. This places every current witness in the broader unequal-arm
class and justifies testing that class. It does not establish that a degree-3
centre is an enriched or causal mechanism; that would require a base-rate
control over the full census.

The arm data give two immediate corrections.

1. **The symmetric gap condition does not survive unequal arms.** Put
   $\delta_i=a_s(A_i)-a_f(A_i)\in\{0,1\}$. Then
   $$
   \alpha(T)=\sum_i a_f(A_i)+\max\left(\sum_i\delta_i,1\right).
   $$
   Six witnesses have gap vector $(1,1,1)$, but the other nine have
   $(0,0,0)$. Thus "some arm has a maximum independent set avoiding its root"
   does not prevent a top-corner break outside full symmetry.
2. **The certificates do not support an arm-size 11--12 restriction.** Their
   component sizes range from 2 to 25, including triples $(2,8,25)$,
   $(4,6,25)$, and $(6,6,23)$. The safe rescope is by the observed rooted-arm
   shapes, not by their average size.

The corrected perturbation script therefore extracts the 17 distinct rooted
arms in the certificates and tests every unordered triple, deduplicating
unrooted trees. Its default run examines 969 arm multisets representing 921
distinct trees, reproduces the exact 15-tree $n=36$ failure set, and finds 126 LC
failures total through $n=76$, with no unimodality alarms:

```bash
python3 scripts/witness_perturbation_20260826.py
```

All 126 breaks remain in the Levit--Mandrescu decreasing zone; their depths range
from 5 to 10. Together with the symmetric sweep, this is finite negative evidence
for #993 and a reason to deprioritize these seed shapes for shallow breaks. It is
not a structural exclusion of either the infinite symmetric family or the full
three-arm class.

Finally, 15 failures among 25,664,800,714 trees is a rate of about
$5.84\times10^{-10}$, or one per 1.71 billion trees. That is 2.8 times the
$n=34$ rate, but the sequence $0,1,0,15$ is too short and parity-sensitive to
support a density trend; the failures remain extremely sparse.
