# PARTIALLY SUPERSEDED — USES WRONG FORMULA IN PLACES

**Some sections use the old polynomial identity** (A = Q_uQ_v + xP_uP_v).
The correct identity is A = P_uP_v + xR_uR_v, proved via the
subdivision-contraction identity I(T') = I(T) + x·I(T/e).
The proved lemma A(x) <= (1+x)I(T;x) is valid but uses the old notation.

**See `notes/subdivision_new_findings.md` for the definitive analysis.**

---

# Subdivision preserves unimodality? (Erdos path)

## Statement (target lemma)

Let T be a tree and T' the tree obtained by subdividing any edge uv
with a new vertex w. If the independence sequence of T is unimodal,
then the independence sequence of T' is unimodal.

If true, this implies any *minimal* counterexample (by number of vertices)
has no degree-2 vertices (is homeomorphically irreducible).

## Exact polynomial identity

Let uv be the subdivided edge. Remove uv to split T into components
A (containing u) and B (containing v). Let (P_u, Q_u) be the rooted pair
for A at u, and (P_v, Q_v) the rooted pair for B at v:

  P = polynomial for sets excluding the root,
  Q = polynomial for sets including the root (with the x factor).

Then the independence polynomials satisfy:

  I(T)  = P_u P_v + P_u Q_v + Q_u P_v,
  I(T') = (P_u + Q_u)(P_v + Q_v) + x P_u P_v
        = I(T) + Q_u Q_v + x P_u P_v.

This identity has been checked by brute force on random trees and edges.

## Why this might work

We need to show the added term

  A(x) := Q_u Q_v + x P_u P_v

cannot introduce a new valley when added to I(T).

Empirical observations (checked for all roots on n <= 10 trees and random
trees up to n=50):
  - P_u and B_u := Q_u/x are log-concave in these small tests.
  - Consequently, P_u P_v and Q_u Q_v are log-concave (Hoggar) in those tests.

The challenge is that a sum of log-concave sequences need not be unimodal.
We likely need a structural inequality relating the modes or first differences
of the summands.

### Closure lemma: adjacent sums preserve log-concavity

Let (b_0, ..., b_m) be log-concave with no internal zeros and define
r_k = b_k + b_{k+1} for 0 <= k <= m-1. Then (r_k) is log-concave.
This is a special case of Hoggar's convolution closure, since r is the
convolution of b with (1, 1). A direct proof is:

  r_k^2 - r_{k-1} r_{k+1}
  = (b_k^2 - b_{k-1} b_{k+1})
    + (b_{k+1}^2 - b_k b_{k+2})
    + (b_k b_{k+1} - b_{k-1} b_{k+2}),

where the last term is nonnegative by ratio monotonicity. Hence r is
log-concave.

Implication for rooted pairs: if B_T := Q_T / x is log-concave, then
(1+x) B_T is log-concave. This helps control terms like (1+x) B_u B_v,
but it still leaves the main obstruction: sums of log-concave sequences
can fail unimodality without an additional mode-ordering or difference
dominance lemma.

### Obstruction: B_T is not log-concave in general

The property “B_T is log-concave for every rooted tree” is false.
If F is any tree whose independence polynomial is not log-concave, construct
T by adding a path r--u--w where w is a vertex of F, and root T at r.
Then T - N[r] = F, so B_T(x) = I(F; x), inheriting any log-concavity failure.

Verified from the two n=26 LC-failure trees in `results/analysis_n26.json`:
constructing T as above yields B_T exactly equal to I(F), and B_T is not
log-concave. This blocks any global closure chain starting with “B_T is LC.”

## A plausible sufficient condition

Write i_k for coefficients of I(T) and a_k for coefficients of A(x).
If we can show, for all k in the descent region of i_k,

  (i_{k+1} - i_k) >= - (a_{k+1} - a_k),

then i_k + a_k stays unimodal. This is the same “difference dominance” pattern
used in the broom s >= p proof.

Given the rooted decomposition, a_k is built from P_u P_v and Q_u Q_v,
so the goal is to express their first differences in terms of the rooted
increment sequences and to bound them against the differences of I(T).

## Empirical evidence

No counterexample found by:
  - exhaustive n <= 10 (all trees, all edges),
  - random trees up to n=40 (sampled edges).

Extended checks (edge subdivisions):
  - exhaustive n <= 17: 1,242,936 subdivisions, no counterexample,
  - exhaustive n = 18: 2,105,739 subdivisions, no counterexample,
  - exhaustive n = 19: 5,723,190 subdivisions, no counterexample,
  - random n = 30..200: 3,240 sampled subdivisions, no counterexample.

Additional empirical facts (n <= 10 exhaustive; random larger tests):
  - The added term A(x) = Q_u Q_v + x P_u P_v is log-concave and unimodal.
  - The first descent index of I(T') is never earlier than that of I(T).
  - A(x) is nondecreasing up to the first descent of I(T).
  - The tail inequality ΔA_k <= -ΔI_k for all k >= d(I)+1 held for all
    subdivisions through n = 19.
  - The boundary inequality at k = d(I) is false (see counterexample below).

These suggest that subdivision “delays” the peak rather than creating a valley,
but this is only evidence; not a proof.

## Tail reduction lemma (subdivision)

Let d = d(I(T)) be the first descent index of I(T), and let A(x)=I(T')-I(T).
If

  (1) d(A) >= d, and
  (2) ΔA_k <= -ΔI_k for all k >= d+1,

then I(T') is unimodal.

Sketch: for k <= d-2, (1) implies ΔA_k >= 0 and unimodality of I(T) gives
ΔI_k >= 0, so Δ(I+A)_k >= 0. For k >= d+1, condition (2) yields
Δ(I+A)_k <= 0. The only free index is k=d-1, which does not create a valley
given the sign pattern on either side.

This reduces subdivision preservation to a tail inequality plus mode-ordering
of A (no boundary condition at k=d is required).

## Tail ratio-monotonicity sufficient condition (new)

Define
  p = P_u P_v,
  q = Q_u Q_v,
  r = P_u Q_v + Q_u P_v,
  i_k = p_k + r_k,
  a_k = q_k + p_{k-1},
  j_k = q_k + r_k.

Then
  Δa_k + Δi_k = (p_{k+1} - p_{k-1}) + Δj_k.   (*)

Lemma (sufficient). Assume nonnegative coefficients and:
  (1) the tail ratios ρ^p_t = p_t/p_{t-1} are nonincreasing for t >= d+1,
  (2) the tail ratios ρ^j_t = j_t/j_{t-1} are nonincreasing for t >= d+2,
  (3) p_{d+2} <= p_d,
  (4) j_{d+2} <= j_{d+1}.

Then Δa_k <= -Δi_k for all k >= d+1.

Sketch: (3) implies ρ^p_{d+2} ρ^p_{d+1} <= 1; with (1), this yields
p_{k+1} <= p_{k-1} for all k >= d+1. Condition (4) and (2) give
Δj_k <= 0 for all k >= d+1. Plug into (*) to get Δa_k + Δi_k <= 0.

This replaces the too-strong requirement “j is log-concave everywhere” with
tail ratio monotonicity plus two boundary checks.

Empirical note (TODO verify): in scans through n<=19, boundary checks (3)(4)
and tail dominance held for all edges, while j sometimes failed log-concavity.
So a proof requiring global log-concavity of j is likely too strong.

New ratio-tail scan (n<=14, all trees/all edges, networkx backend) found no
counterexample to any of the four conditions above:
  results/subdivision_ratio_n14.json
with 66,698 edge checks and no failures of (3), (4), or tail ratio monotonicity
for p and j. This supports the feasibility of the ratio-tail route.

Extended scan (n<=16, backend auto) likewise found no failures:
  results/subdivision_ratio_n16.json
with 464,872 edge checks and no failures of (3), (4), or tail ratio monotonicity.

Finite-core ratio check (b0=6, lambda0=5, networkx cores) also found no failures:
  results/finite_core_ratio_b6_l5.json
with 1,387,276 edge checks and no failures of (3), (4), or tail ratio monotonicity
for p and j. This supports the “finite kernel + ratio-tail” closure route.

Finite-core ratio check (b0=7, lambda0=4) likewise found no failures:
  results/finite_core_ratio_b7_l4.json
with 2,226,670 edge checks and no failures of (3), (4), or tail ratio monotonicity.

Finite-core ratio check (b0=8, lambda0=4) likewise found no failures:
  results/finite_core_ratio_b8_l4.json
with 19,752,622 edge checks and no failures of (3), (4), or tail ratio monotonicity.

### Boundary inequality is false (small counterexample)

The boundary condition ΔA_d <= -ΔI_d fails for T = K_{1,3} when subdividing
a leaf edge:

  I(T) = 1 + 4x + 3x^2 + x^3, so d(I)=1 and ΔI_1 = -1.
  A(x) = x + 3x^2 + x^3, so ΔA_1 = 2.

Thus ΔA_d = 2 \nleq 1 = -ΔI_d. This shows the boundary inequality at k=d
cannot be used as a general criterion.

### New inequality candidate (tested n <= 12 exhaustive)

Define i_k = [x^k] I(T) and a_k = [x^k] A(x) for the subdivided edge.
In exhaustive checks for all trees up to n = 12 (networkx enumeration, all edges),
the following coefficientwise bounds held without exception:

  (1) a_k <= i_k + i_{k-1}
  (2) a_k <= i_{k-1} + i_{k+1}

Also, A(x) was unimodal in all these cases.

If (1) can be proved combinatorially (e.g., by an injection from the “new”
independent sets in T' into independent sets in T of sizes k and k-1),
it gives a concrete handle for the first-difference dominance inequality
needed to show I(T) + A is unimodal.

## Lemma (proved): A(x) <= (1+x) I(T;x) coefficientwise

Write Q_u = x R_u and Q_v = x R_v, where R_u and R_v have nonnegative
coefficients. Note that R_u counts independent sets in A - N[u], while
P_u counts independent sets in A - {u}, so coefficientwise R_u <= P_u.
Likewise, R_v <= P_v.

Compute:

  (1+x) I(T) - A
    = (1+x)(P_u P_v + P_u Q_v + Q_u P_v) - (Q_u Q_v + x P_u P_v)
    = P_u P_v + (1+x)(P_u Q_v + Q_u P_v) - Q_u Q_v
    = P_u P_v + x(P_u R_v + R_u P_v)
      + x^2 (P_u R_v + R_u P_v - R_u R_v).

The first two terms are coefficientwise nonnegative. For the x^2 term,
use R_u <= P_u to get (P_u - R_u) R_v >= 0 coefficientwise, hence
P_u R_v - R_u R_v >= 0, and R_u P_v >= 0. Therefore every coefficient
of (1+x) I(T) - A is nonnegative, so A(x) <= (1+x) I(T;x) coefficientwise.

Equivalently, a_k <= i_k + i_{k-1} for all k.

### Consequences

1) The new term A is dominated by a one-step smoothing of I(T).
2) If the descent of I(T) dominates the local ascent of A, then
   I(T) + A is unimodal.

## What remains (proof skeleton)

It suffices to prove any one of the following local dominance statements:

### Path A: First-difference dominance

Show that for all k in the descent region of I(T),

  (i_{k+1} - i_k) >= - (a_{k+1} - a_k).

This implies I(T) + A is unimodal. The coefficientwise bound
A <= (1+x) I(T) suggests this should be true, but it does not by itself
control differences. A viable route is to bound (a_{k+1} - a_k) by a
nonnegative linear combination of nearby i_j differences.

### Path B: Mode ordering

Prove that the first descent index of A(x) is no earlier than that of I(T).
Empirically, A is nondecreasing up to the first descent of I(T).
If this can be proved, then the sum I(T) + A cannot create a new valley.

### Path C: Injection strengthening

The inequality A <= (1+x)I is now proved algebraically. A purely
combinatorial injection from independent sets counted by A_k into those
counted by I_k and I_{k-1} might give the missing “difference control.”
If the injection is edge-local, it would likely yield Path A.

### Path D: Tail ratio monotonicity (new)

Use the sufficient condition above: prove tail ratio monotonicity for
p = P_u P_v (from t >= d+1) and j = q + r (from t >= d+2), plus the two
boundary checks p_{d+2} <= p_d and j_{d+2} <= j_{d+1}. This is weaker than
global log-concavity and is consistent with current empirical scans.

Sketch of a proof route for Path D:
  1) Interpret p, q, r as independence polynomials of explicit forests:
       p = I((A-u) ⊔ (B-v)),
       q = x^2 I((A-N[u]) ⊔ (B-N[v])),
       r = x I((A-u) ⊔ (B-N[v])) + x I((A-N[u]) ⊔ (B-v)).
     Hence j = q + r is the “at least one root chosen” polynomial for the
     forest F=(A ⊔ B) with the uv edge removed.
  2) Prove a tail log-concavity / tail ratio-monotonicity lemma for forests
     that starts earlier than the LM threshold. For example, show that if a
     forest is leaf-heavy (has a hub with s large compared to the core), then
     its coefficient ratios are nonincreasing from k >= d(I)+1. This should
     be approachable with the leaf-attachment MBI bounds.
  3) Reduce the remaining “leaf-light” cases to a finite kernel class, then
     verify Path D’s boundary checks and ratio monotonicity by enumeration.

This would complete subdivision preservation without requiring global
log-concavity of rooted polynomials. The missing step is a rigorous tail
ratio-monotonicity lemma for forest polynomials that aligns with d(I)+1.

#### Lemma attempt A (PF∞ core, conditional)

If A is PF∞ (real-rooted with nonpositive zeros), then (1+x)^s A is PF∞ for all s
(Hoggar), hence its coefficient ratios are globally nonincreasing. This yields
tail ratio monotonicity for p whenever the core polynomial for p is PF∞. A similar
statement applies to the fixed-core form of j when its core is PF∞. This is too
strong for general trees but may cover leaf-heavy regimes where the effective core
is a small PF∞ polynomial (paths, stars, and their small unions).

#### Lemma attempt B (binomial mixtures, open)

Let U_s(x) = (1+x)^s A(x), with A(x)=∑ a_j x^j, a_j≥0. Define r_k = U_{s,k+1}/U_{s,k}.
We need r_k to be nonincreasing for k ≥ k0 (ideally k0 = d(I)+1).

Observations:
  - r_k is a weighted average of r_j(k) = (s-k+j)/(k+1-j), with weights
    w_j(k) ∝ a_j C(s,k-j).
  - For fixed j, r_j(k) is strictly decreasing in k.
  - The weight family w_j(k) is TP2/MLR in (k,j) because the binomial matrix
    C(s,k-j) is totally positive.
  - Update formula: with weights w_j(k) at level k,
      r_k = E_k[r_J(k)],
    and
      r_{k+1} = E_k[r_J(k) r_J(k+1)] / E_k[r_J(k)].
    Also w_j(k+1) ∝ w_j(k) r_j(k), so w(k+1) is a size-biased tilt of w(k).

Plausible route: show that for k in the right tail, the decrease of r_j(k) in k
dominates the upward shift of the weights w_j(k), yielding r_{k+1} ≤ r_k. This is
an “expectation of a decreasing function under an MLR family” type statement, but
the monotonicity direction is not immediate. A clean TP2 inequality or a discrete
Chebyshev sum bound might settle this.

#### Lemma attempt C (pairwise log-concavity decomposition)

For U_s(x) = (1+x)^s A(x) with A(x)=∑ a_j x^j, define
  U_k = ∑_j a_j C(s,k-j),
  r_j(k) = C(s,k+1-j)/C(s,k-j).
Then
  U_{k+1}^2 - U_k U_{k+2}
    = ∑_{i,j} a_i a_j C(s,k-i) C(s,k+1-j) [r_i(k) - r_j(k+1)].

Diagonal terms (i=j) are nonnegative because each binomial row is log-concave.
Thus a sufficient condition for tail log-concavity at index k is:

  r_i(k) >= r_j(k+1)  for all i,j with a_i,a_j>0.

Since r_j(k) increases in j and decreases in k, the strongest needed check is
  r_{j_{\max}}(k) >= r_{j_{\min}}(k+1).

This is not always true, but in right-tail regimes (k >= s/2 + d + O(1)) it may
hold uniformly, which would imply tail log-concavity of U_s without any LC
assumption on A. Turning this into a clean inequality is open, but it provides
an explicit target for a TP2/MLR argument.

Quick check: the “min >= max” condition reduces to r_0(k) >= r_d(k+1), i.e.
  (s-k)/(k+1) >= (s-k-1+d)/(k+2-d).
This inequality simplifies to (s+1)(1-d) >= 0, so it only holds for d <= 1.
Thus the pairwise-min/max sufficient condition is too strong for general d >= 2
(even in far-right tails). This route is therefore a dead end unless the
effective support of A excludes one of the extremes.

#### Finite-core fallback (kernelization + check)

If leaf-attachment MBI bounds force each hub to be leaf-light (s(v) ≤ λ0), then
the remaining class is finite once subdivision preservation removes degree-2
vertices. In that finite class, the ratio-tail conditions can be checked directly.

Evidence: for (b0=6, λ0=5), the ratio-tail conditions and boundary checks hold for
all edges (1,387,276 checks; `results/finite_core_ratio_b6_l5.json`).

This suggests a viable proof strategy:
  (i) prove leaf-heavy hubs are boundary-safe (Lemma~ref in notes/leaf_attachment_mbi.md),
  (ii) prove subdivision preservation for the finite residual class by direct
       verification of the ratio-tail conditions.

#### Lemma attempt D (fixed-core asymptotic ratio monotonicity)

For fixed core (A,B) in a leaf-attachment model U_s=(1+x)^s A, the ratio expansion
from notes/general_core_asymptotics.md gives, in the central window
k = floor(s/2) + y with |y| = O(1),

  r_k = U_{s,k+1}/U_{s,k} = 1 - (4y + 2 - 4μ)/s + O(1/s^2).

Hence r_{k+1} = r_k - 4/s + O(1/s^2), so r_k is strictly decreasing in k across
the entire descent window once s is large (O_H(1) terms fixed). This yields the
tail ratio monotonicity needed in Lemma 2 for leaf-heavy cases; the remaining
leaf-light cases are finite and can be checked.

## Next step (proof attempt)

1) Prove log-concavity of P_u and Q_u/x for all rooted trees by induction.
2) Use ratio-monotonicity to control modes of P_u P_v vs Q_u Q_v.
3) Prove a first-difference dominance inequality for A(x) vs I(T).
4) Conclude subdivision preserves unimodality.

Alternative line:
  - Express A in terms of B = Q/x and F = P - B:
      A = (1+x) B_u B_v + B_u F_v + F_u B_v + F_u F_v.
    If B and F are log-concave with a mode ordering, a structured
    sum-of-log-concave lemma might apply.

If successful, any minimal counterexample has no degree-2 vertices, which
dramatically reduces the search for extremal structure.
