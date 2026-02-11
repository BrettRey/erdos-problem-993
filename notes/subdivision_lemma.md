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

Additional empirical facts (n <= 10 exhaustive; random larger tests):
  - The added term A(x) = Q_u Q_v + x P_u P_v is log-concave and unimodal.
  - The first descent index of I(T') is never earlier than that of I(T).
  - A(x) is nondecreasing up to the first descent of I(T).

These suggest that subdivision “delays” the peak rather than creating a valley,
but this is only evidence; not a proof.

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
