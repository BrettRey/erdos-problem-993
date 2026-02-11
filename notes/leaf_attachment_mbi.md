# Leaf-Attachment MBI Lemmas (Boundary Control for I = (1+x)^s A + x B)

Setup: let
  I_s(x) = (1+x)^s A(x) + x B(x),
  A(x) = \sum_{j=0}^d a_j x^j,  a_j >= 0,  a_d > 0,
  B(x) = \sum_{j=0}^e b_j x^j,  b_j >= 0.

Define U_s(x) = (1+x)^s A(x), V(x) = x B(x), and
  \Delta F_k = [x^{k+1}]F - [x^k]F.

These lemmas give explicit boundary control of \Delta U_{s,k} near
k = t-2, t-1 (t = ceil((2 alpha_s - 1)/3), alpha_s = s + d).

---

## Lemma 1 (TP Sign-Strip / Binomial Increment)

For all k:

  \Delta U_{s,k} = \sum_{j=0}^d a_j [C(s, k+1-j) - C(s, k-j)].

Hence:

  k <= floor(s/2) - 1   =>  \Delta U_{s,k} >= 0,
  k >= ceil(s/2) + d    =>  \Delta U_{s,k} <= 0.

Reason: each binomial increment has known sign across the midpoint, and
all weights are nonnegative.

MBI link: once k is in the right strip and k > deg(V)=e+1, we have
\Delta V_k = 0 and thus \Delta V_k <= -\Delta U_{s,k} automatically.

---

## Lemma 2 (MLR Envelope for Binomial Mixtures)

Define
  U_{s,k} = [x^k] U_s = \sum_{j=0}^d a_j C(s, k-j).

Then for d <= k <= s-1,

  (s-k)/(k+1) <= U_{s,k+1}/U_{s,k} <= (s-k+d)/(k+1-d).

Proof: U_{s,k+1}/U_{s,k} is a convex combination of
  r_j(k) = C(s, k+1-j)/C(s, k-j) = (s-k+j)/(k+1-j),
and r_j increases in j, so it is bracketed by min/max.

Corollary:
  k >= (s + 2d - 1)/2  =>  U_{s,k+1} <= U_{s,k}, i.e. \Delta U_{s,k} <= 0.

---

## Lemma 2b (Variation-Diminishing / PF Kernel)

The binomial kernel (C(s,j)) is PF_\infty, so convolution is variation-diminishing:

  sc(\Delta U_s) <= sc(\Delta a),

where sc(.) counts sign changes and \Delta a_j = a_{j+1} - a_j.
If A is unimodal (sc(\Delta a) <= 1), then \Delta U_s has at most one sign change.
Since \Delta U_{s,0} = s + a_1 - 1 >= 0, once \Delta U_s turns negative it stays
negative.

MBI link: if the boundary index k_* lies past deg(V), then \Delta V_{k_*} = 0 and
the one-check condition \Delta U_{s,k_*} <= 0 suffices.

---

## Lemma 3 (Boundary Slack Lower Bound)

For k >= ceil(s/2) + d,

  -\Delta U_{s,k} >= a_0 (C(s,k) - C(s,k+1))
                  = a_0 * (2k+1-s)/(k+1) * C(s,k).

At k ~ 2s/3, the RHS is a fixed positive fraction of C(s,k), hence
exponentially large in s.

MBI link: for fixed-core V, |ΔV_k| is O(1), so for large s,
  \Delta V_k <= -\Delta U_{s,k}
at k = t-2, t-1.

---

## Corollary (Leaf-Heavy Cutoff; explicit s0)

Let e = deg B, d = deg A. A crude explicit sufficient condition for MBI at
k = t-2, t-1 is:

  s >= s0(d,e) := max(e+1-d, 2d+11, ceil((3e-2d+13)/2)).

Then:
  \Delta U_{t-2}, \Delta U_{t-1} < 0  and  \Delta V_{t-2} = \Delta V_{t-1} = 0.

Thus MBI holds automatically at those boundary indices.

Notes:
- The constants are sufficient, not optimized.
- The 2d+11 term ensures t_s-2 and t_s-1 lie beyond (s+2d-1)/2 so the MLR
  envelope gives \Delta U_{s,k} <= 0 at both indices.

Interpretation: any minimal boundary violator must have s bounded in terms of
the core degree (leaf-heavy cases are automatically safe).

---

## Probability view (random shift)

Let J have P(J=j)=a_j/A(1) and X ~ Bin(s,1/2). Then
  U_{s,k} / A(1) = P(X+J = k).
The MLR envelope above is the statement that the likelihood ratio
U_{s,k+1}/U_{s,k} is bounded by the min/max binomial ratios, which provides
one-line sign control of \Delta U_{s,k} at boundary indices.

---

## PF∞ One-Check Lemma (optional stronger hypothesis)

If A is PF∞ (equivalently real-rooted with nonpositive zeros), then the ratio
r_k = U_{s,k+1}/U_{s,k} is nonincreasing. Therefore, a single check
r_{k_*} <= 1 implies \Delta U_{s,k} <= 0 for all k >= k_*.

MBI link: with k_* = t-2 (or t-1) and k_* > deg(V), boundary control collapses
to one inequality.

---

## Consequence for minimal-counterexample strategy

The lemmas above show that large leaf attachment forces MBI at k=t-2,t-1.
Therefore, any minimal counterexample to unimodality must be leaf-light at
every high-degree vertex (bounded s relative to the core), reducing the problem
to a finite, small-core class.

This is the structural reduction step needed to make the MBI program viable.
