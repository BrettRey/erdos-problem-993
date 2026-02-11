# Lamport-style transition for unimodality (proof obligation)

Goal: recast unimodality as an inductive safety property under the rooted
composition step used in the tree DP. This avoids any global “B_T is LC”
assumption (which is false).

## Rooted composition step

Let the current parent state be (P,Q), and a child subtree state be (U,V):

- P, Q: parent polynomials (exclude/include root)
- U = I_child = P_child + Q_child
- V = P_child

Then after attaching the child,

  P' = P (U + V),    Q' = Q U,
  I' = P' + Q' = I U + P V,   where I = P + Q.

## Safety property (unimodality)

Let ΔX_k := X_{k+1} - X_k. The sequence is unimodal iff:

  if ΔX_k < 0 at some k, then ΔX_j <= 0 for all j >= k.

Equivalently, once the sequence starts decreasing, it never increases again.

## Sufficient condition for safety preservation

Assume IU is unimodal, and let d = first descent index of IU (Δ(IU)_d < 0).
If for all k >= d,

  Δ(PV)_k <= -Δ(IU)_k,      (difference dominance)

then Δ(IU + PV)_k <= 0 for all k >= d, and I' is unimodal.

A stronger sufficient condition is:

  Δ(PV)_k <= 0 for all k >= d,

since Δ(IU)_k <= 0 for k >= d.

## Empirical evidence (n <= 10)

Checked all trees up to n = 10 (networkx enumeration) and all nodes, with
multiple random child orderings. For each composition step:

- IU was unimodal.
- PV was nonincreasing from d onward.
- The difference dominance inequality held: Δ(PV)_k <= -Δ(IU)_k.

This is only evidence; no proof yet.

## Next possible lemmata

1) **Mode ordering**: show first descent of PV occurs no later than the
   first descent of IU in rooted composition.
2) **Ratio interlacing** for IU and PV to ensure the difference dominance.
3) **Local injection**: map independent sets counted by PV in the tail to
   sets in IU to compare first differences.

If any of these is proved, the inductive safety argument would yield
unimodality for all trees by structural induction on the DP.
