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

## Restricted-case lemma (broom root + leaf child)

Let $T_s=\broom(p,s)$ with $p \ge 2$ and $s \ge p$, rooted at the hub.
Write
  $P_s(x)=(1+x)^s A(x)$, $Q_s(x)=xB(x)$, $I_s=P_s+Q_s$,
where $A(x)=I(P_{p-1};x)$ and $B(x)=I(P_{p-2};x)$.
For Lamport composition with a leaf child $(U,V)=(1,x)$, we have
  $IU=I_s$ and $PV=xP_s$.
Let $d(F)=\min\{k:\Delta F_k<0\}$.

**Claim:** $d(PV) \ge d(IU)$.

**Sketch:** For $k \ge t+1$ where $t=\deg B=\lfloor (p-2)/2\rfloor$, we have
$\Delta I_{s,k}=\Delta P_{s,k}$ since $b_k=b_{k-1}=0$. Also
$\Delta P_{s,k}\ge 0$ for $k \le \lfloor s/2\rfloor-1$ by binomial monotonicity,
hence $d(P_s)\ge \lfloor s/2\rfloor\ge t+1$ when $s\ge p$.
Thus $d(I_s)\le d(P_s)$ and $d(PV)=d(xP_s)=d(P_s)+1\ge d(I_s)$.

This gives **mode ordering** for this specific Lamport step, but only in the
broom-root / leaf-child regime.
