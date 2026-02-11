# Lamport-style transition for unimodality (proof obligation)

Goal: recast unimodality as an inductive safety property under the rooted
composition step used in the tree DP. This avoids any global “B_T is LC”
assumption (which is false).

## Rooted composition step

Let the current parent state be (P,Q), and a child subtree state be (U,V):

- P, Q: parent polynomials (exclude/include root)
- U = P_child (child excluded)
- V = Q_child (child included, with the x factor)

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
For Lamport composition with a leaf child, $(U,V)=(1,x)$ (since $P_{\text{leaf}}=1$,
$Q_{\text{leaf}}=x$), so $IU=I_s$ and $PV=xP_s$.
Let $d(F)=\min\{k:\Delta F_k<0\}$.

**Claim:** $d(PV) \ge d(IU)$.

**Sketch:** For $k \ge t+1$ where $t=\deg B=\lfloor (p-2)/2\rfloor$, we have
$\Delta I_{s,k}=\Delta P_{s,k}$ since $b_k=b_{k-1}=0$. Also
$\Delta P_{s,k}\ge 0$ for $k \le \lfloor s/2\rfloor-1$ by binomial monotonicity,
hence $d(P_s)\ge \lfloor s/2\rfloor\ge t+1$ when $s\ge p$.
Thus $d(I_s)\le d(P_s)$ and $d(PV)=d(xP_s)=d(P_s)+1\ge d(I_s)$.

This gives **mode ordering** for this specific Lamport step, but only in the
broom-root / leaf-child regime.

## Conjectural extension: broom root + path child

Let $T_s=\broom(p,s)$ rooted at the hub with $s \ge p$, and let the child
subtree be a path of length $\ell$ (i.e., $\ell$ edges, $\ell+1$ vertices)
rooted at its endpoint. Then the Lamport step uses
  $U = P_{\text{child}} = I(P_\ell;x)$ and $V = Q_{\text{child}} = x I(P_{\ell-1};x)$.

**Conjecture (mode ordering):** $d(PV) \ge d(IU)$.

**Stronger conjecture (difference dominance):**
for all $k \ge d(IU)$, $\Delta(PV)_k \le -\Delta(IU)_k$.

**Evidence:** exhaustive parameter scan for
  $2 \le p \le 6$, $p \le s \le 10$, $0 \le \ell \le 6$
found no violations of either inequality.

The proof attempt likely mirrors the leaf case:
  - express $P_s = (1+x)^s A$ and $I_s = P_s + x B$,
  - use binomial monotonicity for $P_s$ up to $\lfloor s/2\rfloor$,
  - show the $xB$ term is supported strictly below the first descent once
    $s$ is large relative to $\deg A$ and $\deg U$.

## Lemma (eventual mode ordering for broom + path child)

Fix $p \ge 2$ and path length $\ell \ge 0$. There exists $s_0(p,\ell)$ such that
for all $s \ge s_0(p,\ell)$, the Lamport step for the broom root with path child
satisfies $d(PV) \ge d(IU)$.

**Sketch of proof.**
Write $P_s=(1+x)^s A$, $Q_s=xB$ with $A=I(P_{p-1})$, $B=I(P_{p-2})$.
Let $C=I(P_\ell)$ and $D=I(P_{\ell-1})$, so $U=C$, $V=xD$.

1) The correction term $xBC$ has degree at most
   $1+\deg B + \deg C \le \lfloor (p+\ell)/2 \rfloor$.
   For $s \ge p+\ell+2$, we have $\lfloor s/2 \rfloor > \deg(xBC)$, so
   in the central window the coefficients of $IU=(P_s+xB)C$ agree with those
   of $P_s C$. Hence $d(IU)=d(P_s C)$ for all such $s$.

2) For any fixed polynomial $H$, the first descent of $(1+x)^s H$ occurs at
   $k = \lfloor s/2\rfloor + m_H - 1$ for large $s$, where
   $m_H=\lceil \mu_H + \tfrac12\rceil$ and $\mu_H=H'(1)/H(1)$ (by the
   same expansion used in the asymptotic theorem).
   Applying this to $H=A C$ and $H=A D$ gives
   \[
     d(P_s C) = \lfloor s/2\rfloor + m_{AC} - 1,\quad
     d(P_s D) = \lfloor s/2\rfloor + m_{AD} - 1
   \]
   for all sufficiently large $s$.

3) Since means add under convolution with nonnegative coefficients,
   $\mu_{AC}=\mu_A+\mu_C$ and $\mu_{AD}=\mu_A+\mu_D$.
   For paths, $\mu_n$ is strictly increasing and $\mu_n-\mu_{n-1}<1$:
   writing $F_n=I(P_n;1)$ and $G_n=I'(P_n;1)$, we have
     $F_n=F_{n-1}+F_{n-2}$ and
     $G_n=G_{n-1}+G_{n-2}+F_{n-2}$,
   hence
     $\mu_n=\frac{F_{n-1}\mu_{n-1}+F_{n-2}(1+\mu_{n-2})}{F_n}$.
   Since $\mu_{n-1}\ge\mu_{n-2}$, this gives
     $\mu_n-\mu_{n-1}\le F_{n-2}/F_n<1$.
   Therefore $m_{AC} \le m_{AD}+1$.

4) Therefore $d(P_s C) \le d(P_s D)+1$, and since $PV=xP_s D$ has
   $d(PV)=d(P_s D)+1$, we obtain $d(PV) \ge d(IU)$ for all large $s$.

This establishes mode ordering in the broom-root + path-child regime
for sufficiently large $s$.
