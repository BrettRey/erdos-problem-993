# General Core + Leaves Asymptotic for the Near-Miss Ratio

## Setup

Fix a tree \(H\) and a distinguished vertex \(v \in V(H)\).
For each \(s \ge 0\), let \(H_s\) be the tree obtained by attaching \(s\) new leaves to \(v\).

Write
\[
A(x) := I(H-v; x), \qquad B(x) := I(H-N[v]; x),
\]
so \(A(x)=\sum_{j=0}^d a_j x^j\) with \(a_j \ge 0\), \(A(1)>0\), and \(d=\deg A\) fixed (depends only on \(H,v\)).

Define
\[
\mu := \frac{A'(1)}{A(1)},\qquad
m := \left\lceil \mu + \frac12 \right\rceil,\qquad
C := (4m+2)-4\mu.
\]

Let \(c_k(s):=[x^k]I(H_s;x)\), and define the near-miss ratio
\[
\nm(H_s):=\max_{j\ge j_0}\frac{c_{j+1}(s)}{c_j(s)},
\]
where \(j_0\) is the first descent index (\(c_{j_0}(s)<c_{j_0-1}(s)\)); if no descent occurs, set \(\nm(H_s)=0\).

## Theorem (fixed core, large leaf attachment)

For fixed \((H,v)\), as \(s\to\infty\),
\[
\nm(H_s)=1-\frac{C}{s}+O_{H,v}\!\left(\frac1{s^2}\right),
\qquad
C=(4m+2)-4\mu,
\]
with \(\mu=A'(1)/A(1)\) and \(m=\lceil \mu+\tfrac12\rceil\).

Moreover,
\[
C\in[4,8),
\]
so \(\nm(H_s)<1\) for all sufficiently large \(s\). In particular, \(H_s\) is eventually unimodal.

## Proof

### 1) Closed form and decomposition (fixed-core only)

By include/exclude at the distinguished vertex \(v\),
\[
I(H_s;x)=(1+x)^sA(x)+xB(x).
\]
So with
\[
d_k(s):=[x^k](1+x)^sA(x),\qquad e_k:=[x^k]xB(x),
\]
we have \(c_k(s)=d_k(s)+e_k\), where \(e_k\) is independent of \(s\) and supported on finitely many \(k\) (depends only on \(H,v\)).

### 2) Central-window dominance of the main term (fixed-core only)

For \(k=\lfloor s/2\rfloor+O_{H,v}(1)\), \(d_k(s)\asymp_{H,v}\binom{s}{\lfloor s/2\rfloor}\), while \(e_k=O_{H,v}(1)\).
Hence
\[
\frac{e_k}{d_k(s)}=O_{H,v}(2^{-s}\sqrt{s}),
\]
so exponentially small in \(s\). Therefore, in the central window,
\[
\frac{c_{k+1}(s)}{c_k(s)}
=\frac{d_{k+1}(s)}{d_k(s)}+O_{H,v}(2^{-s}).
\]

Thus the near-miss asymptotics are governed by \(r_k(s):=d_{k+1}(s)/d_k(s)\).

### 3) Ratio expansion for the binomially-smoothed core (fixed-core only)

Write
\[
r_k(s)=
\frac{\sum_{j=0}^d a_j\binom{s}{k-j}\rho_j(k,s)}
     {\sum_{j=0}^d a_j\binom{s}{k-j}},
\qquad
\rho_j(k,s):=\frac{\binom{s}{k+1-j}}{\binom{s}{k-j}}
=\frac{s-k+j}{k+1-j}.
\]

Set \(k=s/2+y\) with \(|y|=O_{H,v}(1)\). For each fixed \(j\),
\[
\rho_j(k,s)=1-\frac{4y+2-4j}{s}+O_{H,v}\!\left(\frac1{s^2}\right),
\]
uniformly in bounded \(y\).

Also, for fixed \(j\), \(\binom{s}{k-j}/\binom{s}{k}=1+O_{H,v}(1/s)\), so the induced weights on \(j\) are
\[
\frac{a_j\binom{s}{k-j}}{\sum_t a_t\binom{s}{k-t}}
=\frac{a_j}{A(1)}+O_{H,v}\!\left(\frac1s\right),
\]
and their mean is \(\mu+O_{H,v}(1/s)\).
Therefore
\[
r_k(s)=1-\frac{4y+2-4\mu}{s}+O_{H,v}\!\left(\frac1{s^2}\right).
\]

### 4) First descent location and near-miss value (fixed-core only)

The linear term in \(y\) has slope \(-4/s\), so the ratio profile is strictly decreasing in \(k\) for large \(s\) (up to \(O(s^{-2})\) error).
Hence the first descent occurs when \(r_k(s)\) first drops below \(1\), i.e. near
\[
y>\mu-\frac12.
\]
Thus the first descent corresponds to \(y=m-1\), where \(m=\lceil \mu+\tfrac12\rceil\), and the near-miss ratio (which starts one step after first descent) is attained at \(y=m\):
\[
\nm(H_s)=1-\frac{(4m+2)-4\mu}{s}+O_{H,v}\!\left(\frac1{s^2}\right).
\]
This is exactly the claimed formula with \(C=(4m+2)-4\mu\).

### 5) Constant range

Let \(\delta:=m-(\mu+\tfrac12)\in[0,1)\). Then
\[
C=(4m+2)-4\mu=4+4\delta\in[4,8).
\]
So \(1-\nm(H_s)=C/s+O(s^{-2})>0\) for all sufficiently large \(s\), proving eventual unimodality.
\(\square\)

## What uses only “fixed core”, and what is path-specific?

### Uses only fixed-core structure

- The decomposition \(I(H_s;x)=(1+x)^sA(x)+xB(x)\).
- \(A,B\) are fixed polynomials independent of \(s\).
- The binomial smoothing ratio expansion in the central window.
- Exponential smallness of the \(xB(x)\) correction versus central binomial scale.
- The discrete threshold argument producing \(m=\lceil\mu+\tfrac12\rceil\).

### Path-specific ingredients

- None in the asymptotic argument above.

The broom case is recovered by taking \(H=P_p\) and \(v\) an endpoint:
\[
A(x)=I(P_{p-1};x),\quad B(x)=I(P_{p-2};x).
\]

## Connection to Galvin–Hilyard eventual unimodality under leaf attachment

Galvin–Hilyard-type eventual unimodality statements for repeated leaf attachment assert qualitative behavior: for fixed core, unimodality holds for all sufficiently large \(s\).

The theorem above is a quantitative refinement in the same regime:
- it gives a first-order asymptotic for the near-miss margin,
\[
1-\nm(H_s)=\frac{C}{s}+O\!\left(\frac1{s^2}\right),
\]
- and an explicit constant \(C\) from the core statistic \(\mu=A'(1)/A(1)\) via \(m=\lceil\mu+\tfrac12\rceil\).

So it explains not just that unimodality eventually holds, but how close the sequence is to violating unimodality and at what \(1/s\) rate the gap closes.

## Minimal edit plan for `paper/main.tex`

1. **Insert one theorem + short proof sketch** in Section `\section{Broom asymptotics}` (after the current asymptotic broom theorem):
   - State the general core theorem for \(H_s\).
   - Give the \(A,B,\mu,m,C\) definitions and the formula
     \(\nm(H_s)=1-C/s+O(1/s^2)\), \(C\in[4,8)\).

2. **Add a one-paragraph remark immediately after**:
   - “Brooms are the endpoint-core specialization.”
   - “No path-specific step is needed for the asymptotic mechanism.”

3. **Add one sentence in Discussion**:
   - Position this as a quantitative strengthening of eventual-unimodality-under-leaf-attachment results (Galvin–Hilyard line).

4. **(Optional) Appendix move**:
   - Keep the full technical proof in an appendix/notes-derived subsection and leave only theorem + proof sketch in the main flow.
