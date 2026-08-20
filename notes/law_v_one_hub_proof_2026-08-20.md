# LAW-V for a trivial hub: a complete proof

Date: 2026-08-20. Companion verifier:
`verify_law_v_one_hub_20260820.py`.

## Result

Use the notation of `notes/c2_single_index_laws_2026-08-11.md`.  For an
arm multiset \(a\), put

\[
 Q=\prod_i P_{a_i},\qquad R=\prod_i P_{a_i-1},
\]

where \(P_j=I(P_j;x)\), with \(P_{-1}=P_0=1\).  If the other hub is
trivial, then

\[
 C=Q+xR,\qquad B=Q+x(Q+R)=C+xQ.
\]

**Theorem.** For every finite multiset of positive arm lengths and every
index \(k\),

\[
 V_k=c_kb_k-c_{k+1}b_{k-1}\geq 0.
\]

Thus G3 of `gpt_attack/law_v_packet_2026-08-12/law_v_packet.py` is
proved.  G1 is its one-arm special case.  By symmetry, the result holds
whenever either hub is trivial.

## Lemma 1: extending every path by one vertex raises coefficient ratios

Write \(p_{n,j}=[x^j]P_n=\binom{n+1-j}{j}\).  Then

\[
 p_{n,j}p_{n-1,j-1}\geq p_{n,j-1}p_{n-1,j}
 \tag{1}
\]

at every index, with out-of-support coefficients interpreted as zero.
On the common positive support this follows because

\[
 \frac{p_{n,j}}{p_{n-1,j}}
 =\frac{n+1-j}{n+1-2j}
\]

is increasing in \(j\).  The support-boundary cases are immediate.
In likelihood-ratio notation, (1) says \(P_{n-1}\preceq_{lr}P_n\).

Likelihood-ratio order is preserved by common convolution with a
log-concave sequence having no internal zeros.  For completeness, place
the two sequences in a two-row array.  The consecutive minors encode
likelihood-ratio order.  Convolution on both rows multiplies this array
by the Toeplitz matrix of the common factor; that matrix is TP2 exactly
when the factor is log-concave.  Cauchy--Binet makes every resulting
minor a sum of products of nonnegative minors.  This is the same TP2
argument as Lemma K2 in
`notes/karlin_lr_order_toolkit_2026-08-12.md`, without the harmless row
shift.

Replace the factors \(P_{a_i-1}\) by \(P_{a_i}\), one at a time.  Every
other factor is a path polynomial and hence log-concave with no internal
zeros.  Common-convolution closure and transitivity give

\[
 R\preceq_{lr}Q,
 \qquad r_jq_{j-1}\leq r_{j-1}q_j.
 \tag{2}
\]

## Lemma 2: the shifted cross term is nonnegative

The path polynomials are real-rooted, so their product \(Q\) is
log-concave.  Ratio monotonicity for \(Q\), followed by (2), gives

\[
 q_kq_{k-1}\geq q_{k+1}q_{k-2},
 \tag{3}
\]

and

\[
 \frac{r_k}{r_{k-1}}
 \leq \frac{q_k}{q_{k-1}}
 \leq \frac{q_{k-1}}{q_{k-2}},
 \qquad
 r_{k-1}q_{k-1}\geq r_kq_{k-2}.
 \tag{4}
\]

Again the claims are immediate if a denominator would be outside the
positive support.  Since \(c_j=q_j+r_{j-1}\), (3)--(4) imply

\[
\begin{aligned}
c_kq_{k-1}-c_{k+1}q_{k-2}
 &=q_kq_{k-1}-q_{k+1}q_{k-2}\\
 &\quad+r_{k-1}q_{k-1}-r_kq_{k-2}\geq0.
\end{aligned}
\tag{5}
\]

## Proof of the theorem

The polynomial \(C=Q+xR\) is the independence polynomial of the spider
with arm multiset \(a\).  Li, Li, Yang, and Zhang, Theorem 3.1
(arXiv:2501.04245; text on disk at
`notes/literature/li_yang_zhang_zhang_2025_symmetric_function_log_concavity.txt`)
prove that every spider independence polynomial is log-concave.  Hence
\(c_k^2-c_{k-1}c_{k+1}\geq0\).  Using \(B=C+xQ\) and (5),

\[
\begin{aligned}
V_k
 &=c_k(c_k+q_{k-1})-c_{k+1}(c_{k-1}+q_{k-2})\\
 &=\bigl(c_k^2-c_{k-1}c_{k+1}\bigr)
   +\bigl(c_kq_{k-1}-c_{k+1}q_{k-2}\bigr)\geq0.
\end{aligned}
\]

This proves the theorem. \(\square\)

## What the proof localizes

The tempting two-hub extension is to write, for general nontrivial
hubs,

\[
 B=C+xH,\qquad
 H=R_aQ_b+Q_aR_b-R_aR_b,
\]

and hope that the second bracket in the proof stays nonnegative by the
same likelihood-ratio argument.  It does not: exact small examples such
as \(a=(1,1)\), \(b=(1,1,5)\) make
\(c_kh_{k-1}-c_{k+1}h_{k-2}<0\), even though the full \(V_k\) remains
positive.  Therefore the general LAW-V genuinely needs quantitative
absorption by the spider Turan gap.  It cannot follow from a direct
two-hub likelihood-ratio order.  The companion verifier includes this
counterexample as a regression.

## Verification status

The verifier checks the path identity and its likelihood-ratio minors,
the product order \(R\preceq_{lr}Q\), every displayed nonnegative
summand, and LAW-V itself for all arm multisets of at most four arms,
individual lengths at most 10, and total arm weight at most 24.  It also
checks the stated failure of the naive two-hub lift.  These computations
audit the proof; they are not used in place of its universal steps.
