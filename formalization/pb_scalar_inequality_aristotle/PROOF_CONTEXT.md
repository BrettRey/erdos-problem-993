# Proof context for Proposition 3.1

Source: Section 4 of *Variance and local log-concavity of Poisson–binomial laws*
(B. Reynolds, 2026), reproduced verbatim below after a roadmap. Lean
definitions: `PBScalar/Statement.lean`.

## Definitions (Sections 2–3 of the paper)

For `0 < δ < 1/4`:

- `K = max{r ∈ ℤ_{≥1} : (r+1)δ < 1}`; in Lean, `(K+1)δ < 1 ≤ (K+2)δ`.
- `a = (1-2δ)/(1-δ)`.
- `R_r = a^r ∏_{j=1}^{r-1}(1 - jδ)` and `L_r = (1-δ)^{-r} ∏_{j=2}^{r+1}(1 - jδ)`, for `1 ≤ r ≤ K`.
- Weights `w_0 = 1`, `w_r = R_r`, `w_{-r} = L_r`.
- `A(δ) = (1/2) ∑_{i,j=-K}^{K} w_i w_j (i-j)^2`.
- Target: `A(δ) ≥ (3+δ)/(4δ^2)`.

## Roadmap

1. **Change of variable.** Put `H = 1/δ - 1`, so `H > 3`. Then:
   - `K = max{r ≥ 1 : r < H}`;
   - `L_r = b_r := ∏_{s=1}^{r}(1 - s/H)`, which is ≥ 0 for `r ≤ K < H`;
   - the target's right-hand side is `Q(H) = (3H+4)(H+1)/4`.
2. **Right weights dominate left.** `R_r ≥ L_r` for `1 ≤ r ≤ K`. Put `C_r = R_r/L_r`; then `C_1 = 1` and `C_{r+1}/C_r = 1 + 2r/((H+1)(H-r-1)) ≥ 1`.
3. **Monotonicity.** All weights are ≥ 0, so `A` is nondecreasing in each weight: `∂A/∂w_ℓ = ∑_j w_j (ℓ-j)^2 ≥ 0`. Hence `A ≥ A_sym`, the value with `w_{±r} = b_r`.
4. **Pairwise identity.** `(1/2)∑∑ w_i w_j (i-j)^2 = (∑ w)(∑ i^2 w) - (∑ i w)^2`. In the symmetric case `A_sym = S·T`, where `S = 1 + 2∑_{r≤K} b_r` and `T = 2∑_{r≤K} r^2 b_r`.
5. **Compact range, `3 < H ≤ 16`** (i.e. `1/17 ≤ δ < 1/4`). There are 13 cells.
   - **Cell `[3,4]` (`K = 3`).** Keep the unsymmetrized weights. The identity `A - Q = P(H)/(4H^5(H+1)^3)` holds with the displayed degree-10 polynomial `P`. The degree-10 Bernstein coefficients of `P(3+t)` on `[0,1]` are all positive.
   - **Cells `[m, m+1]`, `m = 4, …, 15`.**
     - `P_m(H) = 4H^{2m}(S_m T_m - Q(H))` is an integer polynomial of degree `2m+2`, where `S_m` and `T_m` sum to `m`.
     - On `(m, m+1]` we have `K = m`. At `H = m` we have `K = m-1` and `b_m = 0`, so `S_m = S` and `T_m = T` on the closed cell.
     - The degree-`(2m+2)` Bernstein coefficients of `P_m(m+t)` are all positive.
   - **Data.** `data/universal_pb_finite_bernstein_full_certificate_2026-07-16.json` has a `payload.cells` list with 13 entries. Each entry gives `left_endpoint`, `degree`, `expected_denominator`, `numerator_power_coefficients_low_to_high` (the cleared numerator, in powers of `t` where `H = left_endpoint + t`) and `bernstein_coefficients_low_to_high`. Rationals are written as decimal `p` or `p/q`.
   - **Provenance.** Two independent programs generated and rechecked all 275 coefficients, and an independent third-party reconstruction (ChatGPT Pro, 2026-10-03) agreed exactly.
   - Since `binom(d,i) t^i (1-t)^{d-i} ≥ 0` on `[0,1]`, positivity of the coefficients gives positivity of the polynomial.
6. **Large range, `H ≥ 16`** (i.e. `δ ≤ 1/17`).
   - **Covering.** Choose `J ≥ 5` with `J(J+1)/2 ≤ H ≤ (J+1)(J+2)/2`. Then `J ≤ K`.
   - **Lower bounds.** `b_r ≥ λ_r := 1 - r(r+1)/(2H) ≥ 0` for `r ≤ J`. This is the Weierstrass product inequality `∏(1-x_s) ≥ 1 - ∑ x_s` for `x_s ∈ [0,1]`.
   - **Reduction.** The closed forms `S̃_J`, `T̃_J` give `S ≥ S̃_J ≥ 0` and `T ≥ T̃_J ≥ 0`, hence `ST ≥ S̃_J T̃_J`. It remains to show `N_J(H) := H^2(S̃_J T̃_J - Q(H)) > 0`, where `N_J` is the explicit quartic displayed below.
   - **`J = 5`.** Write `H = 16 + 5t`. The degree-4 Bernstein coefficients are `2360, 7500, 25055/2, 254205/16, 31115/2`.
   - **`J ≥ 6`.** Write `H = J(J+1)/2 + (J+1)t` and `J = u+6`. Then the coefficients are `β_i = μ_i π_i(u)/2880` with the printed `μ_i` and `π_i`, where every `π_i` has positive coefficients. These are polynomial identities in `u` (and `t`).
   - **Check.** `data/verify_pb_large_h_range.py` checks all of these exactly in SymPy.
7. **Assembly.** G1 and G2 cover `1/17 ≤ δ < 1/4` and `δ ≤ 1/17`, which gives G0.

## Section 4 of the paper (LaTeX, verbatim)

```latex
\section{The scalar inequality}
\label{sec:scalar}

Set
\begin{equation}\label{eq:H-def}
  H=\frac{1-\delta}{\delta}=\frac1\delta-1>3.
\end{equation}
Section~\ref{sec:compact} treats $3<H\leq16$ cell by cell, and
Section~\ref{sec:large} treats $H\geq16$. Then $K=\max\{r\in\mathbb Z_{\geq1}:r<H\}$, and the bounds $L_r$ of
\eqref{eq:left-mass-bound} become
\begin{equation}\label{eq:b-def}
  L_r=b_r:=\prod_{s=1}^r\left(1-\frac{s}{H}\right),
\end{equation}
which are nonnegative because $r\leq K<H$. For $1\leq r\leq K$, the bounds
$R_r$ of \eqref{eq:right-mass-bound} satisfy $R_r\geq b_r$. Indeed, if
$C_r=R_r/L_r$, then $C_1=1$ and, for $1\leq r<K$ (so that $H-r-1>0$),
\begin{equation}\label{eq:C-ratio}
  \frac{C_{r+1}}{C_r}
  =\frac{(H-1)(H+1-r)}{(H+1)(H-r-1)}
  =1+\frac{2r}{(H+1)(H-r-1)}\geq1.
\end{equation}

Since every weight is nonnegative, the form \eqref{eq:A-def} is
nondecreasing in each weight: with $K$ fixed,
\[
  \frac{\partial A}{\partial w_\ell}
  =\sum_{j=-K}^K w_j(\ell-j)^2\geq0.
\]
Replacing every $R_r$ by $b_r$ therefore does not increase $A$; denote the
resulting value by $A_{\mathrm{sym}}$. For the symmetric weights $w_0=1$ and
$w_{\pm r}=b_r$, put
\begin{equation}\label{eq:S-T}
  S=1+2\sum_{r=1}^K b_r,
  \qquad
  T=2\sum_{r=1}^K r^2b_r.
\end{equation}
The pairwise identity gives
\[
  A_{\mathrm{sym}}
  =\left(\sum_{r=-K}^K w_r\right)
   \left(\sum_{r=-K}^K r^2w_r\right)
   -\left(\sum_{r=-K}^K rw_r\right)^2
  =ST,
\]
because the final sum vanishes by symmetry. Hence
$A(\delta)\geq A_{\mathrm{sym}}=ST$. Under \eqref{eq:H-def}, the right-hand
side of \eqref{eq:scalar-target} is
\begin{equation}\label{eq:Q-def}
  Q(H)=\frac{(3H+4)(H+1)}4.
\end{equation}
For $H\geq4$, it remains to prove $ST\geq Q(H)$.

\subsection{The range \texorpdfstring{$3<H\leq16$}{3 < H <= 16}}
\label{sec:compact}

On $3<H\leq4$ we have $K=3$, and we keep the unsymmetrized bounds $L_r$ and
$R_r$, $r=1,2,3$. A SymPy computation (in the supplement) gives
\begin{equation}\label{eq:asymmetric-P}
  A(\delta)-Q(H)=\frac{P(H)}{4H^5(H+1)^3},
\end{equation}
where
\begin{align*}
  P(H)={}&-3H^{10}-16H^9+750H^8-3676H^7+6613H^6
  -5460H^5\\
  &+800H^4+4696H^3-7176H^2+4272H-912.
\end{align*}
We substitute $H=3+t$ and expand $P(3+t)$ in the degree-ten Bernstein basis
on $0\leq t\leq1$.

We call each interval $[m,m+1]$ a cell. For $m=4,\ldots,15$ and $H$ in the
cell $[m,m+1]$, define
\begin{equation}\label{eq:compact-ST}
  S_m=1+2\sum_{r=1}^m b_r,
  \qquad
  T_m=2\sum_{r=1}^m r^2b_r,
\end{equation}
and
\begin{equation}\label{eq:compact-P}
  P_m(H)=4H^{2m}\bigl(S_mT_m-Q(H)\bigr).
\end{equation}
For $m<H\leq m+1$, we have $K=m$. At $H=m$, we have $K=m-1$ and $b_m=0$.
Hence $S_m=S$ and $T_m=T$ on the whole cell. Each $P_m$ is an integer
polynomial of degree $2m+2$, and we expand $P_m(m+t)$ in the
degree-$(2m+2)$ Bernstein basis.

If $G(t)=\sum_{j=0}^d a_jt^j$ has degree $d$, the Bernstein conversion used
here is
\begin{equation}\label{eq:bernstein-conversion}
  \beta_i=\sum_{j=0}^i a_j
  \frac{\binom{i}{j}}{\binom{d}{j}},
  \qquad
  G(t)=\sum_{i=0}^d \beta_i\binom di t^i(1-t)^{d-i}.
\end{equation}
The expansions of $P(3+t)$ and of the twelve polynomials $P_m(m+t)$,
$4\leq m\leq15$, have $11$ and $264$ coefficients respectively, and all
$275$ are positive rationals. Since the Bernstein basis polynomials are
nonnegative on $[0,1]$, this proves \eqref{eq:scalar-target} for
$3<H\leq16$. A SymPy program computes each numerator and its Bernstein
coefficients. An independently written checker, using only exact rational
arithmetic from the Python standard library, recomputes them from the
defining formulas. Both programs and the
exact coefficients are in the supplement.

\subsection{The range \texorpdfstring{$H\geq16$}{H >= 16}}
\label{sec:large}

We cover $H\geq16$ by the intervals between consecutive triangular numbers:
choose an integer $J\geq5$ such that
\begin{equation}\label{eq:triangular-cell}
  \frac{J(J+1)}2\leq H\leq\frac{(J+1)(J+2)}2.
\end{equation}
Either choice of $J$ is valid at a shared endpoint. The first interval
($J=5$) is checked directly, and the rest ($J\geq6$) uniformly in $J$. If
$J\geq5$, then $J<J(J+1)/2\leq H$, so $J\leq K$. For $1\leq r\leq J$, the elementary product
bound gives
\begin{equation}\label{eq:bonferroni}
  b_r\geq\lambda_r:=1-\frac{r(r+1)}{2H}\geq0.
\end{equation}
Put
\begin{align}
  \tilde S_J&=1+2\sum_{r=1}^J \lambda_r
  =2J+1-\frac{J(J+1)(J+2)}{3H},
  \label{eq:S0}\\
  \tilde T_J&=2\sum_{r=1}^Jr^2\lambda_r
  =\frac{J(J+1)(2J+1)}3-\frac{\sigma_3+\sigma_4}{H},
  \label{eq:T0}
\end{align}
where
\[
  \sigma_3=\sum_{r=1}^Jr^3=\frac{J^2(J+1)^2}{4},
  \qquad
  \sigma_4=\sum_{r=1}^Jr^4=\frac{J(J+1)(2J+1)(3J^2+3J-1)}{30}.
\]
Since $J\leq K$, all omitted weights are nonnegative and $b_r\geq\lambda_r$
for $r\leq J$. Hence $S\geq\tilde S_J$, $T\geq\tilde T_J$, and
$A(\delta)\geq ST\geq\tilde S_J\tilde T_J$. It remains to prove
$\tilde S_J\tilde T_J\geq Q(H)$.

Let $N_J(H)=H^2\bigl(\tilde S_J\tilde T_J-Q(H)\bigr)$. Expanding
\eqref{eq:S0}, \eqref{eq:T0}, and \eqref{eq:Q-def} gives
\begin{align*}
  N_J(H)={}&-\frac34H^4-\frac74H^3
  +\left(\frac{J(J+1)(2J+1)^2}{3}-1\right)H^2\\
  &-\left((2J+1)(\sigma_3+\sigma_4)
    +\frac{J^2(J+1)^2(J+2)(2J+1)}{9}\right)H\\
  &+\frac{J(J+1)(J+2)}{3}(\sigma_3+\sigma_4).
\end{align*}
For $J=5$, only $16\leq H\leq21$ is needed. With $H=16+5t$ and
$0\leq t\leq1$, the degree-four Bernstein coefficients of $N_5(16+5t)$ in $t$
are
\[
  2360,\qquad 7500,\qquad \frac{25055}{2},\qquad
  \frac{254205}{16},\qquad \frac{31115}{2}.
\]

For $J\geq6$, write
\[
  H=\frac{J(J+1)}2+(J+1)t,
  \qquad J=u+6,
  \qquad 0\leq t\leq1.
\]
The degree-four Bernstein coefficients in $t$ of
$N_J\bigl(J(J+1)/2+(J+1)t\bigr)$ are $\beta_i=\mu_i\pi_i(u)/2880$ for
$i=0,\ldots,4$, where
\[
  \begin{gathered}
    \mu_0=J^2(J+1)^2,\qquad \mu_1=J(J+1)^2,\qquad \mu_2=(J+1)^2,\\
    \mu_3=(J+1)^2(J+2),\qquad \mu_4=(J+1)^2(J+2)^2,
  \end{gathered}
\]
and
\begin{align*}
  \pi_0(u)&=121u^4+2474u^3+17431u^2+46014u+25400,\\
  \pi_1(u)&=121u^5+3442u^4+37967u^3+199405u^2+480475u+384510,\\
  \pi_2(u)&=121u^6+4410u^5+65863u^4+512764u^3
    +2172926u^2+4668316u+3829440,\\
  \pi_3(u)&=121u^5+3684u^4+43735u^3+249961u^2+671859u+644400,\\
  \pi_4(u)&=121u^4+2958u^3+25579u^2+88782u+91440.
\end{align*}
Each $\mu_i$ is positive, and each $\pi_i$ has positive coefficients. Thus
$N_J(H)>0$ and $\tilde S_J\tilde T_J>Q(H)$ for $H\geq16$. A short SymPy
program in the supplement checks every expansion in this subsection in exact
arithmetic. This settles $H\geq16$; with Section~\ref{sec:compact}, it
proves Proposition~\ref{prop:scalar}.

```
