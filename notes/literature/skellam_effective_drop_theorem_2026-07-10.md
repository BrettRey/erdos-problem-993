# Skellam Effective-Drop Theorem

Date: 2026-07-10
Status: proved as a limiting-family theorem; not yet a proof for finite signed Poisson-binomial laws

## Statement

Let

\[
Z\sim\operatorname{Pois}(\lambda)-\operatorname{Pois}(\eta),\qquad
V=\lambda+\eta\ge 1,
\]

with mass sequence \(c_z=\Pr(Z=z)\). Let \(M\) be the largest mode and
\(D=M+1\). Strict Bessel Turán positivity makes the interior Skellam mass
sequence strictly log-concave, so this \(D\) is the first strict descent.
Whenever \(c_D>0\), define

\[
\Delta_{\mathrm{eff}}(D)
=1-\frac{c_{D+1}c_{D-1}}{c_D^2}.
\]

Then

\[
\boxed{V\Delta_{\mathrm{eff}}(D)>\frac14.}
\]

The qualifier \(c_D>0\) excludes only the terminal descent for the reflected
\(\operatorname{Pois}(1)\) endpoint. If a terminal descent is assigned
\(\Delta_{\mathrm{eff}}=1\), it also satisfies the quarter bound.

## Inputs

Write

\[
\mu=\lambda-\eta,\qquad
x=2\sqrt{\lambda\eta},\qquad
V^2=\mu^2+x^2.
\]

For \(\lambda,\eta>0\), the Skellam mass is

\[
c_z=e^{-V}\left(\frac{\lambda}{\eta}\right)^{z/2}I_{|z|}(x).
\]

Two external results are used.

1. J. N. Darroch's
   [mode theorem](https://doi.org/10.1214/aoms/1177703287) for
   Poisson-binomial laws: every mode lies at distance less than one from the
   mean (Theorem 4).
2. Baricz's strict Turán lower bound, valid for \(x>0\) and
   \(\nu\ge-\tfrac12\):

   \[
   I_\nu(x)^2-I_{\nu-1}(x)I_{\nu+1}(x)
   >
   \frac{\nu+\tfrac12}{\nu+1}
   \frac{I_\nu(x)^2}
        {\sqrt{x^2+(\nu+\tfrac12)^2}}.
   \]

The second result is Theorem 1, equation (3.11), in Á. Baricz,
[Bounds for Turánians of Modified Bessel
Functions](https://arxiv.org/abs/1202.4853), arXiv:1202.4853v2.

## Closed Darroch Localization for Skellam Modes

Every Skellam mode \(m\) satisfies

\[
|m-\mu|\le1.
\]

To see this, approximate \(\operatorname{Pois}(\lambda)\) and
\(\operatorname{Pois}(\eta)\) by binomials
\(\operatorname{Bin}(n,\lambda/n)\) and
\(\operatorname{Bin}(n,\eta/n)\). After shifting the difference by \(n\),
this is a Poisson-binomial law with mean \(n+\mu\). Darroch localizes each of
its modes to the open unit interval about that mean. Coefficient convergence
then gives the closed interval for every unique limiting mode.

If the Skellam law has two adjacent tied modes, an arbitrarily small
exponential tilt \(c_z\mapsto e^{\theta z}c_z\) selects either one. Such a tilt
is again Skellam after changing the two Poisson rates, and its mean tends to
\(\mu\) as \(\theta\to0\). Taking the two one-sided limits gives the same
closed localization for both tied modes.

Consequently, for \(D=M+1\),

\[
\mu\le D\le\mu+2,
\qquad\text{and hence}\qquad
|D|\le|\mu|+2.
\]

Only this last orientation-free inequality is needed below.

## Interior Proof

Assume \(\lambda,\eta>0\). Put

\[
\nu=|D|,\qquad a=\nu+\frac12,\qquad b=\nu+1.
\]

The geometric prefactor in the Skellam mass cancels from
\(c_{D+1}c_{D-1}/c_D^2\). The absolute-value cases at \(D=0\), at positive
\(D\), and at negative \(D\) all reduce to the order-\(\nu\) Bessel
Turánian. Baricz therefore gives

\[
\Delta_{\mathrm{eff}}(D)
>
\frac{a}{b\sqrt{x^2+a^2}}.
\]

It is enough to prove

\[
\frac{Va}{b\sqrt{x^2+a^2}}\ge\frac14.
\]

After squaring and using \(x^2=V^2-\mu^2\), this is equivalent to

\[
F_\nu(V,\mu)
=(16a^2-b^2)V^2+b^2\mu^2-b^2a^2\ge0.
\]

The coefficient \(16a^2-b^2\) is positive. For
\(\nu=0,1,2\), the unrestricted lower bounds \(V\ge1\) and \(\mu^2\ge0\)
give

\[
F_0\ge\frac{11}{4},\qquad
F_1\ge23,\qquad
F_2\ge\frac{139}{4}.
\]

For \(\nu\ge3\), localization gives

\[
V\ge|\mu|\ge\nu-2.
\]

The coefficients of both \(V^2\) and \(\mu^2\) in \(F_\nu\) are
positive, so

\[
\begin{aligned}
F_\nu
&=(16a^2-b^2)V^2+b^2\mu^2-b^2a^2\\
&\ge (16a^2-b^2)(\nu-2)^2+b^2(\nu-2)^2-b^2a^2\\
&=a^2\bigl(16(\nu-2)^2-b^2\bigr)\\
&=a^2(3\nu-9)(5\nu-7)\\
&\ge0.
\end{aligned}
\]

In every interior case the non-strict algebraic bound combines with Baricz's
strict inequality to give

\[
V\Delta_{\mathrm{eff}}(D)>\frac14.
\]

## Zero-Rate Endpoints

If \(\eta=0\), then \(Z\sim\operatorname{Pois}(V)\). With
\(m=\lfloor V\rfloor\), the largest mode is \(m\), including the integer tie,
so \(D=m+1\) and

\[
\Delta_{\mathrm{eff}}=\frac1{m+2},\qquad
V\Delta_{\mathrm{eff}}=\frac{V}{m+2}.
\]

At \(V=1\), this gives \(1/3\); for large \(V\) it tends to one and is
always above \(1/4\). On
\(V\in[m,m+1)\), its minimum is \(m/(m+2)\), equal to \(1/3\) at \(m=1\).

If \(\lambda=0\), then \(Z=-\operatorname{Pois}(V)\). For nonintegral \(V\)
with \(m=\lfloor V\rfloor\), the largest mode is \(-m\), \(D=1-m\), and

\[
\Delta_{\mathrm{eff}}=\frac1m,\qquad
V\Delta_{\mathrm{eff}}=\frac Vm>1.
\]

For integral \(V=m\ge2\), the largest mode is \(1-m\), so \(D=2-m\) and
\[
V\Delta_{\mathrm{eff}}=\frac m{m-1}>1.
\]
At \(V=1\), the first strict descent is terminal and \(c_D=0\).

## Computational Falsification Harness

`scripts/verify_skellam_effective_drop.py`:

- evaluates the exact Skellam/Bessel formulas at 80 decimal digits;
- scans both rate orientations, the one-sided boundaries, small and large
  variances, and imbalances close to zero and to the boundaries;
- detects every observed interior mode change and bisects its transition
  curve (the tie-selection jump immediately inside the reflected
  \(\operatorname{Pois}(1)\) endpoint is checked separately);
- probes both sides of each tie;
- checks the closed Darroch window, the Baricz lower bound, and the quarter
  bound independently.

Its JSON output is
`results/skellam_effective_drop_verify_2026-07-10.json`.

## What This Does and Does Not Establish

This proves the target quarter-scale reserve inequality for the complete
Skellam family and
therefore removes the most natural diffuse-limit obstruction to the signed
Poisson-binomial program.

This Skellam argument does **not by itself** prove the finite signed theorem.
That finite theorem was subsequently proved by a different cubic-curvature
argument in
`notes/literature/universal_poisson_binomial_effective_drop_2026-07-10.md`.
A direct finite-to-Skellam comparison remains false in sampled boundary
regimes; the successful proof instead exploits finite Poisson-binomial cubic
inequalities.
