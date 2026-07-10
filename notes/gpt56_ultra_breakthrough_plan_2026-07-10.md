# GPT-5.6 Ultra Breakthrough Plan

Date: 2026-07-10

## Outcome

All three planned stages closed, and the factor-by-factor bridge was then
superseded by a stronger result. Hillion--Johnson cubic curvature propagation,
combined with explicit mass windows and exact Sturm/Bernstein certificates,
proves the quarter-scale effective drop and raw reserve for every finite
Poisson-binomial law with `V>=1` at a supported first descent. The active
bottleneck is no longer signed stability; it is perturbation of the certified
product term `A=(1+x)^sQ` by the hub-included term `xR`.

## Recommendation

Concentrate this research push on the signed Poisson-binomial reserve, but split it into three auditable stages:

1. seal the July 4 audit findings with exact certificates and plateau-safe definitions;
2. prove the effective-drop bound for the full Skellam limiting family;
3. attack the first finite stability bridge: adding one reflected Bernoulli to a one-sided low-probability Poisson-binomial law.

This is the best bounded analytic route. The mode--mean/tie-fugacity route has greater global leverage, but it currently lacks a comparable algebraic handle. ECMS remains potentially global but needs a genuinely tree-realizable invariant, not another abstract coefficient-shape closure conjecture.

The current manuscript should remain unchanged during these stages. The signed-reserve work is private proof development, and the active audit ledger explicitly bars manuscript promotion before external audit.

## What Discovery Changed

### Exact side-bound disproof

The corrected conditional side-bound cannot carry the signed theorem by itself. There is a compact exact witness:

\[
X=\operatorname{Bernoulli}(1/4)+\operatorname{Bin}(6,1/2),
\qquad Y=\operatorname{Bin}(4,1/2).
\]

At its first descent, \(D=2\) and \(V=43/16\), while

\[
V\max(B_X,B_Y)=\frac{346967}{1481760}<\frac14.
\]

The actual reserve is large:

\[
V\Delta_{\rm eff}=\frac{19393}{24696},
\qquad
V(1-R_+)=\frac{559}{588}.
\]

This exact row is a cleaner theorem-level disproof than the floating-point \(\operatorname{Bin}(500,1/2)\) row. The latter is also exact-certifiable and belongs to an infinite symmetric family whose corrected side-bound decays to zero.

### The one-sided theorem survives audit

After deleting zero parameters and proving that the post-descent coefficient remains inside support, the existing chain is valid:

\[
p_i\le\frac12,\quad V\ge1
\quad\Longrightarrow\quad
\Delta_{\rm eff}\ge\frac1{4V}.
\]

The missing support step is small but should be explicit. If \(m\) is the number of positive parameters and \(D\) the first strict descent, then \(D\le m-1\), so all three coefficients used in \(\Delta_{\rm eff}\) are positive.

### Strict-descent continuity is false

Let \(X\sim\operatorname{Bin}(5,1/2)\) and \(Y\sim\operatorname{Bernoulli}(q)\). The strict first descent moves from \(4\) to \(3\) for every \(q>0\), and

\[
\Delta_{\mathrm{eff}}(X-Y)\longrightarrow\frac12
\quad\text{while}\quad
\Delta_{\mathrm{eff}}(X)=\frac35.
\]

Thus no uniform vanishing-error comparison at the strict first-descent index can underwrite the perturbative route.

The repair is to introduce the first weak descent

\[
D_{\le}=\min\{k:a_k\le a_{k-1}\}.
\]

The same Newton/localization proof gives the plateau-safe endpoint

\[
\Delta_{\mathrm{eff}}(D_{\le})\ge\frac1{4V}.
\]

For \(\operatorname{Bin}(5,1/2)\), this value is \(1/2\), exactly matching the vanishing reflected perturbation.

### Two status surfaces are stale

The current manuscript has the right dependency map, but `notes/one_private_status.md`, `notes/conjecture_A_analysis.md`, and `README.md` retain stronger claims that are no longer defensible:

- Transfer plus Conjecture A does not by itself prove the Case-B hub comparison; the manuscript correctly keeps a separate Case-B hub bound.
- PNP does not by itself prove unimodality.
- \(d_{\rm leaf}\le1\) log-concavity is false in the July Ramos--Sun stress corpus.
- Positive log-concavity alone does not imply mode \(\le\lceil\mu\rceil\).

These surfaces should receive small status corrections, without rewriting the submitted manuscript.

## Stage 1: Audit-Hardening Patch

### Changes

1. Add an exact rational regression for the small half-heavy witness above.
2. Replace the large-float caveat in the ledger with the compact exact disproof; optionally retain the symmetric \(M\)-family formula as an asymptotic explanation.
3. Add the support-domain lemma to the one-sided note.
4. Record the strict-descent discontinuity and the weak-descent repair.
5. Correct the stale PNP, Case-B, log-concavity, and mode--mean statements in the private status notes and README.

### Verification

- `python3 -m unittest test_all.py -v`
- exact `Fraction` checks for all displayed rational values
- `git diff --check`

### Stopping rule

Do not proceed if the exact regression disagrees with the current analyzer definitions or if an external audit finds a defect in the one-sided support repair.

## Stage 2: Skellam Effective-Drop Theorem

### Target

Let

\[
Z\sim\operatorname{Pois}(\lambda)-\operatorname{Pois}(\eta),
\qquad V=\lambda+\eta\ge1,
\]

and let \(D\) be the first strict descent of its mass sequence \(c_z\). Prove

\[
V\left(1-\frac{c_{D+1}c_{D-1}}{c_D^2}\right)>\frac14.
\]

### Proof skeleton already pressure-tested

Set \(\mu=\lambda-\eta\). For \(\lambda,\eta>0\), write

\[
x=2\sqrt{\lambda\eta},
\qquad
c_z=e^{-V}\left(\frac\lambda\eta\right)^{z/2}I_{|z|}(x).
\]

Strict Bessel Turan positivity makes the Skellam mass sequence strictly
log-concave, hence it has either one mode or two adjacent tied modes. Let \(M\)
be the largest mode; then \(D=M+1\). Approximation by shifted
Poisson-binomial laws, with an arbitrarily small Skellam exponential tilt at
a tie, gives the closed Darroch localization

\[
|M-\mu|\le1.
\]

Thus \(|D|\le|\mu|+2\). Baricz's lower Turan bound gives, with
\(\nu=|D|\), \(a=\nu+1/2\), and \(b=\nu+1\),

\[
\Delta_{\rm eff}
>
\frac{a}{b\sqrt{x^2+a^2}}.
\]

Since

\[
V^2=\mu^2+x^2,
\]

the desired \(V\Delta_{\mathrm{eff}}>1/4\) is equivalent to positivity of

\[
F=(16a^2-b^2)V^2+b^2\mu^2-b^2a^2.
\]

For \(\nu=0,1,2\), the unrestricted worst values \(V=1,\mu=0\) give

\[
F=11/4,\quad 23,\quad 139/4,
\]

respectively. For \(\nu\ge3\), use
\(V\ge|\mu|\ge\nu-2\) to obtain
\(F\ge a^2(3\nu-9)(5\nu-7)\ge0\); Baricz's bound is strict. The
zero-rate endpoints are elementary reflected-Poisson calculations, and the
positive orientation has minimum \(V\Delta_{\mathrm{eff}}=1/3\) at
\(V=1\).

A first 80-digit scan over 216 variance/imbalance rows found minimum \(1/3\), at the one-sided \(\operatorname{Pois}(1)\) endpoint, and no \(1/4\) failure.

### Required hardening

1. Verify the exact statement and parameter range of the cited Baricz inequality from the primary paper.
2. Write the Skellam mode-localization limit carefully, including uniqueness and the one-sided boundary.
3. Add a 100-bit verification script that scans mode-transition curves rather than only a rectangular parameter grid.
4. Seek an independent line-by-line audit before calling this theorem-level.

### Falsifiers

- a high-precision Skellam row with \(V\Delta_{\mathrm{eff}}\le1/4\);
- a failure of the Baricz bound at integer order \(D\);
- a Skellam mode outside the closed Darroch window;
- an unhandled two-mode or zero-rate boundary.

## Stage 3: One-Reflected-Bernoulli Bridge

### Target

For a one-sided low-probability PB law \(X\) with \(V_X\ge1\), and

\[
Y\sim\operatorname{Bernoulli}(q),
\qquad 0\le q\le\frac12,
\]

prove at the first strict descent of \(Z=X-Y\) that

\[
\bigl(V_X+q(1-q)\bigr)\Delta_{\mathrm{eff}}(Z)\ge\frac14.
\]

The exact coefficient transform is

\[
c_k=(1-q)a_k+qa_{k+1}=a_k(1-q+qr_k),
\]

\[
\frac{c_k}{c_{k-1}}
=r_{k-1}\frac{1-q+qr_k}{1-q+qr_{k-1}}.
\]

Split the proof according to whether the reflected perturbation selects the original strict descent or the preceding weak descent. This is the smallest plateau-aware finite stability problem and the right test of whether the Skellam endpoint can be lifted back to finite signed PB laws.

### Breaker pass

- exact rational binomial and grouped-PB grids;
- explicit fair-plateau and near-plateau families with \(q\downarrow0\);
- optimization of \(V\Delta_{\mathrm{eff}}\), recording descent-index jumps as first-class data;
- symbolic checks of the two descent-selection cases.

### Pivot rule

If the one-reflected-Bernoulli lemma fails, use the exact witness to redesign the signed theorem before attempting multiple reflected variables. Do not proceed to a vague total-variance perturbation claim.

## Deferred Route: Mode--Mean / Tie Fugacity

The globally highest-leverage conjecture remains

\[
\operatorname{mode} I(T)\le\lceil\mu(T)\rceil
\]

or the narrower tie condition \(\mu(\lambda_m)\ge m-1\) for \(d_{\rm leaf}\le1\) trees. A future search should optimize the exact tie margin directly. It should not optimize log-concavity defects or \(n/3-\mu\), because neither is the missing implication.

This route is deferred during Stages 1--3 so that the current push produces a bounded theorem or a precise counterexample rather than another broad proof sketch.

## Success Criterion for This Push

Minimum success:

- exact audit-hardening patch merged locally;
- externally auditable Skellam theorem note with primary-source grounding and high-precision verification;
- a proved or exactly refuted one-reflected-Bernoulli lemma.

None of these outcomes should be described as a proof of the signed PB theorem, the Case-B hub bound, PNP, or Erdős Problem 993.
