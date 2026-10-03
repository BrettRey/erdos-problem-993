<!-- Pasted by Brett into the Claude Code session on 2026-10-03: ChatGPT Pro's cross-family referee check of the paper, run on the brief submission/external-check-brief-2026-10-03.md. It reviewed the second-revision PDF (pre-retitle), without the supplement. Saved verbatim by the parent session as text. The report's own attachment (turan-independent-checks-2026-10-03.zip, standard-library Python reconstructing all 275 coefficients) lives in ChatGPT's sandbox and is not in this repository. Its two Pitman corrections and three clarifications were verified (Pitman (20) read from the rendered PDF page 8; the balanced-family degeneration recomputed numerically) and applied; see DECISIONS 2026-10-03 (night). -->

## Referee verdict

**I would recommend acceptance after minor revision.** I found no mathematical error in Theorem 1.1, the new corollaries, the modal upper bound, or the reduction and certificate proof in Sections 2–4. I independently reconstructed **all 275 finite-range Bernstein coefficients**, the remaining polynomial identities, and the CUE threshold calculation. They passed.

There are **two corrections to make in the Pitman discussion**, several smaller clarifications, and the unresolved supplement link. The first Pitman correction also corrects part of my earlier advice: I was too categorical about what his bounds provide.

I checked the revised manuscript against your ten-point brief. The certificate supplement mentioned in the brief was not attached, so my arithmetic checks were reconstructed from the manuscript's defining formulas. I have not checked the original supplement's manifest, programs, or Lean files. variance-scaled-turan-first-descent(1).pdf external-check-brief-2026-10-03.md

The independent programs, full rational coefficients, and verification logs are here:

:chatgpt-content-reference{index="7"}[turan-independent-checks-2026-10-03.zip](sandbox:/workspace/scratch/960bd6a33166/output/turan-independent-checks-2026-10-03.zip)

Both programs use only the Python standard library.

## Changes I would require

### 1. Page 2: the claim about Pitman is too strong

**Classification: an error in the description of earlier work.**

The introduction says that Darroch's theorem, Pitman's ratio bounds, and maximal-mass bounds "by themselves give no lower bound on \(\delta_k\)." Pitman's sharper equation (20) does give a positive lower bound.

Write \(\theta(x)\) for the multiplicative tilt whose tilted mean is \(x\). Pitman gives
\[
\theta\!\left(k+\frac1{k+2}\right)
\le \frac{f_k}{f_{k+1}}
\le
\theta\!\left(k+1-\frac1{n-k+1}\right).
\]
His equation (21) is a coarser consequence. [stat.berkeley.edu](https://www.stat.berkeley.edu/~pitman/453.pdf?utm_source=chatgpt.com)

Combining equation (20) at consecutive indices, I obtain
\[
\delta_k
\ge
1-
\frac{
\theta\!\left(k-\frac1{n-k+2}\right)
}{
\theta\!\left(k+\frac1{k+2}\right)
}
>0,
\qquad 1\le k\le n-1.
\]
The strict positivity follows because \(\theta\) is strictly increasing.

**Fix:** narrow the claim to the actual gap: these results do not directly supply a positive lower bound on \(V\delta_D\) that is uniform over the class under consideration. Do not say that they give no deficit bound.

There is a useful strengthening available here. **Even the bound obtained directly from Pitman's sharper equation (20) degenerates in your existing balanced variance-one family.** Thus acknowledging it does not weaken the central motivation.

For that family, with \(m\) parameters equal to \(\varepsilon_m\), \(m\) equal to \(1-\varepsilon_m\), and
\[
2m\varepsilon_m(1-\varepsilon_m)=1,
\]
the centered tilted mean is exactly
\[
\mu_m(t)-m
=
\frac{\sinh t}{1+(\cosh t-1)/m}.
\]
This follows by differentiating
\[
\log \mathbb E e^{tW_m}
=
mt+m\log\!\left(1+\frac{\cosh t-1}{m}\right).
\]

At \(D=m+1\), the directly combined Pitman bound is
\[
1-
\frac{
\theta_m\!\left(m+1-\frac1{m+1}\right)
}{
\theta_m\!\left(m+1+\frac1{m+3}\right)
}.
\]
Both logarithmic tilts tend to \(\operatorname{arsinh}(1)\), so this expression tends to zero. More precisely, the expression is asymptotic to \(\sqrt2/m\).

That is a precise comparison worth using: **a specified consequence of the earlier theorem degenerates at fixed variance, while your theorem remains bounded below by \(1/4\).** It avoids making an unnecessarily broad claim about every possible use of Pitman's work.

### 2. Page 4: "stronger than (1.5)" needs a common quantity

**Classification: an incorrect comparison as written; the ratio inequality itself is correct.**

The statement
\[
\frac{f_{k+1}}{f_k}<\frac{V}{V+k-\mathbb EW}
\]
bounds an adjacent mass ratio. Equation (1.5) bounds a Turán deficit. A direct assertion that one is stronger than the other leaves the comparison unspecified.

**Fix:** either remove that sentence and describe the ratio bound's behavior, or compare it with a ratio consequence of (1.5).

For example, writing \(q_j=f_j/f_{j-1}\), equation (1.5) and \(q_D<1\) imply
\[
\frac{f_{k+1}}{f_k}
=q_{k+1}
<
\prod_{j=D}^{k}
\left(1-\frac1{4V+j-D}\right)
=
\frac{4V-1}{4V+k-D},
\qquad D\le k<n.
\]
Pitman's ratio estimate can then be compared with this expression. It is indeed smaller sufficiently far to the right.

I would use the shorter fix. The paragraph's useful point is already clear: the stated consequence of equation (21) approaches one when \(D-\mathbb EW\) approaches zero, whereas (1.6) gives a uniform variance-scaled gap.

### 3. Page 4: specify which bound fails for ULC laws

**Classification: a scope clarification.**

"No variance-scaled bound holds for ULC laws in general" is broader than Proposition 1.5 establishes.

Your example establishes the failure of a positive uniform lower bound on \(V\delta_D\), and consequently the failure of the corresponding support-wide conclusion. Its deficit at the mode behaves differently:
\[
\delta_m
=
1-\frac{m^2}{4(m+1)^2}
\longrightarrow \frac34.
\]
Thus this particular family does not establish failure of the modal lower bound.

**Fix:** specify a "variance-scaled lower bound at the first descent." That is exactly what the proposition proves and exactly what this part of the argument needs.

### 4. Section 2.1 and Example 1.6: small notation clarifications

**Classification: presentation issues.**

For differences of Bernoulli sums, define the first descent on the translated law's actual support. The literal restriction \(1\le k\le n\) in (1.2) should not be carried over to a law with possibly negative support. The complement-and-translation argument itself is correct.

In Example 1.6, I would describe \(\Gamma\) as the Gram matrix having the same nonzero eigenvalues as the restricted kernel operator. This makes the \(T^*T\)/\(TT^*\) relationship explicit and avoids a basis or Fourier-sign convention becoming a distraction. It does not change the eigenvalues or any calculation.

### 5. Supplementary Material: replace the placeholder

**Classification: a reproducibility requirement.**

`ZENODO-DOI-PENDING` must become a working location for the stated archive before the paper is finalized. My reconstruction verifies the arithmetic underlying the proof, but it does not verify that the intended release contains the promised files or that its commands work.

The description of the Lean component is appropriately limited: it describes conditional formalization of specified steps. I would preserve that precision.

## Mathematical checks against the brief

### 1. Reciprocal bound (2.5), including the endpoints — confirmed

The stated cubic inequalities agree with the Hillion–Johnson inequalities, and the divisions leading to
\[
\delta_{k-1}(1-\delta_k)\le\delta_k,
\qquad
\delta_{k+1}(1-\delta_k)\le\delta_k
\]
are legitimate for \(1\le k\le n-1\): all the mass factors being divided by are positive. [arxiv.org](https://arxiv.org/pdf/1303.3381?utm_source=chatgpt.com)

For an interior adjacent pair, applying the reciprocal inequalities in both directions gives
\[
\left|\frac1{\delta_{j+1}}-\frac1{\delta_j}\right|\le1.
\]

The endpoints also work. At the left endpoint, the recurrence with \(k=1\) and \(\delta_0=1\) gives
\[
1-\delta_1\le\delta_1,
\]
hence
\[
1\le\frac1{\delta_1}\le2.
\]
This proves the required bound against \(1/\delta_0=1\). The right endpoint is identical.

**There is no missing endpoint assumption or division by a zero boundary mass.**

### 2. Lemma 2.1 and its induction — confirmed

The induction is sound, including the step that excludes an endpoint before continuing.

The relevant sequence is:

1. The preceding index lies in the interior.
2. Therefore the next index lies in the closed support.
3. The reciprocal estimate gives
   \[
   \frac1{\delta_{D\pm r}}\ge\frac1\delta-r>1.
   \]
4. An endpoint has deficit exactly one, so the next index is actually interior.
5. Inversion yields
   \[
   \delta_{D\pm r}\le\frac{\delta}{1-r\delta}<1.
   \]

The strict cutoff \((r+1)\delta<1\) is doing the needed work. It also handles values of \(\delta\) whose reciprocal is an integer; there is no off-by-one error in the definition of \(K\).

### 3. Corollary 1.3 — confirmed

For the first inequality, the direction is correct:
\[
\frac1{\delta_k}
\le
\frac1{\delta_D}+|k-D|
\le
4V+|k-D|.
\]
Taking reciprocals gives the stated lower bound.

For the mean-centered version, Darroch gives
\[
|c-\mathbb EW|<1,
\qquad c=D-1,
\]
so
\[
|D-\mathbb EW|<2.
\]
The triangle inequality therefore produces the second denominator in (1.5), with the correct direction.

For (1.6), with
\[
\eta_-=\frac{f_D}{f_{D-1}}<1,
\qquad
\eta_+=\frac{f_{D+1}}{f_D}\ge0,
\]
we have
\[
1-\eta_+
\ge
1-\frac{\eta_+}{\eta_-}
=
\delta_D.
\]
The endpoint case \(D=n\) is also covered because \(f_{n+1}=0\).

The stated consequence on a window of width \(O(\sqrt V)\) is valid: throughout such a window the lower bound retains order \(1/V\).

### 4. Proposition 1.4: modal upper bound — confirmed

I checked all four points singled out in the brief.

**Tilted variance.** For a single summand,
\[
v_i(t)
=
\frac{p_i(1-p_i)e^t}{(1-p_i+p_ie^t)^2}.
\]
When \(t\ge0\), the denominator before squaring is at most \(e^t\), giving
\[
v_i(t)\ge e^{-t}v_i(0).
\]
When \(t\le0\), it is at most one, giving
\[
v_i(t)\ge e^t v_i(0).
\]
Summing proves \(V(t)\ge e^{-|t|}V\).

**Exactly two modes.** Exponential tilting preserves strict log-concavity and multiplies each adjacent ratio by \(e^t\). At \(t_-\), the ratio between \(c-1\) and \(c\) becomes one; all earlier ratios are larger and all later ratios smaller. Thus those are exactly the two modes. The argument at \(t_+\) is the same.

This remains valid when \(t_-=0\), meaning the original distribution already has two modes.

**Strict Darroch step.** Applying the strict distance bound to both modes gives
\[
\mu(t_-)\in(c-1,c),
\qquad
\mu(t_+)\in(c,c+1).
\]
The intersection of the two restrictions is important, and you have used it correctly.

**Integral.** Consequently,
\[
2>
\mu(t_+)-\mu(t_-)
=
\int_{t_-}^{t_+}V(t)\,dt
\ge
V\left(2-\frac1{q_c}-q_{c+1}\right).
\]
Finally,
\[
\left(2-\frac1{q_c}-q_{c+1}\right)
-
\left(1-\frac{q_{c+1}}{q_c}\right)
=
\left(1-\frac1{q_c}\right)(1-q_{c+1})\ge0.
\]

Hence \(\delta_c<2/V\), with the strict inequality justified. Combining this with the lower bound gives the stated
\[
\frac15\le V\delta_c<2.
\]

### 5. Proposition 1.5: ULC counterexample — confirmed

The normalized sequence
\[
\frac{g_{m+j}}{\binom{2m}{m+j}}
\propto 2^{-|j|}
\]
is log-concave. Multiplication by the strictly log-concave binomial row makes \(g\) strictly log-concave.

Symmetry then gives the unique mode \(m\) and first descent \(m+1\).

The cancellation of the geometric factors at the first descent is correct:
\[
\delta_{m+1}
=
1-
\frac{\binom{2m}{m}\binom{2m}{m+2}}
{\binom{2m}{m+1}^2}
=
\frac{2m+1}{m(m+2)}.
\]

For the variance limit, the centered binomial ratios are bounded by one and tend to one for each fixed integer \(j\). The stated dominating envelope is summable. The limiting centered distribution has weights proportional to \(2^{-|j|}\), with
\[
\sum_{j\in\mathbb Z}2^{-|j|}=3,
\qquad
\sum_{j\in\mathbb Z}j^2\,2^{-|j|}=12.
\]
Its variance is therefore \(12/3=4\).

The observation immediately after the proof is also useful and correct: the reciprocal deficit can jump from a bounded value at the mode to an arbitrarily large value one step later. That directly exhibits behavior excluded by the Bernoulli-sum recurrence.

### 6. Example 1.6: CUE half-circle count — confirmed

**Bernoulli representation.** The cited determinantal-process theorem gives the distribution of the count as a sum of independent Bernoulli variables with parameters given by the eigenvalues of the restricted operator. The cited Meckes–Meckes proposition applies this representation to the unitary eigenvalue count. [arxiv.org](https://arxiv.org/pdf/math/0503110?utm_source=chatgpt.com)

**Parameters in \((0,1)\).** Both \(\Gamma\) and \(I-\Gamma\) are positive definite: a nonzero trigonometric polynomial cannot vanish on an interval. Thus all \(N\) parameters are strictly between zero and one.

**Exact variance.** The diagonal entries are \(1/2\); off the diagonal,
\[
|\Gamma_{jk}|^2=
\begin{cases}
1/(\pi^2(j-k)^2),&j-k\text{ odd},\\
0,&j-k\text{ even}.
\end{cases}
\]
Counting the two matrix diagonals at distance \(d\) gives
\[
V_N
=
\frac N4-\frac2{\pi^2}
\sum_{\substack{1\le d<N\\d\text{ odd}}}\frac{N-d}{d^2}.
\]
The conversion to
\[
V_N=\frac2{\pi^2}
\sum_{\substack{d\ge1\\d\text{ odd}}}
\frac{\min(d,N)}{d^2}
\]
is correct. This representation also proves strict monotonicity in \(N\).

**Asymptotic.** For even \(N\), my expansion gives the slightly stronger statement
\[
V_N=
\frac{\log N+\gamma+\log2+1}{\pi^2}
-\frac1{6\pi^2N^2}
+O(N^{-4}).
\]
Your \(O(N^{-2})\) remainder is therefore correct. In particular, there is no missing \(1/N\) term.

**Threshold.** I checked the two neighboring even values using exact rational sums and rigorous rational bounds on \(\pi\), rather than ordinary floating-point evaluation:

| Quantity | Value, rounded for display |
|---|---:|
| \(V_{1996}\) | \(0.99996543523160400122908765\) |
| \(V_{1998}\) | \(1.00006690864230001837548359\) |
| \(V_{10000}\) | \(1.16323843886831071519829809\) |
| \(1/(4V_{10000})\) | \(0.21491724451886196499928543\) |
| Newton bound at \(N=10000\) | \(0.00039988004798080767692923\) |

The enclosed intervals certify
\[
V_{1996}<1<V_{1998}.
\]
Together with strict monotonicity, this proves that **1998 is the first even threshold**, as claimed.

**First descent.** Rotation by \(\pi\) exchanges the two half circles. Eigenvalues on their boundary have probability zero, so \(X_N\) and \(N-X_N\) have the same law. Symmetry and strict log-concavity give the unique mode \(N/2\), hence \(D=N/2+1\).

### 7. The displayed Pitman ratio estimate — confirmed

The derivation on page 4 is correct.

For \(t\ge0\),
\[
V(t)\le e^tV,
\]
because each denominator \(1-p_i+p_ie^t\) is at least one. If \(\mu(t_k)=k>\mu(0)\), then
\[
k-\mu(0)
=
\int_0^{t_k}V(t)\,dt
\le V(e^{t_k}-1).
\]
Combining this with the strict ratio inequality supplied by equation (21) gives the displayed bound.

Your binomial example also has the stated mode and first descent when \(\varepsilon\) is sufficiently small and \(p<\varepsilon\). The resulting upper bound tends to one while the variance tends to \(\ell+1\).

The needed corrections concern the interpretation and comparison of the estimates, as detailed above.

### 8. Sections 2–4: main proof and certificates — confirmed

#### Standing reductions

Deleting deterministic-zero summands and translating away deterministic-one summands preserve the variance and corresponding deficits.

The argument proving the existence of \(D\) is correct:
\[
\frac{f_n}{f_{n-1}}
=
\left(\sum_i\frac{1-p_i}{p_i}\right)^{-1}.
\]
If this ratio were at least one, the displayed comparison with \(V\) would contradict \(V\ge1\).

Likewise, \(D=1\) would make zero a mode, force \(\mathbb EW<1\), and hence force \(V<1\). Thus \(2\le D\le n\).

The extension to independent differences follows correctly by complementing the subtracted Bernoulli variables and translating.

#### Mass bounds

I checked the indexing and telescoping products in (2.9) and (2.10). They yield
\[
f_{c+r}\ge MR_r,\qquad f_{c-r}\ge ML_r.
\]
Lemma 2.1 guarantees that all the masses used actually lie in the support.

The asymmetry of the two products is necessary and has been handled correctly. In particular, the left-hand argument uses \(q_D<1\) in the right direction.

#### Variance reduction

The pairwise identity gives
\[
V\ge M^2A(\delta).
\]
The supplied proof of
\[
M^2\ge\frac1{1+12V}
\]
is valid. The comparison with the uniform density of height \(M\) has the correct sign both inside and outside its support.

One can see the final implication directly by combining these estimates:
\[
V(1+12V)\ge A(\delta)
\ge\frac{3+\delta}{4\delta^2}.
\]
Since \(v\mapsto v(1+12v)\) is strictly increasing for \(v\ge0\), and
\[
\frac1{4\delta}\left(1+\frac3\delta\right)
=
\frac{3+\delta}{4\delta^2},
\]
this forces \(V\ge1/(4\delta)\). Your contradiction proof expresses the same valid argument.

#### Symmetrization and range coverage

The calculation \(R_r\ge L_r=b_r\) is correct. The quadratic form \(A\) is coordinatewise nondecreasing on nonnegative weights, so replacing \(R_r\) by \(L_r\) indeed gives a lower bound. The symmetric form equals \(ST\).

I checked the range boundaries:

- On \(3<H\le4\), the unsymmetrized calculation uses \(K=3\).
- At the left endpoint of a cell \([m,m+1]\), the extra term is harmless because \(b_m=0\).
- The triangular-number intervals cover the full remaining range.
- The \(J=5\) calculation covers \(16\le H\le21\).
- The uniform \(J\ge6\) calculation begins at \(H=21\).

Thus the introduction's cutoff \(\delta=1/22\) is consistent with Section 4. It is not a stale threshold.

#### Independent exact reconstruction

The checks went beyond spot-checking:

| Part of the proof | Independent result |
|---|---|
| Equation (4.6), unsymmetrized \(K=3\) numerator | Exact identity confirmed |
| Bernstein coefficients on \(3\le H\le4\) | All 11 positive |
| Twelve cell polynomials, \(4\le H\le16\) | All 264 coefficients positive |
| Bernstein conversions | Reconstructed the original polynomials exactly |
| \(J=5\) quartic | All five printed coefficients confirmed |
| \(J\ge6\) | All five symbolic factorizations confirmed |
| Finite-sum formulas | Confirmed by polynomial differences and base values |

The \(J=5\) coefficients are exactly
\[
2360,\quad
7500,\quad
\frac{25055}{2},\quad
\frac{254205}{16},\quad
\frac{31115}{2}.
\]

All five factorizations \(\beta_i=\mu_i\pi_i(u)/2880\) agree with reconstruction from the defining formulas, and every coefficient of every printed \(\pi_i\) is positive.

**I found no arithmetic gap in Proposition 3.1.**

## 9. Novelty

I did not identify an earlier statement of any of the three requested conclusions. That is the result of a targeted literature check, not a proof of priority.

The distinctions I would make in the paper are these:

| Result | Assessment |
|---|---|
| \(V\delta_D\ge1/4\) | No earlier result found; this is the principal contribution |
| \(\delta_k\ge1/(4V+|k-D|)\) | No earlier statement found; an immediate and useful consequence of the principal result and the known cubic inequalities |
| \(1/(4V+1)\le\delta_c<2/V\) | No earlier matching statement found; the lower bound follows from the main result, and the upper bound is an elementary tilt argument |

The reciprocal Lipschitz estimate should continue to be presented as a consequence of Hillion–Johnson, rather than as an independent new structural theorem.

Johnson's earlier curvature condition does yield the odds-sum bound you quote. Its translation into your deficit notation is correct; it does not become a variance-only lower bound by reversing a functional inequality. [arxiv.org](https://arxiv.org/pdf/1507.06268?utm_source=chatgpt.com)

The recent Marsiglietti–Melbourne treatment concerns concentration and related consequences of log-concavity. I did not find your central variance–deficit statement there. [arxiv.org](https://arxiv.org/html/2205.08293v4?utm_source=chatgpt.com)

I would keep the restrained claim that you are unaware of an earlier proof that \(\kappa_\star>0\). It accurately identifies the novelty that matters. There is little to gain from making separate priority claims for every short corollary.

## 10. ECP suitability and the strongest "so what?" objection

**The strongest objection is the relationship between the payoff and the proof's length.** A skeptical referee could reasonably ask whether a substantial certificate argument for one local statistic, with the optimal constant still unresolved, warrants publication.

That objection has three concrete components:

- The main constant is not sharp: the current interval remains \(1/4\le\kappa_\star\le1/3\).
- The support-wide bound follows quickly once the main theorem is proved, and the modal upper bound is elementary.
- The CUE example demonstrates a stronger general lower bound than Newton's inequality; it does not establish a new CUE asymptotic or solve a separate random-matrix problem.

Those are matters of significance, and the paper should answer them directly.

### The strongest answer is already mathematical

The paper establishes a **uniform relationship between global variance and the normalized local curvature near the mode**, even when the number of summands and their individual probabilities give poor information about the variance scale.

The balanced variance-one family makes that distinction concrete. As \(m\to\infty\), the lower bounds under discussion behave as follows:

| Bound evaluated at \(D=m+1\) | Behavior in the balanced \(V=1\) family |
|---|---:|
| Newton/ULC bound | \(\sim 2/m\) |
| Bound \(1/(D+1)\) | \(\sim 1/m\) |
| Quoted Johnson bound | \(O(m^{-2})\) |
| Deficit bound obtained directly from Pitman (20) | \(\sim\sqrt2/m\) |
| Theorem 1.1 | At least \(1/4\) |

The entries are consequences of the formulas above and in the manuscript. This is a sharper explanation of the contribution than a general assertion that variance is the natural scale.

Proposition 1.5 then shows why the result contains structural information: ultra-log-concavity permits the corresponding deficit to collapse at bounded variance. The Bernoulli-sum structure rules that behavior out.

Finally, the modal conclusion
\[
\frac15\le V\delta_c<2
\]
and the bound throughout a standard-deviation window show that the theorem has consequences for the central shape of the whole distribution. Those consequences help explain why the first-descent statistic is worth studying.

### How I would weight the motivation

I would give the greatest emphasis to:

1. **The fixed-variance obstruction to the existing bounds.** This makes the gap explicit and testable.
2. **The separation from ULC laws.** This identifies a property that uses more of the real-rooted/Bernoulli-sum structure.
3. **The central-window and modal conclusions.** These explain what the theorem says about distributional shape.

The CUE example is a useful illustration of the difference in scales: \(1/\log N\) versus \(1/N\). It can remain compact. The case for the paper is stronger when the reader understands the structural result before reaching that example.

I would not demand that you prove the optimal constant as a condition of acceptance. Establishing a positive uniform constant is a meaningful first result if, as the literature check suggests, positivity itself was not previously known. The short variance-one binomial example also gives the open constant problem a definite range.

Nor would I raise the current twelve-page length as an automatic eligibility objection: the 2025 ECP editorial report describes manuscripts within twelve, at most thirteen, pages inclusive in ECP format. [imstat.org](https://imstat.org/wp-content/uploads/2025/09/2025-all-reports.pdf?utm_source=chatgpt.com)

**My recommendation is to make the narrow corrections above and submit.** The proof survived independent reconstruction. The most valuable remaining improvement is an exact account of the gap—especially the corrected, stronger Pitman comparison—so that a referee can see immediately what variance alone now guarantees.
