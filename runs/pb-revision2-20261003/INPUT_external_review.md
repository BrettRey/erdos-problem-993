<!-- Pasted by Brett into the Claude Code session on 2026-10-03 without attribution; its source links carry utm_source=chatgpt.com, so it is probably ChatGPT output. Saved by the parent session as the input that prompted the second revision. Every mathematical claim in it was checked in exact arithmetic or by hand before use (scripts/verify_pb_corollaries_20261003.py); external sources it cites were read before citing (HKPV Theorem 7; Meckes-Meckes Proposition 3). Its section headings are kept; display math is kept as LaTeX. -->

Yes. I'd motivate it around a structural question: **how tightly does the global variance of a Bernoulli sum constrain the local shape of its mass function?** Your result supplies a restriction that even ultra-log-concavity doesn't supply. That's a stronger reason to care than the absence of a previously published bound for this particular deficit.

I found several concrete ways to support that motivation:

- Your theorem and the existing recurrence give a lower bound on the deficit at **every support index**, so the first descent can serve as the starting point for a broader result.
- There's an explicit family of ultra-log-concave distributions with bounded variance and first-descent deficit tending to zero. This establishes a substantive distinction between the distributional classes.
- Eigenvalue counts in random unitary matrices give a natural example where your variance bound has scale \(1/\log N\), while Newton's degree bound has scale \(1/N\).
- A short companion upper bound makes the deficit at a mode comparable to \(1/V\), turning the proposed interpretation into a precise statement about scale.

The current introduction does a careful job of explaining what the cited results don't prove. But its progression from definitions to a literature gap leaves the referee to supply much of the reason for studying that gap. The apparently special choice of the first-descent index, followed by a fairly weak geometric decay corollary, makes that harder than it needs to be. variance-scaled-turan-first-descent.pdf

Here's how I'd strengthen the mathematical case.

## 1. Make local shape versus global spread the central question

The normalized Turán deficit measures something quite natural once it's connected explicitly to the curvature of the log mass function. At an interior support index, put

\[ \mathcal C_k = 2\log f_k-\log f_{k-1}-\log f_{k+1}. \]

Then

\[ \mathcal C_k=-\log(1-\delta_k). \]

Thus \(\delta_k\) is a bounded transformation of the discrete curvature of \(\log f\). Ordinary log-concavity says this curvature is nonnegative. Your theorem puts a quantitative lower bound on it near the mode using the distribution's actual variance.

The Gaussian benchmark explains why \(1/V\) is the natural scale. If \(g\) is a normal density of variance \(V\), then its values at three consecutive integers satisfy

\[ 2\log g(k)-\log g(k-1)-\log g(k+1)=\frac1V, \]

and hence

\[ 1-\frac{g(k-1)g(k+1)}{g(k)^2} = 1-e^{-1/V} \sim \frac1V. \]

This suggests a worthwhile finite-sample question: how much of that relationship survives for arbitrary heterogeneous Bernoulli sums, including distributions for which the variance is small and a normal approximation isn't accurate?

Your result supplies a uniform answer on one side. A Poisson–binomial mass function can't have arbitrarily weak central log-concavity while its variance stays bounded above and bounded away from zero. The number of summands, their individual probabilities, and the location of the distribution can vary freely.

That's also why variance matters more than nominal degree here. Adding almost deterministic summands can greatly increase the degree or shift the location while barely changing the distribution's spread. The variance records how much randomness those summands actually contribute.

There's a current literature context for this question. Marsiglietti and Melbourne's 2026 paper investigates the quantitative consequences of different degrees of log-concavity through concentration, moment, and entropy inequalities. Your result fits that broader investigation, with a focus on how actual spread constrains pointwise shape within the more restricted class of Bernoulli sums. [arxiv.org](https://arxiv.org/html/2205.08293v4?utm_source=chatgpt.com)

The distinction from ultra-log-concavity below makes that positioning exact.

## 2. The first-descent theorem gives a bound throughout the support

This is the first addition I'd make. It follows almost immediately from material already in the paper.

Your recurrence (2.4), applied in both directions, gives

\[ \left| \frac1{\delta_{k+1}}-\frac1{\delta_k} \right| \le 1 \]

for adjacent support indices. The endpoint cases follow using \(\delta_0=\delta_n=1\). In words, **the reciprocal deficit can change by at most one per lattice step**.

This is a particularly informative way to present the consequence of the Hillion–Johnson cubic inequalities already used in your proof. Their inequalities contain the extra local regularity that ordinary Newton inequalities don't express. [arxiv.org](https://arxiv.org/pdf/1303.3381?utm_source=chatgpt.com)

Iterating gives

\[ \frac1{\delta_k} \le \frac1{\delta_D}+|k-D|. \]

Your main theorem supplies \(1/\delta_D\le4V\), so

\[ \boxed{ \delta_k\ge \frac1{4V+|k-D|} \qquad (0\le k\le n). } \]

This deserves to be a named corollary immediately after the main theorem.

It explains why the first descent is useful: a bound there determines a lower bound on the curvature across the support. A referee no longer has to regard \(D\) as an isolated index whose importance is asserted by the title.

Several useful formulations follow.

At the rightmost mode \(c=D-1\),

\[ \delta_c\ge\frac1{4V+1}. \]

That's a slightly stronger and more transparent statement than the current \(V\delta_c\ge1/5\), though the latter remains a convenient consequence.

Using Darroch's theorem, \(|c-\mu|<1\), where \(\mu=\mathbb EW\). Therefore \(|D-\mu|<2\), and

\[ \boxed{ \delta_k \ge \frac1{4V+|k-\mu|+2}. } \]

This version doesn't mention the first descent at all. It gives a direct relationship between the variance, distance from the mean, and local log-concavity.

For example, throughout a central window

\[ |k-\mu|\le A\sqrt V, \]

it yields

\[ \delta_k\ge \frac1{4V+A\sqrt V+2}. \]

For fixed \(A\), the lower bound stays on the \(1/V\) scale across a window measured in standard deviations. And because \(-\log(1-\delta_k)\ge\delta_k\), the same lower bound applies to \(\mathcal C_k\) at interior indices.

These are deductions from your theorem and recurrence, rather than additional difficult results. They substantially broaden what the reader sees the theorem doing.

They also clarify the mechanism of the proof. A very small central deficit forces nearby deficits to be small. Near the mode, that in turn forces a broad region of substantial mass, which requires a large variance. The paper's lengthy scalar calculation supplies the explicit constant in that relationship.

## 3. A short upper bound makes the central scale comparison precise

There's a useful companion estimate:

\[ V\delta_c<2 \]

at a mode \(c\). Combined with the preceding lower bound at the rightmost mode, this gives

\[ \boxed{ \frac1{4V+1}\le\delta_c<\frac2V, \qquad V\ge1. } \]

Consequently,

\[ \frac15\le V\delta_c<2. \]

That allows a precise formulation of the motivation: **for Bernoulli sums with variance at least one, the normalized Turán deficit at a mode is comparable to the reciprocal variance, with universal constants.**

Here's a short proof of the upper estimate, using exponential tilting and Darroch's theorem, both already part of your literature context.

Write \(q_k=f_k/f_{k-1}\), and let

\[ f_t(k)=\frac{e^{tk}f_k}{\sum_j e^{tj}f_j}. \]

The tilted law remains Poisson–binomial. Its mean and variance satisfy

\[ \mu'(t)=V(t), \qquad V(t)\ge e^{-|t|}V(0). \]

The latter follows directly from the tilted Bernoulli probabilities, or from \(|V'(t)|\le V(t)\).

At a mode \(c\), let

\[ a=-\log q_c\le0, \qquad b=-\log q_{c+1}\ge0. \]

At tilt \(a\), the adjacent modes are \(c-1,c\); at tilt \(b\), they're \(c,c+1\). Darroch therefore places the first tilted mean in \((c-1,c)\) and the second in \((c,c+1)\), giving

\[ \mu(b)-\mu(a)<2. \]

Meanwhile,

\[ \begin{aligned} \mu(b)-\mu(a) &=\int_a^b V(t)\,dt\\ &\ge V\left(2-e^a-e^{-b}\right)\\ &=V\left(2-\frac1{q_c}-q_{c+1}\right)\\ &\ge V\left(1-\frac{q_{c+1}}{q_c}\right) =V\delta_c. \end{aligned} \]

The last inequality is simply

\[ (1-1/q_c)(1-q_{c+1})\ge0. \]

Under \(V\ge1\), a mode can't be a support endpoint, so the necessary ratios are defined.

I'd use this as a brief comparison lemma, without making a novelty claim for the upper estimate. It lets you distinguish the substantive lower bound from the easier upper bound, while stating a more recognizable result about scale.

One terminological qualification matters: the displayed comparison is for the normalized deficit \(\delta_c\). The exact log curvature is \(-\log(1-\delta_c)\). They become equivalent as the deficits become small; they shouldn't be identified literally throughout the whole finite-variance range.

## 4. Show that ultra-log-concavity itself can't give the conclusion

This is the second addition I'd prioritize. It's stronger than showing that the familiar Newton bound becomes uninformative.

Fix \(q\in(0,1)\), and define a probability mass function on \(\{0,\ldots,2m\}\) by

\[ g^{(m)}_{m+j} = \frac1{Z_m} \binom{2m}{m+j}q^{|j|}, \qquad -m\le j\le m, \]

where \(Z_m\) normalizes the weights.

This distribution is ULC\((2m)\), since division by the binomial coefficients leaves the log-concave sequence \(q^{|j|}\). It also has full support and is strictly log-concave as an ordinary sequence. Symmetry gives its unique mode at \(m\), so its first descent is \(D=m+1\).

At that index, the geometric factors cancel:

\[ \begin{aligned} \delta_D &= 1- \frac{ \binom{2m}{m}\binom{2m}{m+2} }{ \binom{2m}{m+1}^{\,2} }\\ &= \frac{2m+1}{m(m+2)} \longrightarrow0. \end{aligned} \]

Yet the variances remain bounded. For every fixed \(j\),

\[ \frac{\binom{2m}{m+j}}{\binom{2m}{m}} \longrightarrow1, \]

and this ratio is at most one. After centring at \(m\), dominated convergence with the geometric envelope \(q^{|j|}\) gives convergence of the second moments to those of the symmetric geometric law. In particular,

\[ \operatorname{Var}(g^{(m)}) \longrightarrow \frac{2q}{(1-q)^2}. \]

Taking \(q=1/2\), we obtain

\[ V_m\longrightarrow4, \qquad V_m\delta_D\longrightarrow0. \]

So **no positive universal lower bound of your kind holds for ULC distributions in general**, even with full support, strict ordinary log-concavity, and variance converging to four.

That establishes a meaningful boundary between classes. The Bernoulli-sum or real-rootedness assumption contributes a restriction that ULC alone doesn't encode.

The shape of the counterexample also explains what the extra restriction does. ULC permits a pronounced change in slope at the mode followed by a descending shoulder that becomes almost geometric. The reciprocal-deficit recurrence prevents that abrupt change in local curvature for Bernoulli sums.

This example would give a referee a clear answer to a natural objection: "Isn't this just another consequence of ultra-log-concavity?" It also makes the role of the cubic inequalities intelligible before the technical proof begins.

## 5. Give one natural example where variance and degree are far apart

The best concrete example I found is the number of eigenvalues of a Haar random unitary matrix in a semicircle.

Let \(U_N\) be Haar distributed on \(U(N)\), with \(N\) even, and let \(X_N\) count its eigenvalues in the upper semicircle. Although the eigenvalues are dependent, the count has the distribution of a sum of \(N\) independent Bernoulli variables. Meckes and Meckes state this explicitly for region counts in Proposition 3 of their paper on empirical spectral measures. [arxiv.org](https://arxiv.org/pdf/1612.08100?utm_source=chatgpt.com)

The general reason is the Bernoulli representation of determinantal counts: the Bernoulli parameters are the eigenvalues of the restricted kernel. Their variance is

\[ V=\operatorname{tr}(K-K^2) =\sum_i\lambda_i(1-\lambda_i). \]

That makes the distinction between spectral dimension and actual fluctuations explicit. Many spectral contributions can be almost deterministic. [arxiv.org](https://arxiv.org/pdf/math/0503110?utm_source=chatgpt.com)

For the semicircle count, rotation by \(\pi\) exchanges the two semicircles, so \(X_N\) and \(N-X_N\) have the same distribution. Strict log-concavity then gives

\[ c=N/2, \qquad D=N/2+1. \]

The variance grows only logarithmically. From the Diaconis–Evans trace covariance formula, using the Fourier coefficients of the semicircle indicator, one obtains the exact expression

\[ V_N= \frac2{\pi^2} \sum_{\substack{j\ge1\\j\text{ odd}}} \frac{\min(j,N)}{j^2}, \]

and, for even \(N\),

\[ V_N= \frac{\log N+\gamma+\log2+1}{\pi^2} +O(N^{-2}). \]

These expressions are deductions from their covariance formula. [statistics.berkeley.edu](https://statistics.berkeley.edu/sites/default/files/tech-reports/577.pdf?utm_source=chatgpt.com)

For sufficiently large even \(N\), your theorem therefore gives

\[ \delta_D\ge\frac1{4V_N}, \]

with the lower bound asymptotic to

\[ \frac{\pi^2}{4\log N}. \]

Newton's degree-based lower bound at the same index is

\[ \delta_D \ge \frac{N+1}{(N/2+2)(N/2)} \sim\frac4N. \]

The difference is an unbounded factor, in a standard probabilistic model.

For a concrete illustration, evaluating the exact variance formula at \(N=10{,}000\) gives:

| Quantity | Value |
|---|---:|
| Variance \(V_N\) | \(1.16324\) |
| Your lower bound \(1/(4V_N)\) | \(0.21492\) |
| Newton's lower bound | \(0.00039988\) |

Your lower bound is about 537 times larger.

The bound throughout the support strengthens the example further:

\[ \delta_k \ge \frac1{4V_N+|k-N/2-1|}. \]

Throughout any fixed number of standard deviations around the mean, this remains a lower bound of order \(1/\log N\), while Newton's bound remains of order \(1/N\).

The appropriate claim is that this familiar example exposes the scale that degree-based inequalities miss. Establishing an improvement over the best CUE-specific estimates would require a separate comparison. The example already does useful work without that claim: it shows that the distinction between degree and variance arises naturally and can matter greatly.

## 6. Anticipate the two strongest "so what?" objections

### "Doesn't normal approximation already tell us this?"

It suggests the scale, but a usual local central limit theorem doesn't directly establish the finite-sample curvature inequality.

For example, Auld and Neammanee give explicit bounds on individual Poisson–binomial masses with absolute error of order \(V^{-1}\). [Springer Nature Link](https://link.springer.com/article/10.1186/s13660-024-03143-z?utm_source=chatgpt.com) Near the mode, the relevant scales are:

\[ f_k\asymp V^{-1/2}, \qquad \delta_k\asymp V^{-1}, \qquad f_k^2-f_{k-1}f_{k+1}\asymp V^{-2} \]

for the Gaussian benchmark.

Substituting separate mass errors of order \(V^{-1}\) into the Turán determinant allows errors of order \(V^{-3/2}\), larger than the \(V^{-2}\) quantity whose positivity you need to quantify. The cancellation between three adjacent probabilities matters.

Refined expansions controlling those cancellations could establish asymptotic curvature results. Your theorem instead gives an explicit inequality uniformly down to \(V=1\). That's the distinction I'd state.

### "Is the payoff just another tail bound?"

This is where I'd be restrained, because the literature comparison makes a strong tail-based motivation vulnerable.

Your recurrence does improve the present geometric corollary. The bound throughout the support gives

\[ q_{D+j} \le \frac{4V-1}{4V+j-1}, \]

and hence

\[ \frac{f_{D+r}}{f_D} \le \prod_{j=1}^{r} \frac{4V-1}{4V+j-1} \le \exp\!\left[ -\frac{r(r+1)}{2(4V+r-1)} \right]. \]

This reaches the standard-deviation scale; the manuscript's \(\exp[-r/(4V)]\) bound becomes uninformative on that scale as \(V\) grows.

However, Pitman's tilted-ratio inequalities already yield, by a short variance calculation,

\[ \frac{f_{D+r}}{f_D} \le \prod_{j=1}^{r}\frac{V}{V+j-1}. \]

That's often stronger farther into the tail. Its first factor is one, so it doesn't provide your uniform immediate drop after the first descent. But it does mean that Gaussian-scale decay shouldn't carry the novelty claim. [statistics.berkeley.edu](https://statistics.berkeley.edu/sites/default/files/tech-reports/453.pdf?utm_source=chatgpt.com)

The direct pointwise deficit bound, its propagation across the support, and its failure for ULC laws make a better case.

## 7. How I'd reorganize the motivation

I'd give the introduction the following progression.

1. **Start with the relationship between local shape and spread.** Explain the Gaussian \(1/V\) benchmark and ask whether Bernoulli sums satisfy a uniform finite-variance counterpart.

2. **Explain why degree can lose the relevant information.** Keep your existing variance-one example, but move the underlying point earlier. Mention determinantal counts as a natural setting where the discrepancy occurs.

3. **State the first-descent theorem, then immediately state the bound throughout the support.** The latter tells the reader why proving a result at \(D\) has wider significance. Promote the modal estimate from its current position as a remark.

4. **Include the ULC counterexample as a short proposition or example.** This establishes the role of the stronger distributional assumption. It's one of the most persuasive additions available.

5. **Give the semicircle count as a worked application.** A short subsection should suffice. The Bernoulli representation, logarithmic variance, and comparison with Newton are the essential parts.

6. **Then position the result against the literature.** Organize the comparison by what the results control: location of the mode, ratios, maximal mass, local normal approximation, and pointwise curvature. Your existing references can largely stay; the reader will now know why those distinctions matter.

7. **Present the proof mechanism before the scalar machinery.** Display the reciprocal-deficit inequality and explain how small central deficit forces broad mass and therefore large variance. That makes the eventual exact arithmetic calculation serve a visible mathematical purpose.

I'd also broaden the title once the bound throughout the support is stated prominently. Two possibilities are:

- *Variance and local log-concavity of Poisson–binomial laws*
- *Variance bounds for Turán deficits of Bernoulli sums*

The original title could remain appropriate if you prefer to emphasize the main quantitative theorem, but it currently makes the scope look narrower than the consequences warrant.

My assessment is that the strongest revision would be a more convincing short paper on inequalities and distributional shape. Its contribution would be a universal relationship between variance and central log-concavity, a bound extending across the support, and an explicit demonstration that the relationship requires more than ULC. The random-matrix example supplies a setting in which the choice of variance changes the scale of the answer.

**I'd prioritize the bound throughout the support and the ULC counterexample.** Both are short, both sharpen the mathematical content, and together they give the referee a concrete reason to care about the theorem before reaching the proof.
