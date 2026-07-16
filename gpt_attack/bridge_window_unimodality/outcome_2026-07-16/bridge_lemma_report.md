REFUTED(T3)

# Bridge lemma report: zero locations do not force window unimodality

## 1. Scope of the verdict

The literal tree statement BRIDGE is not refuted here. What is refuted, in the precise sense of sub-target T3, is the proposed implication from pointwise zero geometry to coefficient-window unimodality.

For the fixed degree parameter \(d=24\), I construct normalized polynomials \(P_s\) with all coefficients positive such that:

1. \(\alpha_s=\deg P_s\to\infty\);
2. their zeros obey the complete common degree-\(25\) zero-free open set, a unique dominant negative zero, the exact \(0.90\) positive-axis sector, and the quantified cardioid-cusp upper envelope;
3. nevertheless \(P_s\) has a valley at indices \(2s-1,2s,2s+1\), with \(2s/\alpha_s\to 1/2\).

Thus this is not a fixed-prefix artifact. The construction uses \(\alpha\to\infty\) structurally and places the valley proportionally inside the requested window.

## 2. The zero hypothesis being refuted

Let \(\mathcal Z_{24}\) be the set of complex numbers that occur as zeros of independence polynomials of finite graphs of maximum degree at most \(25\), and let

\[
 \mathcal U_{24}=\mathbb C\setminus\overline{\mathcal Z_{24}}
\]

be the full common zero-free open set, in the notation of Bencs--Buys--Peters. For a polynomial with a unique smallest-modulus zero \(-\beta\), put

\[
 r(z)=\frac{|z|}{\beta},
 \qquad
 \phi(z)=\pi-|\operatorname{Arg}z|
\]

for nonreal \(z\), where \(\operatorname{Arg}z\in(-\pi,\pi]\). Define the cardioid-cusp profile

\[
 F(r)=\frac{(2(r-1))^{3/2}}6
 \qquad (r\geq 1).
\]

The hypothesis \(\mathsf H_{24}\) is:

- H1 (dominance): the smallest-modulus zero is unique, simple, and negative real;
- H2 (all rigorous common zero-free geometry): no zero lies in \(\mathcal U_{24}\);
- H3 (sector): every nonreal zero satisfies \(|\operatorname{Arg}z|>9/10\);
- H4 (cusp upper envelope): every nonreal zero satisfies \(\phi(z)\leq F(r(z))\).

H4 is the precise direction of Conjecture B' in the supplied note real_collar_conjecture_2026-07-16.md. In the measured band \(1\leq r\leq2\), its profile is exactly the prediction

\[
 \phi\leq \frac{(2(r-1))^{3/2}}6
\]

recorded in data.md. Requiring this displayed formula for every \(r\geq1\), rather than only in the sampled band, makes the hypothesis fully explicit. For robustness, the construction will also be shown to satisfy the additional capped lower collar

\[
 \phi(z)\geq
 g(r(z))
 :=
 \frac{\{2(\min\{r(z),2\}-1)\}^{3/2}}6.
\]

H2 contains all the rigorous K4 regions. Precise supplied-source references are Bencs--Csikvári--Srivastava--Vondrák, Theorems 1.1 and 8.1 in arxiv_2204_04868.txt, and Bencs--Buys--Peters, Theorems 1.5 and 1.7 and Lemma 4.5 in arxiv_2111_06451.txt. The finite-size sector is deliberately fixed at \(9/10\), although Bencs--Buys--Peters, Remark 5.5, explains why no such fixed sector can hold for all trees asymptotically.

## 3. The construction

Set

\[
 Q(x)=(1+x)^6+20x
     =1+26x+15x^2+20x^3+15x^4+6x^5+x^6
\]

and, for every integer \(s\geq3\), set

\[
 \boxed{P_s(x)=Q(x)(1+x^2)^{2s}.}
\]

Then

\[
 P_s(0)=1,
 \qquad
 \alpha_s=\deg P_s=4s+6.
\]

### Lemma 1 (the base is the certified split-graph control)

\(Q\) is the independence polynomial of \(K_{20}\vee E_6\), whose maximum degree is \(25\).

#### Proof

An independent set in the join cannot use vertices from both sides. It may contain at most one of the \(20\) clique vertices, while it may contain an arbitrary subset of the six vertices in the edgeless side. Hence

\[
 I(K_{20}\vee E_6;x)=(1+x)^6+20x=Q(x).
\]

A clique vertex has degree \(19+6=25\), and an edgeless-side vertex has degree \(20\). Thus the maximum degree is \(25=d+1\) for \(d=24\). \(\square\)

### Lemma 2 (a hand-checkable dominant-root collar)

\(Q\) has exactly one zero in \(|z|<2/5\). That zero is simple and negative real, say \(-\beta\), with \(0<\beta<1/20\). Every other zero of \(Q\) has modulus greater than \(2/5\).

#### Proof

On \(|z|=2/5\),

\[
 |(1+z)^6|
 \leq (7/5)^6
 =\frac{117649}{15625}
 <8
 =|20z|.
\]

Rouché's theorem therefore says that \(Q(z)=20z+(1+z)^6\) and \(20z\) have the same number of zeros, counted with multiplicity, in the disk: exactly one.

Moreover,

\[
 Q(0)=1>0,
 \qquad
 Q(-1/20)=-\frac{16954119}{64000000}<0.
\]

The intermediate value theorem supplies a negative real zero in \((-1/20,0)\). It must be the unique zero in the disk, and its multiplicity is one. The remaining zeros lie strictly outside \(|z|=2/5\). \(\square\)

### Lemma 3 (the exact \(0.90\) sector, without root-finding)

Every nonreal zero \(z\) of \(Q\) satisfies

\[
 |\operatorname{Arg}z|>9/10.
\]

#### Proof

By conjugation it is enough to exclude a zero \(z=re^{i\theta}\) with \(0<\theta\leq9/10\). Since \(9/10<\pi/2\), both \(z\) and \(1+z\) lie in the open right half-plane. Write

\[
 1+z=\rho e^{i\psi},
 \qquad 0<\psi<\theta.
\]

The equation \(Q(z)=0\) is

\[
 (1+z)^6=-20z.
\]

Both \(6\psi\) and \(\theta+\pi\) lie strictly between \(0\) and \(2\pi\), so equality of arguments gives

\[
 6\psi=\theta+\pi.
\]

Since \(\psi<\theta\), this implies \(\theta>\pi/5\). Elementary sine-rule geometry in the triangle with vertices \(0,1,1+z\) gives

\[
 r=\frac{\sin\psi}{\sin(\theta-\psi)},
 \qquad
 \rho=\frac{\sin\theta}{\sin(\theta-\psi)}.
\]

Taking moduli in \((1+z)^6=-20z\) therefore yields

\[
 \sin^6\theta
 =20\sin\psi\,\sin^5(\theta-\psi). \tag{3.1}
\]

Now

\[
 0<\theta-\psi
 =\frac{5\theta-\pi}{6}
 <\frac14,
\]

using \(\theta\leq9/10\) and the elementary bound \(\pi>3\). Therefore the right side of (3.1) is less than

\[
 20(1)(1/4)^5=\frac5{256}.
\]

On the other hand, \(\pi/5<\theta<\pi/2\), so

\[
 \sin^6\theta>\sin^6(\pi/5)>\frac1{27}.
\]

The last inequality follows from

\[
 \sin^2(\pi/5)=\frac{5-\sqrt5}{8}>\frac13,
\]

because \(7>3\sqrt5\). Finally \(1/27>5/256\), contradicting (3.1). \(\square\)

### Lemma 4 (a hand-checkable cusp certificate and an exact stronger collar)

Every nonreal zero of \(P_s\) has \(r>8\), and consequently satisfies H4. Exact rational root counts strengthen this to \(r>25\) and show that every negative-axis deviation exceeds \(\pi/6\); hence the additional capped lower collar displayed in Section 2 also holds.

#### Proof

By Lemma 2, a nonreal zero of \(Q\) has modulus greater than \(2/5\), while \(\beta<1/20\), so its dominance ratio is greater than \(8\). The added roots \(\pm i\) have ratio \(1/\beta>20\). Therefore, for every nonreal zero,

\[
 F(r)
 =\frac{(2(r-1))^{3/2}}6
 >\frac{14^{3/2}}6
 >7
 >\pi
 \geq\phi.
\]

The middle inequality uses \(\sqrt{14}>3\), and the next uses \(\pi<4\). This proves H4 entirely by hand.

For the advertised stronger collar, exact Gaussian-rational Cauchy-index root counts put one zero of \(Q\) in each of the following pairwise-disjoint intervals or rectangles:

\[
\begin{gathered}
 (-3,-2),\qquad (-1/25,-2/51),\\
 [-2,-3/2]\times[-2,-3/2]i,
 \qquad [-2,-3/2]\times[3/2,2]i,\\
 [1/5,1/4]\times[-5/4,-1]i,
 \qquad [1/5,1/4]\times[1,5/4]i.
\end{gathered}
\]

The six exact counts sum to \(\deg Q=6\), so the boxes are exhaustive. In particular,

\[
 \frac2{51}<\beta<\frac1{25},
\]

while every nonreal zero has modulus greater than \(1\). Thus \(r>25\).

For the left-half-plane pair,

\[
 \tan\phi
 =\frac{|\operatorname{Im}z|}{|\operatorname{Re}z|}
 \geq\frac{3/2}{2}
 =\frac34
 >\frac1{\sqrt3},
\]

so \(\phi>\pi/6\). For the right-half-plane pair, \(\phi>\pi/2\). The extra roots of \(P_s\) are \(i\) and \(-i\), each with multiplicity \(2s\); they have \(r=1/\beta>25\) and \(\phi=\pi/2\).

Since \(r>25\), the optional capped lower-collar value is

\[
 g(r)=\frac{2^{3/2}}6=\frac{\sqrt2}{3}<\frac\pi6;
\]

the last inequality follows from \(\pi>3>2\sqrt2\). Hence the additional lower collar holds as well.

The root counts are certified by the replayable script [bridge_lemma_exact_check.py](/tmp/claude-agent-output/bridge_lemma_exact_check.py). It invokes exact SymPy Poly.count_roots calls on Gaussian-rational rectangles; it never invokes nroots or performs floating-point root-finding. \(\square\)

### Lemma 5 (the complete \(d=24\) common zero-free set is respected)

No zero of \(P_s\) lies in \(\mathcal U_{24}\).

#### Proof

By Lemma 1, every zero of \(Q\) is a zero of the independence polynomial of a graph of maximum degree \(25\). Hence no such zero can lie in the common zero-free set.

It remains to check the added roots \(\pm i\). Bencs--Buys--Peters define

\[
 \mathcal C_d=
 \left\{
 -\frac{d^d u}{(u+d)^{d+1}}:|u|<1
 \right\}.
\]

Their Lemma 4.5 in arxiv_2111_06451.txt proves

\[
 \mathcal U_d\subseteq\mathcal C_d.
\]

For \(d=24\) and \(|u|<1\),

\[
 \left|
 \frac{24^{24}u}{(u+24)^{25}}
 \right|
 <\frac{24^{24}}{23^{25}}
 <\frac12. \tag{5.1}
\]

Here is an exact check of the last inequality. Termwise comparison of the binomial expansion gives

\[
 \left(1+\frac1{23}\right)^{23}
 <\sum_{k=0}^{\infty}\frac1{k!}<3,
\]

where the final bound follows from

\[
 \sum_{k=0}^{\infty}\frac1{k!}
 =2+\sum_{k=2}^{\infty}\frac1{k!}
 <2+\sum_{k=2}^{\infty}\frac1{2^{k-1}}
 =3.
\]

Consequently,

\[
 \frac{24^{24}}{23^{25}}
 =\frac1{23}\left(1+\frac1{23}\right)^{24}
 <\frac1{23}\left(3\cdot\frac{24}{23}\right)
 <\frac4{23}<\frac12.
\]

Thus \(\mathcal U_{24}\subset\{|z|<1/2\}\), whereas \(|i|=|-i|=1\). The added roots are outside it. \(\square\)

### Lemma 6 (positive coefficients and a proportional window valley)

Every coefficient of \(P_s\) is a positive integer. If

\[
 p_{s,k}=[x^k]P_s(x),
\]

then

\[
 p_{s,2s-1}>p_{s,2s}<p_{s,2s+1}. \tag{6.1}
\]

Moreover,

\[
 2s+1\leq
 L_s=\left\lceil\frac{2\alpha_s-1}{3}\right\rceil,
 \qquad
 \frac{2s}{\alpha_s}\longrightarrow\frac12.
\]

#### Proof

The factor \((1+x^2)^{2s}\) has positive coefficients at every even degree from \(0\) through \(4s\). Since every coefficient of \(Q\) from degree \(0\) through \(6\) is positive, their convolution is positive at every degree from \(0\) through \(4s+6\): for \(k\leq4s\), use the \(Q\)-term of degree \(0\) or \(1\) according to the parity of \(k\); for \(k>4s\), use the \(Q\)-term of degree \(k-4s\).

Put

\[
 B_j=\binom{2s}{s-j}
 \qquad (0\leq j\leq3).
\]

At \(b=2s\), direct convolution gives

\[
\begin{aligned}
 p_{s,b}&=B_0+15B_1+15B_2+B_3,\\
 p_{s,b-1}&=26B_1+20B_2+6B_3,\\
 p_{s,b+1}&=26B_0+20B_1+6B_2.
\end{aligned}
\]

Therefore

\[
 p_{s,b-1}-p_{s,b}
 =-B_0+11B_1+5B_2+5B_3>0,
\]

because

\[
 \frac{B_1}{B_0}=\frac{s}{s+1}>\frac1{11}.
\]

Also \(B_0\geq B_1\geq B_2\geq B_3>0\), so

\[
\begin{aligned}
 p_{s,b+1}-p_{s,b}
 &=25B_0+5B_1-9B_2-B_3\\
 &\geq15B_0+5B_1>0.
\end{aligned}
\]

This proves (6.1).

Finally,

\[
 \alpha_s=4s+6,
 \qquad
 L_s=\left\lceil\frac{8s+11}{3}\right\rceil,
\]

and \(2s+1\leq(8s+11)/3\). Hence

\[
 a=2s-1<b=2s<c=2s+1\leq L_s.
\]

Also

\[
 \frac{b}{\alpha_s}
 =\frac{2s}{4s+6}
 \longrightarrow\frac12,
 \qquad
 \frac12-\frac{b}{\alpha_s}
 =\frac3{4s+6}.
\]

\(\square\)

### Theorem 7 (T3 refutation)

The sequence

\[
 P_s(x)=\bigl((1+x)^6+20x\bigr)(1+x^2)^{2s},
 \qquad s\geq3,
\]

has all positive integer coefficients, satisfies \(\mathsf H_{24}\), has degree \(\alpha_s=4s+6\to\infty\), and has a window valley at

\[
 2s-1<2s<2s+1\leq L_s.
\]

Therefore the pointwise zero-location hypotheses H1--H4 do not imply BRIDGE.

#### Proof

The roots of \(P_s\) are the six roots of \(Q\), together with \(i\) and \(-i\), each repeated \(2s\) times. H1 follows from Lemma 2 because the added roots have modulus \(1\). H2 is Lemma 5. H3 follows from Lemma 3 for the \(Q\)-roots and from \(|\operatorname{Arg}(\pm i)|=\pi/2\) for the added roots. H4 is Lemma 4. Lemma 6 supplies positivity, the proportional window valley, and \(\alpha_s\to\infty\). \(\square\)

## 4. Exact instances

For \(s=3\),

\[
 P_3(x)=Q(x)(1+x^2)^6
\]

has \(\alpha=18\), \(L=\lceil35/3\rceil=12\), and ascending coefficient sequence

\[
\begin{aligned}
(&1,26,21,176,120,516,336,856,546,880,\\
 &546,576,336,236,120,56,21,6,1).
\end{aligned}
\]

Thus the exact valley is

\[
 p_5=516>p_6=336<p_7=856,
\]

with \(7\leq L=12\).

For \(s=4\), the corresponding exact valley is

\[
 p_7=2064>p_8=1338<p_9=3108,
\]

with \(\alpha=22\) and \(L=15\).

The replay script checks the six exact root boxes, positivity, and these inequalities for \(3\leq s\leq100\). Lemma 6, not this finite run, proves the coefficient claim for every \(s\geq3\).

## 5. The missing graph/tree-only input

The construction isolates the failure. Since \((1+x^2)^{2s}\) has no linear term,

\[
 [x]P_s(x)=26
\]

for every \(s\), even though \(\alpha_s=4s+6\to\infty\). Equivalently, if \(z_1,\ldots,z_{\alpha_s}\) are the roots with multiplicity, logarithmic differentiation of

\[
 P_s(x)=\prod_{j=1}^{\alpha_s}\left(1-\frac{x}{z_j}\right)
\]

at \(x=0\) gives

\[
 -\sum_{j=1}^{\alpha_s}\frac1{z_j}
 =\frac{P_s'(0)}{P_s(0)}
 =26. \tag{7.1}
\]

The \(2s\) copies of \(i\) and the \(2s\) copies of \(-i\) cancel exactly in this inverse-root moment.

For a normalized independence polynomial of any graph,

\[
 [x]I(G;x)=|V(G)|=n\geq\alpha(G).
\]

Thus an actual tree sequence with \(\alpha\to\infty\) necessarily has

\[
 -\sum_j z_j^{-1}=n\geq\alpha.
\]

A successful bridge therefore needs at least multiplicity-sensitive control of an inverse-root moment, or an equivalent tree recurrence or coefficient condition. Radial tightness is not enough: the added roots here remain fixed at modulus \(1\); it is their growing multiplicity and exact angular cancellation that defeat every location-only hypothesis.

## 6. An additional proved positive step: saddle localization

Although T3 is already complete, the following step is useful for any repaired bridge and explicitly covers the high-degree-hub regime.

### Lemma 8 (degree-uniform tree saddle localization)

Let \(T\) be any tree, let \(\alpha=\alpha(T)\geq3\), and put

\[
 I(x)=\sum_{j=0}^{\alpha}i_jx^j,
 \qquad
 \mu(x)=\frac{xI'(x)}{I(x)}.
\]

For each \(1\leq k<\alpha\), there is a unique \(x_k>0\) with \(\mu(x_k)=k\). If

\[
 k\leq L=\left\lceil\frac{2\alpha-1}{3}\right\rceil,
\]

then

\[
 x_k\leq256.
\]

If \(k\leq2\alpha/3\), then \(x_k\leq64\).

#### Proof

Under the hard-core probability law

\[
 \Pr_x(J=j)=\frac{i_jx^j}{I(x)},
\]

we have \(\mu(x)=\mathbb E_xJ\) and

\[
 \frac{d}{d\log x}\mu(x)=\operatorname{Var}_x(J)>0.
\]

The variance is positive because the support contains both \(0\) and \(1\). Thus \(\mu\) increases continuously from \(0\) to \(\alpha\), proving existence and uniqueness.

Let \(f(t)=\log I(e^t)\). Convexity and \(f'(t)=\mu(e^t)\) give, for \(t>0\),

\[
 \mu(e^t)\geq\frac{f(t)-f(0)}t.
\]

There is at least one maximum independent set, so \(f(t)\geq\alpha t\). A tree is bipartite, hence \(\alpha\geq n/2\), and \(I(1)\leq2^n\); therefore

\[
 f(0)=\log I(1)\leq n\log2\leq2\alpha\log2.
\]

Consequently,

\[
 \mu(e^t)
 \geq\alpha-\frac{2\alpha\log2}{t}. \tag{8.1}
\]

For \(\alpha\geq3\), checking the three residue classes modulo \(3\) gives

\[
 L\leq\frac{3\alpha}{4}.
\]

Taking \(t=8\log2\) in (8.1) yields

\[
 \mu(256)\geq\frac{3\alpha}{4}\geq L,
\]

and monotonicity gives \(x_k\leq256\). Taking \(t=6\log2\) yields \(\mu(64)\geq2\alpha/3\), proving the sharper assertion.

The proof uses only bipartiteness and coefficient positivity, so it includes hub-bouquet trees without a hidden bounded-degree or low-hub assumption. \(\square\)

## 7. Verification map

| Item | How a human checks it |
|---|---|
| Lemma 1 | Count independent sets in \(K_{20}\vee E_6\) by hand; the two vertex degrees are \(25\) and \(20\). |
| Lemma 2 | Check the displayed rational Rouché inequality and the exact value of \(Q(-1/20)\). |
| Lemma 3 | Follow the argument equation, the two sine-rule identities, and the rational comparison \(1/27>5/256\). No root approximation is used. |
| Lemma 4 | H4 itself follows by hand from Lemma 2 and \(F(8)>7>\pi\). For the optional stronger collar, run the linked exact script: SymPy Poly.count_roots uses exact Gaussian-rational Cauchy-index/argument-principle calculations, and the six counts sum to degree \(6\). Then check the angular inequalities from the rational boxes. |
| Lemma 5 | Use Bencs--Buys--Peters, Lemma 4.5 in arxiv_2111_06451.txt, then check (5.1) by the displayed binomial estimate. The named K4 regions are supported by Bencs--Csikvári--Srivastava--Vondrák, Theorems 1.1 and 8.1 in arxiv_2204_04868.txt, and Bencs--Buys--Peters, Theorems 1.5 and 1.7. |
| Lemma 6 | Expand the three coefficients using \(B_j=\binom{2s}{s-j}\); each difference has the displayed strictly positive lower bound. |
| Theorem 7 | Combine Lemmas 2--6 and the factorization of \(P_s\). |
| Exact instances | Multiply \(Q\) by the integer binomial coefficients of \((1+x^2)^6\), or run the linked big-integer convolution. |
| Missing moment | Logarithmically differentiate the root factorization at \(0\); for a graph, count singleton independent sets. |
| Lemma 8 | Check the Gibbs variance identity, the convex secant inequality, and \(L\leq3\alpha/4\) in the three residue classes of \(\alpha\) modulo \(3\). |

No claim in this report depends on floating-point root-finding. Every numerical constant in the refutation is effective.

## 8. What I tried that failed

### Far-negative padding

The first obstruction was \(Q(x)(1+x/N^3)^N\), which preserves the split-graph valley and even the \(0.90\) sector. Its valley remains at \(k=2\), so I rejected it because the packet separately declares the prefix \(k<\varepsilon\alpha\) out of scope.

### Ordinary saddle/local-CLT route

The entropy argument proves the effective saddle bound in Lemma 8, but a standard relative-error local central limit theorem is too coarse to control the discrete second difference of adjacent saddle masses. The difficulty is real at critical hub saddles: Bencs--Buys--Peters, Corollary 5.4 and Remark 5.5 in arxiv_2111_06451.txt, locate the cusp boundary and the loss of a uniform positive-axis sector.

### Spider route

Li--Yang--Zhang--Zhang, Theorem 3.1 in li_yang_zhang_zhang_2025_symmetric_function_log_concavity.txt, proves that every spider has a strongly log-concave independence polynomial, so a structural T1 conclusion is already available. Their proof uses two-row Schur positivity rather than K3--K5 zero geometry and does not use \(\alpha\to\infty\), so I did not relabel it as the requested spectral bridge.

### Root-gap route

Prakash--Sharma, Theorem 1.1 in arxiv_2510_09197.txt, isolates the smallest zero only at the stated scale \((\beta/n)^{O(n)}\). The exponent and constants are non-explicit in the supplied \(O(n)\) statement, so that statement does not furnish the packet's required effective constants; in any case, the separation does not control the multiplicity-sensitive inverse-root moment exposed by (7.1).

## 9. Self-grade

Grade: full success at T3, not T1/T2/T4. The construction is exact, effective, fixed-\(d\), has a proportional non-prefix valley, satisfies a stronger-than-sampled sector and a quantified cusp, and identifies the missing graph/tree-only inverse moment. It refutes zero-location-only BRIDGE hypotheses; it does not claim a counterexample among trees.
