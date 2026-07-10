# GPT-5.6 Ultra Breakthroughs

Date: 2026-07-10

## Headline Result

The finite signed reserve problem is solved at the target quarter scale.

Let $W$ be any finite Poisson-binomial law with variance $V\ge1$, and let

\[
D=\min\{k\ge1:f_k<f_{k-1}\}
\]

be its first supported strict descent. Then

\[
V\left(1-\frac{f_{D-1}f_{D+1}}{f_D^2}\right)\ge\frac14.
\]

Writing $R_-=f_D/f_{D-1}<1$ and $R_+=f_{D+1}/f_D$, we have

\[
1-R_+\ge 1-\frac{R_+}{R_-}
=1-\frac{f_{D-1}f_{D+1}}{f_D^2}.
\]

Consequently,

\[
V\left(1-\frac{f_{D+1}}{f_D}\right)\ge\frac14.
\]

Every finite signed Bernoulli sum $X-Y$ becomes an ordinary
Poisson-binomial law after a deterministic shift, so this is the full
nonterminal finite signed theorem. If the largest mode is the upper support
endpoint, the next descent has zero mass; the effective quotient is undefined,
but the raw reserve is immediately one.

The proof has passed an independent line audit and an exact replay of every
computer-assisted certificate. That is sufficient for the present theorem
claim. A Lean/Aristotle pass would be an additional audit, especially valuable
because the endpoint convention and descent indexing are delicate.

## Why the Proof Works

The successful route is neither perturbative continuity nor an iteration over
reflected Bernoulli factors.

Hillion--Johnson's cubic inequalities imply that the normalized Turán
curvature at neighboring indices cannot grow too quickly. If the curvature
$\delta$ at the first descent were small, those inequalities would force an
explicit two-sided window of masses comparable to the modal mass.

Two independent variance bounds then trap the modal mass:

\[
V\ge p_{\max}^2A(\delta),
\qquad
V\ge\frac{p_{\max}^{-2}-1}{12}.
\]

The theorem reduces to the scalar lower bound

\[
A(\delta)\ge\frac{3+\delta}{4\delta^2}.
\]

That scalar inequality is proved by:

- one exact asymmetric certificate on $3<H<4$;

- twelve exact symmetric certificates on $4\le H\le16$;

- an analytic Bonferroni reduction and nonnegative Bernstein coefficients for
  $H\ge16$.

The original finite certificates use Sturm root counts. An independent
conversion subsequently found a simpler Lean-facing route: after translating
each of the 13 finite cells to $0\le t\le1$, all 275 full-degree Bernstein
coefficients are strictly positive. Exact reconstruction passed for every
cell, and the generated Lean 4.28 identity file compiles.

The endpoint convention $\delta_0=\delta_n=1$ is essential. An independent
audit found a truncated abstract sequence that defeats the recurrence if the
terminal cubic inequality is omitted; the genuine Poisson-binomial endpoint
law rules it out.

## Supporting Theorems Found En Route

Before the universal proof emerged, three bounded results were proved and
audited independently:

1. The complete Skellam family satisfies
   $V\Delta_{\mathrm{eff}}>1/4$.

2. Subtracting one reflected Bernoulli preserves the quarter bound, with an
   exact plateau-selection threshold.

3. Subtracting two reflected Bernoullis preserves the quarter bound through a
   six-term curvature identity and exact $5\times5$ Bernstein matrices.

A separate agent also closed the three-reflected-factor case directly through
a ten-term curvature identity and finite Bernstein tensors. That proof is now
redundant as a dependency. It should not be sent for mechanical formalization
until its identity and tensors have first been written to disk.

## Exact Verification

| Artifact | Exact/high-precision checks | Result |
| --- | ---: | --- |
| Universal theorem replay | 1 asymmetric cell, 12 compact cells, full analytic Bernstein tail | pass |
| Lean-friendly finite replay | 13 cells, 275 exact Bernstein coefficients | 275 positive |
| Universal finite scan | 11,320 PBD vectors | 0 failures |
| Skellam harness | 2,691 rows, 305 refined transitions | 0 failures |
| One reflected Bernoulli | 47,850 exact rows | 0 failures |
| Two reflected Bernoullis | 97,488 exact rows, 584,928 component bounds | 0 failures |
| Project tests | 55 tests | pass |

The universal replay is a proof certificate. The finite grids are
falsification checks only.

## What This Closes

For the hub product term

\[
A(x)=(1+x)^sQ(x),
\qquad
V_A=\frac{s}{4}+V_Q,
\]

the theorem gives, for $s\ge4$,

\[
\max_{j\ge D_A}\frac{A_{j+1}}{A_j}
\le 1-\frac1{4V_A}
=1-\frac1{s+4V_Q}.
\]

Thus the general signed bridge and the product-term reserve criterion in issue
#5 are closed.

## What Remains Open

The full hub-bouquet polynomial is

\[
I(x)=A(x)+B(x),
\qquad
B(x)=xR(x).
\]

The new theorem supplies a concrete post-descent perturbation budget:

\[
B_{j+1}-B_j\le\frac{A_j}{4V_A}.
\]

A proof must also prevent $A+B$ from descending prematurely before $D_A$.

For fixed arms, $B$ has fixed support and should reduce to a finite threshold
certificate once $s$ is large. For growing broom arms, $B$ overlaps the mode
and remains the real analytic obstruction.

This does not prove:

- the full hub-bouquet theorem or issue #5;

- the general Case-B mode bound for multi-hub trees;

- Conjecture A or PNP;

- Erdős Problem 993.

The submitted manuscript was deliberately left unchanged.

## Recommended Next Move

The best next proof target is the fixed-arm $A+xR$ perturbation theorem,
because the new reserve gives an explicit budget and the remaining support
comparison is finite. In parallel, the universal theorem itself is now
substantial enough to consider as a standalone note or a deliberate manuscript
revision.

Aristotle is useful at the formalization stage. The higher-value target is the
universal theorem, because it now carries the argument and has a durable proof
note and replayable certificate. A sensible staged formalization is:

- first formalize the endpoint-aware curvature-propagation lemma, the forced
  mass window, and the raw-from-effective corollary;

- then replay the finite scalar certificates in Lean rather than importing
  them as assumptions;

- finally formalize the Hillion--Johnson cubic input for a fully end-to-end
  theorem.

The first minimal Lean 4.28 packet at
`formalization/pb_effective_drop_aristotle/` is now complete. Aristotle project
`8d59a353-5c7b-4071-b837-9ab7bf561be3` filled all nine proof holes covering
curvature propagation, endpoint exclusion, the crossing ratio, and the
raw-from-effective implication. Independent local replay builds without
`sorry` or new axioms. This milestone remains conditional on the normalized
Hillion--Johnson recurrence and therefore does not yet constitute an
end-to-end formalization of the universal PB theorem.

## Principal New Files

- `notes/literature/universal_poisson_binomial_effective_drop_2026-07-10.md`

- `scripts/verify_universal_pb_effective_drop.py`

- `results/universal_pb_effective_drop_certificate_2026-07-10.json`

- `notes/literature/universal_pb_finite_bernstein_replacement_2026-07-10.md`

- `scripts/verify_universal_pb_finite_bernstein.py`

- `results/universal_pb_finite_bernstein_certificate_2026-07-10.json`

- `formalization/pb_effective_drop_aristotle/`

- `notes/literature/skellam_effective_drop_theorem_2026-07-10.md`

- `notes/literature/one_reflected_bernoulli_effective_drop_2026-07-10.md`

- `notes/literature/two_reflected_bernoulli_effective_drop_2026-07-10.md`
