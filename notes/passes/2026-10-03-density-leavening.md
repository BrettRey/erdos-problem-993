# Density leavening: PB paper
<!-- SUMMARY: density-leavening pass on variance-local-log-concavity-poisson-binomial.tex, judged against the redundancy pass · status: 9 edits applied, paper now 13 pp (Brett: fine) · updated: 2026-10-03 -->

Manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`.
Before: sha256 `b5bc1d8285c5521a…` (commit `74bf546`). After: `5258b16330fec5fc…`.
PDF after: `08abe1f8f589640a…`, **13 pages**, with references [7]–[11] on page 13.
Brett, 2026-10-03: "13pp is fine." The journal's limit, as recorded in
`submission/venue-decision-2026-07-16.md` l. 44, is "within 12 pages, at most
13". The 12-page cap was this project's own safety margin.

Brett's brief: "The point is to make the paper readable and accessible." The
reader assumed throughout is a probabilist who doesn't work on log-concavity or
Bernstein certificates.

## Findings by signature

### 1. Terminology load in the opening

- **D7** l. 61–64. Before Theorem 1.1, the reduction sentence named "the
  deficits defined below" and "the first-descent index" before either was
  defined. It now reads: "…summands with $p_i=1$ only translate the law,
  which changes nothing below except by a shift of index (Section 2.1)." The
  sentence is shorter and has no forward-referenced terms; §2.1 keeps the
  details. This is the only edit inside the cold-read opening. Checked against
  the venue record's four tasks for pages 1–2 (l. 39): define the deficit,
  explain the variance scale, state the 1/4 theorem and the 1/3 obstruction,
  distinguish prior work. The edit touches none of them.

### 2. Abstraction with no observable

None found that leavening would fix. The paper's claims are inequalities, each
stated with its constants. The motivating Gaussian benchmark, the binomial
example and Table 1 already give concrete cases.

### 3. Missing intermediate steps

This is where the paper asked most of its reader.

- **D1** (proof outline). The target $(3+\delta)/(4\delta^2)$ appeared
  without explanation. Added: "Together they give $V(1+12V)\geq A(\delta)$.
  Since $V(1+12V)$ increases with $V$ and equals $(3+\delta)/(4\delta^2)$ at
  $V=1/(4\delta)$, the theorem reduces to…". Check: $V\geq M^2A$ and
  $M^{-2}\leq1+12V$ give $V(1+12V)\geq A$. Sympy confirms
  $(1/(4\delta))(1+3/\delta)-(3+\delta)/(4\delta^2)=0$.
- **D1b** (heuristic paragraph). It said the masses within distance about
  $\delta^{-1/2}$ are comparable to $M$, without saying why. Added: "the lower
  bounds on $f_{c\pm r}/M$ are about $e^{-r^2\delta/2}$ when $r\delta$ is
  small", and the closing step now uses the exact $V(1+12V)\geq A(\delta)$.
  Check: $\log R_r=r\log a+\sum_{j<r}\log(1-j\delta)=-r^2\delta/2+O(r\delta)+O(r^3\delta^2)$.
  At $\delta=10^{-3}$, $r=30$ this gives $\log R_r\approx-0.469$ and $\log L_r\approx-0.470$ against
  $-r^2\delta/2=-0.45$. The paragraph stays labelled "Heuristically".
- **D2** (Proposition 1.4 proof). The tilting argument came with no statement
  of what it does. Added, phrased as the idea so it can't be read as part of
  the argument: tilting moves the mean at a rate equal to the tilted variance;
  the tilts at which $c-1,c$ and then $c,c+1$ are the modes differ by
  $-\log(1-\delta_c)$; and by Darroch's theorem the mean moves by less than 2
  between them. Check: $t_+-t_-=\log q_c-\log q_{c+1}=-\log(q_{c+1}/q_c)=-\log(1-\delta_c)$
  by (2.2), with $c\leq n-1$ since $D\leq n$.
- **D3** (Pitman (20) paragraph). It went from "the arguments differ by $2/m$"
  and "slope $\sqrt2$" straight to "the bound is $\sqrt2/m$". Added: the
  tilted mean is $m+1+O(1/m)$ at $t=\operatorname{arsinh}1$; since
  $\log\theta$ is the inverse of the tilted mean, the two values of
  $\log\theta$ differ by $\sqrt2/m+O(m^{-2})$, "and so does the bound". Check
  (sympy series): mean $=m+1+(1-\sqrt2)/m+O(m^{-2})$, slope $=\sqrt2+O(1/m)$.
  The bound $1-e^{-\Delta}=\Delta+O(\Delta^2)$.
- **D6** (Example 1.6). "Hence $V_N=\operatorname{tr}\Gamma-\operatorname{tr}\Gamma^2$"
  left out the reason. It now reads "Hence $V_N$, the sum of $p(1-p)$ over the
  eigenvalues $p$ of $\Gamma$, is $\operatorname{tr}\Gamma-\operatorname{tr}\Gamma^2$".
- **D9** (§4.1). It didn't say why $3<H\leq4$ keeps the unsymmetrized bounds.
  Added: "because the symmetrized bound fails there: $ST<Q(H)$ at $H=7/2$,
  for example." Check (exact rationals): $ST-Q=-3237589/1882384$ at
  $H=7/2$. It's also negative at $H=3.01$ and $3.2$, and positive at $3.8$
  and $4$.
- **D5** (§4.2). "The elementary product bound" named an inequality the
  reader might not recognize, and the triangular-number cells had no
  explanation. Now: "the inequality $\prod_s(1-x_s)\geq1-\sum_sx_s$ for
  $x_s\in[0,1]$ gives (4.10), where $\lambda_r\geq0$ because
  $r(r+1)/2\leq J(J+1)/2\leq H$; this is why the cells begin at triangular
  numbers." Check: $x_s=s/H\in[0,1]$ since $s\leq J\leq K<H$, and
  $\sum_{s\leq r}s/H=r(r+1)/(2H)$.

### 4. Stacked premodifiers

None to unpack. "Maximal-mass bound", "first-descent index", "pairwise
variance identity", "normalized Turán deficit" and "degree-four Bernstein
coefficients" are terms of art for this readership, each defined or standard.

### Reorders (no content change)

- **D4** In the Pitman ratio paragraph, the paper's own result
  ($f_{D+1}/f_D\leq1-1/(4V)$ from (1.7)) now comes first, so the reader knows
  why Pitman's ratio bound is being derived. Before, the point arrived in a
  final "By contrast" sentence.
- **D8** The §4 roadmap ("Section 4.1 treats… Section 4.2 treats…") moved
  from between the definition of $H$ and its consequence ("Then $K=\ldots$")
  to the end of the §4 preamble, where it introduces the subsections.

## Judged against the redundancy pass

The test is the redundancy registry's own criterion: would "says again what
the paper has already said" cut this sentence?

| Edit | Restates anything? | Verdict |
|---|---|---|
| D1 | No. $V(1+12V)\geq A$ and the value at $V=1/(4\delta)$ appear nowhere else; §3 argues by contradiction without showing where the constant comes from. | Keep. |
| D1b | No. The decay rate is new. Its last clause uses the D1 inequality instead of re-saying "the two bounds force", so it's shorter than a restatement. | Keep. |
| D2 | **Partly.** A statement of the idea previews the proof that follows. It doesn't repeat any sentence in the paper, and it gives the reader the shape of an argument that is otherwise a chain of inequalities. | Keep, for readability. This is the one place where leavening and redundancy pull in opposite directions, and readability decides it. |
| D3, D5, D6, D9 | No. Each adds a step or reason the text left out. | Keep. |
| D7 | Removes words. | Consistent with redundancy. |
| D4, D8 | Reorders only. | Neutral. |

None of the three redundancy cuts (R1–R3) is undone: no restated closer, κ⋆
consequence or doubled closure has returned.

## Re-verification of passes made stale by this round

Scoped to the diff, as for the previous round
(`2026-10-03-post-edit-reverification.md`).

| Pass | Check on this diff | Result |
|---|---|---|
| adversarial-cold-read | Only D7 falls in the extracted opening (`coldread.py extract --prompt-only`); it's checked against the venue record's four tasks above. | Round-9 verdict stands. No new readers. |
| reader-pass, coherence-cohesion | Read every hunk in the rendered PDF, including the reordered Pitman paragraph and the new position of the §4 roadmap. | Better; no breaks. |
| level-category-audit | New predications: "tilting moves the mean", "$\log\theta$ is the inverse of the tilted mean", "the symmetrized bound fails". Each is a mathematical relation correctly attributed. The tilt idea is marked as the idea, not as a step of the proof. | No issue. |
| terminological-hygiene | New terms: "tilted variance" (already used in the proof), "symmetrized bound" (matching "symmetrization" in the outline). | No issue. |
| numbers-audit | Every new number or identity is checked above: $7/2$, $\sqrt2/m$, $m+1+O(1/m)$, $(3+\delta)/(4\delta^2)$, $e^{-r^2\delta/2}$. | Correct. |
| negative-claims-audit | No negative claim added. "The symmetrized bound fails there" is a computed counterexample, not a claim about the literature. | No issue. |
| source-reread, external-review-triage | No cited claim changed. D3 adds an inference inside the existing Pitman paragraph; the citation still covers equation (20) only. No review point is affected. | No issue. |
| build-integrity | pdflatex/bibtex/pdflatex ×2: 0 undefined references, 0 overfull boxes, 13 pages. | Clean. |
| proofread, de-ai-prose, editorial-scar-tissue | Read every hunk. There's one semicolon clause in D5 and no colon reveals. The only colon (D9, "fails there: $ST<Q(H)$ at $H=7/2$") introduces the witness. Nothing reads as a correction of earlier wording. | No issue. |
| redundancy, lakoffian-metaphor, paragraph-opening-audit, rhetoric-and-humour | Redundancy is judged above. No new frames except "tilting moves the mean", which is literal. D4 changes one paragraph opening: it now opens with the paper's own bound. | No issue. |

## Side effect outside the manuscript

D5, D8 and D9 change §4's wording, so the Aristotle packet's "§4 verbatim" in
`formalization/pb_scalar_inequality_aristotle/PROOF_CONTEXT.md` differs further
in wording. The mathematics is unchanged.
