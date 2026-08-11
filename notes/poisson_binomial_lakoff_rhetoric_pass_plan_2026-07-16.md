# Poisson--binomial paper: conceptual-metaphor and rhetoric pass
Target: `paper/poisson_binomial/main.tex`
## What this pass assumes
1. “Lakoff-style” means inspecting the conceptual mappings that organize the prose, not imitating Lakoff's sentence style or adding decorative imagery.
  
2. The paper's primary model is an indexed probability profile: the masses rise to a modal peak, the first strict descent follows the largest mode, and local curvature records the change in adjacent ratios. This spatial model fits the mathematics and should remain.
  
3. Conventional mathematical terms such as _mode_, _first descent_, _curvature_, _left/right_, _compact range_, and _geometric envelope_ carry technical content. Replacing them wholesale would make the paper less precise.
  
4. The pass is prose-only. It won't change a statement, inequality, hypothesis, proof dependency, citation, certificate claim, or novelty qualifier.
  
5. The voice will remain lean and journal-appropriate, with sparse mathematical _we_ and natural contractions.
  

The approach fails if a revision changes a mathematical entailment, makes the prose more figurative, or replaces a standard technical expression with an idiosyncratic one. Those changes will be rejected even if they sound smoother.
## Pass 1: make the conceptual system coherent
Keep the useful profile/topography model, but remove optional crossings into unrelated source domains.

- Replace the economic metaphor in “variance discounts nearly deterministic Bernoulli factors” with the literal claim that nearly deterministic factors contribute little to variance.
  
- Replace “furnish the propagation input” with the specific mathematical action: the cubic inequalities yield the neighboring-deficit recurrence.
  
- Introduce “modal mass windows” once as lower bounds for masses on either side of the modal index. After that definition, the shorthand is earned.
  
- Replace force language such as a window “forced by” propagation with direct dependency language: the recurrence implies the bounds, and the two variance inequalities become incompatible under the contrary assumption.
  
- Replace “certificates handle the compact interval” with “certificates prove positivity on the compact interval.”
  
- Replace “on the lower-bound side,” “fall,” and other dispensable directional language with literal statements where the spatial model adds nothing.
  
- Retain _first descent_, _curvature_, _left/right_, _mode_, and _geometric envelope_: each corresponds directly to an indexed sequence or a standard bounding object.
  
## Pass 2: strengthen the argument's rhetoric
- Keep the page-one order: limitation of known results, variance-scaled question, affirmative theorem.
  
- Compress the three defensive optimality paragraphs after Proposition 1.2 into one economical account of what the example proves, where the lower-bound proof loses sharpness, and why the paper claims only the bracket.
  
- Make the signed extension's introduction direct: translation extends the result to differences of Bernoulli sums.
  
- Remove the standalone sentence that merely restates the geometric-envelope corollary.
  
- Break the long related-work paragraph into units organized by what each result controls: mode location/stability, adjacent ratios, and single-atom bounds.
  
- Rewrite the proof roadmap as an explicit implication chain: cubic recurrence → two-sided modal lower bounds → variance/max-atom squeeze → scalar positivity, with exact certificates only on the finite range.
  
- Replace “For completeness” and “For clarity” frames with direct topic sentences.
  
- Rename “Variance synthesis and scalar reduction” to “Variance bounds and scalar reduction”; _synthesis_ adds an architectural/chemical metaphor without helping the reader.
  
- Make the audit status agentive and plain: the paper won't be submitted until an independent human proof and certificate audit is complete.
  
## Preservation and verification gates
- Compare every changed sentence against the surrounding displayed mathematics.
  
- Run the house-style and terminology checks, then a read-only proofreading audit.
  
- Rebuild with the ECP class and inspect all affected pages.
  
- Require no undefined references or citations, no new box warnings, and a clean `git diff --check`.
  
- Re-run the mathematical test suite and exact certificate checker if any wording touches a computational claim.
