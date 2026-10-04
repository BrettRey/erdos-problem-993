# Pasted review: proportion, completeness, Remark 5.1 (2026-10-03)

Pasted into the session by Brett on 2026-10-03 with no instruction attached. The source links carry `utm_source=chatgpt.com`, so it is treated as an AI-generated ChatGPT review (model unknown). Saved verbatim below by the parent session (Claude Code, Opus 5.5). The parent's triage is in DECISIONS.md (2026-10-03, "Pasted ChatGPT review").

---

Yes. I’d cut Remark 5.1 from the main text, provided the added references do identifiable work. It’s a sensible economy, and it doesn’t leave a hole in the argument. I’d also move some of the detailed AI-use history into the supplement; that would make substantially more room.

My overall judgment of variance-local-log-concavity-poisson-binomial.pdf is that it reads like a normal specialist probability paper. The mathematical exposition is coherent, its technical density is mostly appropriate, and the contribution is clearly stated. It needs a few adjustments of proportion and completeness, rather than a major reorganization.

## What cutting Remark 5.1 would lose

There are two distinct Pitman comparisons here, and the numbering makes them easy to confuse:

- **Equation (5.1), the preceding discussion, and Table 1** compare a lower bound on the normalized deficit with your main result. They substantiate the introduction’s claim that existing bounds can degenerate while variance stays fixed. Keep these.
- **Remark 5.1** compares Pitman’s direct adjacent-mass-ratio bound with your subsidiary consequence (1.6). It shows that this consequence remains informative in a situation where that particular Pitman estimate approaches one.

The remark therefore has mathematical value, but it’s peripheral to the paper’s main contribution. No subsequent proof depends on it. Removing it leaves the central comparison with Pitman intact, together with the evidence for the advantage of variance scaling.

If you preserve its derivation in the supplement, it’ll still be available to a referee who asks whether (1.6) follows from a familiar ratio estimate. I wouldn’t feel obliged to retain a numbered remark in the article just to anticipate that question.

It’s also fairly short—about seven lines in the PDF. Cutting it alone might not pay for three or four complete bibliography entries plus additional prose.

## Does the paper read normally?

Yes, particularly in the things that matter most.

The introduction establishes a precise question: whether local curvature near the mode admits a uniform variance-scaled lower bound. Theorem 1.1 gives a clear answer, and Proposition 1.2 immediately indicates the remaining uncertainty about the constant. The whole-support corollary and the estimate at the rightmost mode then explain why the first-descent formulation has broader significance.

Two parts do especially useful work for the reader. Proposition 1.6 identifies a meaningful boundary of the result: ultra-log-concavity doesn’t suffice. The balanced family and Table 1 give a controlled comparison with previous bounds. Together, these supply much more convincing motivation than a general assertion that variance is a natural parameter.

The proof overview on page 4 is also good. It explains the mechanism:

A small deficit forces nearby deficits to be small; that forces enough mass to remain near the maximum over a sufficiently wide region; the resulting spread imposes a variance lower bound.

The following Gaussian heuristic helps the reader understand the scale before encountering the products and auxiliary weights. I’d protect those paragraphs when making cuts. The opening explanation of the exponential-tilting argument in §5.2 is another example of effective mathematical exposition: it says what the calculation accomplishes before performing it.

I don’t see a need for an elementary background section, a separate conclusion, or an expanded catalogue of applications. The paper already gives readers enough to understand its question and follow the strategy.

## Where the density could be improved

### The introduction is busy

Section 1 has nine numbered theorem, proposition, corollary, remark, and example items. That’s recognizable mathematical practice, but it gives several relatively minor observations almost the same visual prominence as the main results.

Remarks 1.4 and 1.8 could readily become ordinary prose. Remark 1.7 could be integrated with the discussion following Proposition 1.6. This would improve the flow without removing mathematical content.

The abstract is similarly comprehensive: it includes the curvature interpretation, functional-inequality motivation, the main bound, the obstruction to a better constant, the support-wide consequence, the ULC counterexample, the Haar example, and formal verification. It’s understandable, but it asks the reader to absorb almost the whole paper’s inventory at once. I’d treat simplifying that inventory as optional polish, not a requirement for submission.

### The technical density in §4 is justified

The most demanding pages are the scalar certification argument, especially §4.2. I wouldn’t describe that as inappropriate density. The polynomial expansions and Bernstein coefficients are part of the proof’s evidence, and you explain the main choices: why symmetrization helps, why it isn’t used in the small range, and why the intervals begin at triangular numbers.

There’s a reasonable division of labour between the conceptual explanation on page 4 and the explicit certificate in §4. Readers can understand the design before deciding how closely to inspect the arithmetic.

I’d be more reluctant to remove those coefficients or the explanation of the interval checks than to remove Remark 5.1. Their presence makes the computer-assisted argument inspectable.

### The AI disclosure is the clearest departure from ordinary article style

The approximately 300-word disclosure on page 13 partly reads as a history of successive review sessions. Details such as which review prompted a rounding correction or compared the exposition with papers in the target journal are useful provenance, but they don’t all need to occupy space in the article.

The substantive disclosure should remain clear: AI systems contributed to proof search, proposed results and arguments, assisted with computation and drafting, and produced formalizations. The scope of formal verification and your responsibility also belong in the paper.

The sequence of sessions and the editorial changes they prompted fit naturally in the archive. Moving that material would make the article more conventional while preserving the full account of how the work was done. I’d retain the specific acknowledgment that Aristotle’s proof supplied the monotonicity argument, since that’s a substantive mathematical contribution.

## Which background references would earn their space?

The current eleven references aren’t, by themselves, a sign of inadequate scholarship. You already cite a recent Poisson–binomial survey, the results used in the proof, and the closest bounds you discuss. I’d add references to strengthen particular connections, without aiming at a numerical target.

Three plausible additions have distinct roles:

- **Saumard and Wellner, *Log-concavity and strong log-concavity: A review* (2014).** This would give readers an entry into the broader discrete and continuous curvature literature behind the opening motivation. It covers connections with concentration and functional inequalities, so it fits naturally around the discussion of quantitative log-concavity. [arxiv.org](https://arxiv.org/abs/1404.5886?utm_source=chatgpt.com)
- **Liggett, *Ultra Logconcave Sequences and Negative Dependence* (1997).** This would give the ULC discussion a foundational reference beyond the sources from which you quote the particular inequality. Its central convolution result also helps situate ULC as a substantive class of distributions. It’s useful if you add a sentence that develops that context. [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0097316597927904?utm_source=chatgpt.com)
- **Brändén, *Unimodality, log-concavity, real-rootedness and beyond*.** This would support the combinatorial setting behind Remark 1.8. Its treatment explicitly connects real-rootedness with normalized coefficient log-concavity and discusses probabilistic consequences and Stirling-polynomial examples. That gives it a clear job in the paper. [arxiv.org](https://arxiv.org/html/1410.6601v1?utm_source=chatgpt.com)

I wouldn’t add a fourth unless you’ve identified another specific gap. Two well-placed background citations could do more for the reader than four citations appended to a broad opening sentence.

For further work, you already have the right question: the value of κ⋆. I’d develop the existing sentences after (1.4) to ask whether the variance-one binomial sequence is asymptotically extremal. That adds something to the statement that the optimal constant is unknown, without requiring a new concluding section.

I’d be more cautious about using the final sentence to promise consequences for log-Sobolev inequalities. Those inequalities provide legitimate motivation, but the paper doesn’t derive a new one. In particular, the local comparison with Johnson’s quantity at D doesn’t itself establish a global curvature constant. A question about how to use the spatially varying bound in functional inequalities would be appropriately open.

## Two concrete things to fix

**The supplement needs an explicit access point.** On page 12, the supplement description ends with the bold phrase *Code and exact certificate data*, but I couldn’t find an embedded archive link or an archival identifier in this PDF. If the files accompany a journal submission separately, the editor may have everything needed. For the circulating paper, though, readers should be able to reach the exact certificate and formalization from the article itself.

The description of the formalization’s scope is otherwise commendably precise: it distinguishes the main theorem from the other introductory results. Keep that distinction when shortening the surrounding material. I’ve assessed the description here; I haven’t run the supplementary Lean code.

**There’s a small boundary issue in the introductory Johnson comparison.** Page 2 writes Johnson’s quantity at D as

\[
(f_D/f_{D+1})\delta_D.
\]

That expression needs a restriction to \(D<n\), or an explicit convention for a zero denominator. The case \(D=n\) really can occur with \(V\ge1\): for \(\operatorname{Bin}(5,0.7)\), the variance is 1.05, the mode is 4, and the first descent is at 5, so \(f_{D+1}=0\). Your main argument handles this boundary through \(\delta_n=1\); the explanatory comparison just needs the same care.
