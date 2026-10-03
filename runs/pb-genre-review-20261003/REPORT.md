# Codex genre and style review of the PB paper (2026-10-03)

Saved verbatim by the parent session (Claude Code, Opus 5.5) from the stdout of `tools/codex/bin/codex-readonly` (model printed in the session log: `gpt-6.1-sol`). Prompt: `prompt.md`. The comparator texts were arXiv versions of seven 2026 ECP papers, extracted to a session scratch folder (not committed; listed in the manifest).

---

# Conventions (with evidence)

**The paper already reads like a mathematics paper suitable in register for ECP. Its main weakness is the organization of the introduction, not the author’s phrasing or lack of mathematical sophistication.** It needs clearer separation between announcing the contribution, establishing its significance, and proving subsidiary results.

The amount of exemplification is defensible. The examples establish different aspects of the contribution; they are not elementary exercises added for a novice audience. Their placement and elaboration deserve revision.

The observations below concern the seven supplied arXiv texts. They establish a range of practice, not an ECP house style or binding editorial rules. I have ignored layout and have not re-reviewed the mathematics.

## 1. Introduction: recognizable moves, flexible ordering

The common expectation is that a reader can recover:

- the probabilistic setting;
- the question or limitation in existing knowledge;
- the contribution and its relationship to prior results;
- the main method;
- where to find the details.

There is **no mandatory order for all six moves in your question**.

| Comparator | Organization |
|---|---|
| **Niles-Weed, *Reweighted information inequalities*** | §1 develops the Langevin setting, explains the limitation of existing inequalities, formulates the question, works through a Gaussian-mixture example, and previews the results. General statements wait until §3, after “Preliminaries and prior work.” |
| **Chewi–Pooladian–Zhang, *Stability of the Kim–Milman flow map*** | §1 starts with the stability question, compares existing transport maps, introduces the alternative map, previews bounds and applications, and mentions future extensions. Detailed background and related work occupy §2; the principal theorem appears in §3. |
| **Wang, *Infinite divisibility of α-Cauchy distributions*** | §1 moves from general background and criteria to the particular distribution, the previously posed question, Theorem 1.1, and the proof methods. Further connections follow the theorem. |
| **Samani, *An exponential concentration inequality …*** | A short §1 gives the setting, classical and recent results, the new finite-sample contribution, and a section outline. Results, proof sketch, lemmas, and proofs receive separate sections. |
| **Austern–Mackey, *Bounding Hellinger distance with Stein’s method*** | §1 is exceptionally short and contribution-first: “We introduce a new, explicit bound.” §2 then contains the theorem, interpretive remarks, applications, and related work. |
| **Athaide et al., *A correlation bound …*** | The introduction is explicitly divided into the model, existing correlation results, and a summary of new results and methods. |
| **Cai, *Rates of bulk convergence …*** | §1 gives the setting and question, then substantial definitions and a kernel table before the contribution and strategy. Formal results occupy §2. |

The transferable convention is **functional separation**, not a fixed template. A long opening section becomes easier to use when its parts have distinct purposes. Athaide et al. make those purposes explicit with subsections; Austern–Mackey put most of the material outside the nominal introduction.

Cai is evidence that a lengthy introductory development is possible. Niles-Weed is evidence that an example can precede the general theorem. Neither supports a rule that everything must be stated immediately or that introductions must be almost notation-free.

## 2. Theorems, proofs, remarks, and examples

**Formal results are prominently identified, and their proof destinations are recoverable.** Immediate proof is not obligatory.

- Wang states Theorem 1.1 in §1 and proves its two directions in §§2–3.
- Samani states the principal results in §2, explicitly points to §5 for proofs, and places preparatory lemmas in §4.
- Austern–Mackey assemble results and applications in §2, with dedicated proof sections afterward.
- Athaide et al. state the main and auxiliary results in §1.3 and prove them in §2; an introductory reduction lemma is proved in Appendix A.
- Niles-Weed usually places proofs near their statements in §3.
- Chewi–Pooladian–Zhang combine statements, proofs, remarks, and worked applications within their main-results sections.

**Remarks interpret a result, clarify assumptions, or identify extensions.** Austern–Mackey’s Remarks 2.2–2.4 explain parameter selection and integrability. Niles-Weed’s Remark 3.2 gives a probabilistic interpretation of the construction. Chewi–Pooladian–Zhang’s Remark 3.4 explains how their estimate strengthens a familiar inequality.

Thus mathematical interpretation is welcome. It need not be disguised as another proposition.

**Examples are placed where they perform argumentative work.** They can motivate the question, demonstrate applicability, or explain a hypothesis. They are not required to occupy a separate “Examples” section.

## 3. Motivation and exemplification

The sample contains substantial variation:

- **Niles-Weed** develops a motivating two-Gaussian mixture in §1 and later a metastability application in §4.
- **Chewi–Pooladian–Zhang** work through Gaussian mixtures and potentials that are strongly convex outside a ball in §§3.1.1–3.1.2, then revisit those settings in §4.
- **Austern–Mackey** devote §§2.1–2.2 to applications, including a specific chaos statistic.
- **Wang, Samani, Athaide et al., and Cai** have no comparable sequence of worked application examples. Their model or distribution is already the concrete object of investigation.

Accordingly, **there is no supported numerical quota for examples**. A short paper can have none, one substantial application, or several related instantiations. What matters is whether each illustration changes the reader’s understanding of the theorem.

Motivation is usually mathematical: a known asymptotic result lacks a rate; a familiar inequality fails in an important regime; a known question remains unanswered. Broad social or practical motivation is optional. Wang and Cai include some; Austern–Mackey barely do.

No comparator suggests that a paper needs a plot, simulation, or elementary numerical walkthrough merely to seem accessible.

## 4. Voice, phrasing, sentences, and notation

**“We” is conventional even for a sole author.** Niles-Weed, Wang, Samani, and Cai all use it in mathematical exposition. “I” appears in personal acknowledgements, notably Wang and Samani. Your use of “we” is entirely normal.

The prevailing register combines:

- direct assertions;
- procedural verbs such as *define*, *put*, *fix*, *consider*, and *apply*;
- short sentences identifying the next objective;
- denser sentences when assumptions or qualifications must be stated precisely.

There is no uniform short-sentence style. Chewi–Pooladian–Zhang and Niles-Weed write substantial explanatory paragraphs; Austern–Mackey tolerate very dense notation. Wang and Samani sometimes give more elementary explanation than your paper does.

The useful distinction is between **mathematical density** and **discourse density**. A theorem can carry many symbols because its function is clear. A paragraph becomes harder to read when it simultaneously introduces notation, derives a bound, changes regime, and evaluates novelty.

**Signposting usually names an operation or purpose.** Samani separates a “Proof Sketch” from the lemmas. Athaide et al. explain the order in which the component results will yield the main theorem. Cai’s §3 explains the reduction from a probabilistic distance to an operator norm. These are useful guides, not decorative transition sentences.

## 5. Introducing results from the literature

The comparators generally identify **what a cited result supplies and why it matters**.

- Chewi–Pooladian–Zhang introduce a lemma from Wibisono with a precise theorem reference immediately before using it.
- Wang states Kelker’s results in the form needed for the argument and explains the connection.
- Samani explicitly states the DKW and Laurent–Massart inequalities as proof inputs.
- Athaide et al. distinguish standard random-current machinery from the newer extension to signed interactions.
- Austern–Mackey explain which assumptions or losses distinguish existing CLTs from their result.
- Niles-Weed separates background, antecedents, and the closest comparable contribution in §2.1.

Detailed citations to theorem or equation numbers are normal. A literature discussion should distinguish a source’s stated theorem from a consequence the present author derives from it.

Qualified novelty statements are also normal. Chewi–Pooladian–Zhang and Austern–Mackey use knowledge-qualified formulations; Niles-Weed similarly qualifies the apparent absence of precise statements. Your “We are not aware …” is appropriate.

## 6. Computer-assisted steps and acknowledgements

**This sample does not establish a convention for computer-assisted proofs.** None presents a computational certificate as a necessary component of its proof in the way yours does. Athaide et al. mention numerical algorithms in the literature, but that is a different function.

I can assess your presentation for clarity, but cannot claim that a particular certificate format is standard among these comparators.

One relevant comparison is **Niles-Weed’s brief acknowledgement of ChatGPT and Claude assistance**. It establishes that acknowledging AI assistance is not itself alien to this sample. It does not establish how much disclosure is appropriate when assistance contributed substantially to proof discovery or formalization.

For your paper, the essential reader-facing distinction is between:

1. mathematical arguments printed in the paper;
2. exact computations required to complete those arguments;
3. additional checks and formal verification.

The reader should understand those categories without reconstructing the project’s development history.

## 7. Ending the paper

**None of the seven comparators has a conventional conclusion section summarizing the paper.**

They finish with a proof, a final result or remark, or technical material, followed by references and acknowledgements. Samani and Athaide et al. end their mathematical exposition with appendices; Austern–Mackey end with additional lemmas.

Open questions need not be saved for the end:

- Chewi–Pooladian–Zhang identify future extensions in §1.
- Wang identifies questions left for future research in §1.
- Cai states the conjectured optimality of the rates in §1.

Your unresolved optimal constant belongs naturally beside the bound on that constant. A concluding section would add little.

# Section-by-section comparison

Line numbers below refer to the supplied `.tex`; theorem numbers refer to the supplied PDF text.

## Title and abstract: lines 3–33

The title is restrained, descriptive, and appropriate. It identifies the objects and the relationship being studied without making a promotional claim.

The abstract is recognizably mathematical: setting, quantity, main bound, limitation on the constant, consequence, scope distinction, application, and method. Its density is within the comparator range. Samani’s abstract is also explicit about its inequality; Cai’s includes rates; Niles-Weed’s emphasizes interpretation.

Two details deserve attention:

- **Lines 25–26:** the relative clause “which for a normal density …” has *discrete curvature* as its intended antecedent. A hurried reader could attach it to the deficit or its transform. This is a small but consequential expository ambiguity.
- **Line 33:** “The proof” should say **“The proof of the main theorem.”** The supplement explains that other introductory results are not formalized. The abstract should not invite a broader reading.

The abstract announces many secondary features. Trimming one would be a matter of taste, not a prerequisite for mathematical credibility.

## Introduction, opening and main theorem: lines 46–127

This is the strongest genre fit in the manuscript.

The progression is clear:

1. define the Bernoulli sum and notation;
2. introduce the local deficit;
3. motivate the variance scale through curvature;
4. identify what existing bounds fail to provide;
5. state the theorem.

The normal comparison at lines 76–89 is exactly the kind of concise mathematical motivation that a probability reader can use. The distinction between three adjacent masses and a global variance makes the question intelligible.

The literature paragraph at lines 91–109 is useful because it compares **dependencies of bounds**, rather than merely listing papers. It resembles the functional comparisons in Chewi–Pooladian–Zhang and Austern–Mackey.

“The number of summands is the wrong scale for this question” is emphatic but acceptable. The qualification *for this question* prevents it from becoming an indiscriminate dismissal of existing work.

The main theorem arrives on PDF page 2. That is prompt enough. The problem is not a delayed theorem.

A minor refinement would make Theorem 1.1 more portable: “Under the standing notation above” asks someone reading only the theorem to retrieve the setting. This is common enough to be acceptable, but a short restatement would help.

## Constant, corollary, and modal bound: lines 129–242

All this material is relevant:

- the Bernoulli example establishes why a variance restriction is necessary;
- Proposition 1.2 limits the possible universal constant;
- Corollary 1.3 supplies the support-wide consequence;
- Proposition 1.4 returns explicitly to curvature at the mode.

The open question at lines 174–176 is well placed and appropriately qualified.

**The friction comes from four full proofs interrupting the introductory progression**: the three propositions in §1 and the corollary. Immediate proof is not improper, but it creates a hybrid opening section: introduction, results catalogue, and technical section.

The proof of Proposition 1.4 is especially substantial for this location. It introduces exponential tilting, tilted parameters, means and variances, and a multi-stage integral argument. Its opening explanation at lines 220–224 is good; the problem is placement, not lack of explanation.

Corollary 1.3’s forward reference to the reciprocal bound is acceptable. Samani also states results using quantities developed later. But accumulated forward references make the opening less self-contained than its uninterrupted presentation suggests.

A referee would probably notice the cost of navigating this section before objecting to any individual proof.

## ULC distinction and comparison with previous bounds: lines 244–333

This material protects the contribution from an obvious objection: why does ordinary ultra-log-concavity not already give the result?

Proposition 1.5 is therefore valuable. The observation at lines 277–279 gives it a second purpose by identifying the failed reciprocal regularity. That is an effective mathematical explanation.

The balanced family and Table 1 also deserve to remain. They show that the distinction survives **within Bernoulli sums**, rather than only in a broader class.

The least readable passage is **lines 281–301**. It moves through:

- construction of the balanced family;
- its variance and mode;
- Johnson’s odds-sum bound;
- a bound derived from Pitman;
- a new inverse-tilt notation;
- asymptotic evaluation;
- a limiting distribution and Bessel notation.

This is too much discourse work for one paragraph. The issue is not that an ECP reader cannot follow the formulas; it is that the paragraph keeps changing its task.

Similarly, **lines 324–333** contain a derived ratio estimate and a separate limiting family in compressed prose. The comparison matters, but it reads like a short technical remark embedded in the introduction.

Table 1 is helpful. Its final row distinguishes the actual deficit from lower bounds; that distinction is worth preserving.

## Polynomial consequences and Haar-unitary example: lines 335–376

The brief real-rooted-polynomial consequence is conventional and economical.

Example 1.6 is a good probability application. It demonstrates a change of order, from \(1/N\) to \(1/\log N\), in a recognizable model. It also provides a useful bridge from a Bernoulli-sum inequality to a broader probability audience.

It is, however, **closer to a small application subsection than a casual example**. It develops the compressed kernel, a Gram matrix, strict parameter bounds, a variance identity, asymptotics, symmetry, monotonicity, and an exact applicability threshold.

The comparison at \(N=10^4\) is a useful numerical illustration. The exact threshold \(N=1998\) has a different function: it checks the hypothesis \(V\geq1\). Both are legitimate, but the precise threshold attracts disproportionate attention.

A referee could reasonably ask for a shorter main-text presentation or clearer internal separation. There is no comparator-based reason to remove the application altogether.

## Proof strategy and heuristic: lines 378–411

The strategy paragraph at lines 378–390 is strong. It identifies the propagation mechanism, the mass lower bounds, the variance reduction, and the scalar target.

The heuristic at lines 392–396 should stay. It explains why the scale is plausible and why the proof architecture has a chance of working. Samani devotes an entire section to comparable explanatory work.

**The strategy arrives too late:** the introduction runs through PDF page 6, and the reader encounters several subsidiary proofs before receiving this account of the main proof.

Lines 398–411 then overfill the roadmap. The exact parameter ranges belong naturally in §4. The final sentence introduces \(H,R_r,L_r,b_r,\lambda_r\), none yet developed at this point. It gives a ledger of possible losses before the reader understands the quantities being compared.

A short statement that the argument does not determine the optimal constant would suffice here.

## §2, recurrence and standing reductions: lines 413–581

This section reads like standard mathematical exposition.

Its opening states its purpose. The reductions are explicit, the external inequalities are attributed precisely, and the reciprocal bound is extracted transparently. Lemma 2.1 is followed immediately by its proof, then by the mass estimates it enables.

The treatment of deterministic summands in §2.1 repeats the introduction’s informal reduction more formally. That is acceptable: previewing a reduction and later justifying it is normal practice.

The extension to differences of independent Bernoulli sums is short and useful. It does not need an additional worked example.

The left- and right-mass calculations are appropriately explicit. Shortening them merely to appear more sophisticated would harm the exposition.

## §3, variance bounds and scalar reduction: lines 583–642

This is the clearest structural section.

It tells the reader what the mass bounds are for, identifies the literature input, states the remaining scalar proposition, and proves the main theorem conditional on that proposition.

That proof order is entirely conventional. Wang and Samani likewise organize proofs around intermediate inputs rather than following statement order mechanically.

The short proof of the maximal-mass bound at lines 605–616 earns its space: it keeps a central input understandable and makes the reduction easy to inspect. “For completeness, here is a short proof” is normal mathematical phrasing.

No substantial restructuring is needed here.

## §4.1, compact range: lines 644–703

The monotonicity argument is explained rather than hidden behind software. The interval partition and exact-arithmetic requirement are also stated.

The opening parenthesis at lines 670–672 introduces \(S,T\) from the next subsection while explaining why a particular simplification fails. It is accurate as signposting, but mildly awkward in dependency order.

The twelve decimal ratios at lines 692–696 are less useful than the partition rule and exact checks. They give reassurance about margins, but are not themselves the certificate. A sentence giving the minimum margin would communicate the relevant fact more efficiently.

The essential exact computations are in the supplement. The pending archive link is consequently a substantive submission issue, not merely unfinished copyediting.

## §4.2, large range: lines 705–826

The main chain is well motivated:

- symmetrize the weights;
- factor the quadratic expression;
- introduce a positivity test;
- use triangular cells to obtain nonnegative lower bounds;
- exhibit positive coefficients.

The explanation at lines 757–760 of why the cells start at triangular numbers is particularly good.

The polynomial displays are dense, but their function is obvious: they are inspectable positivity certificates. Cai’s operator expansions and Samani’s derivative appendix show that explicit algebra is compatible with this genre.

This subsection ends by closing the proposition and identifying the supplementary expansion checker. A separate conclusion is unnecessary.

## Supplement and acknowledgements: lines 828–895

The supplement gives an unusually careful account of what is required and what is additional. Its distinction between the scalar proposition, the deduction of the main theorem, and unformalized results is valuable.

The difficulty is packaging: **one long paragraph serves as archive catalogue, proof-dependency statement, verification report, and formalization inventory**. A reader should be able to find the required files without processing every supplementary assurance.

The acknowledgements are the most conspicuously unusual part of the paper. They record model labels, effort settings, dates, successive reports, contributions to individual results, and corrections to previous drafts. This reads like a development and audit log.

Some detail may be warranted by the extent of assistance. Nevertheless, naming AI-generated documents as “referee reports” risks blurring their status with journal peer review. The sentence about the author’s responsibility is appropriate and should remain.

# Prioritized changes

## An ECP referee or editor would likely ask for these

These are predictions about reader-facing clarity and submission readiness, not claims that every ECP referee would demand the same remedy.

| Priority | Location | Concrete change | Why |
|---|---|---|---|
| **1** | Lines 831–836; §4.1 | Replace `ZENODO-DOI-PENDING` with an accessible archive reference and identify the required checking programs unambiguously. | Necessary computations must be inspectable. Unlike stylistic preferences, the present placeholder prevents the described proof materials from being located. |
| **2** | Lines 129–411 | Give the opening material a visible hierarchy. Keep the setting, gap, main theorem, principal consequences, and a short strategy together; place subsidiary proofs and detailed comparisons in clearly named subsequent subsections or a section. | Athaide et al., Austern–Mackey, and Samani make the functions of their opening material easy to distinguish. Your theorem is early, but its surrounding section becomes difficult to navigate. |
| **3** | Lines 281–301 and 324–333 | Separate the construction of the balanced family from the derivation and asymptotic evaluation of the cited bounds. Display the derived Pitman estimate and define its inverse-tilt notation locally. | The comparators explain literature inputs before using them. These passages currently combine too many operations and require the reader to reconstruct their hierarchy. |
| **4** | Abstract, lines 25–26 and 33 | State the normal-curvature comparison in a separate clause or sentence; specify that the **main theorem’s proof** is formalized. | The abstract should distinguish curvature from its transform and avoid implying that all announced results are formalized. |
| **5** | Lines 831–859 | Divide the supplement description into required proof computations, additional checks, and formalization scope. Put optional or earlier verification routes after the required route. | The reader needs a clear proof-dependency account. The present information is useful but buried in a catalogue. |

## Matters of taste, with a strong editorial preference where indicated

| Priority | Location | Concrete change | Why |
|---|---|---|---|
| **6** | Lines 867–894 | Shorten the main-text AI acknowledgement to substantive contributions and responsibility. Put dates, effort settings, successive checking reports, and correction history in an archive provenance document. Call synthetic reviews “AI-generated reviews.” | Niles-Weed supplies a concise disclosure precedent. Your fuller contribution requires fuller disclosure, but the chronological log is disproportionate to the paper’s mathematical exposition. |
| **7** | Lines 378–411 | Move the conceptual proof sketch earlier. Retain the heuristic; replace the parameter-by-parameter certificate roadmap with one short description of the exact and symbolic stages. | Samani and Cai give readers the strategy before sustained technical execution. The current roadmap introduces unnecessary future notation. |
| **8** | Example `ex:cue`, lines 343–376 | Retain the Bernoulli representation, variance growth, modal index, order comparison, and one numerical illustration. Move the exact first-applicable \(N\) calculation to the supplement, or isolate it in a short final note. | Applications are conventional, but their detail should serve the application’s main point. The exact threshold distracts from the change of scale. |
| **9** | Lines 689–702 | Replace the list of twelve rounded margins with their minimum and the explanation of how the intervals are checked; retain full values in the archive if useful. | Readers need the reduction and certificate, rather than a long display of diagnostic output. |
| **10** | Theorem `thm:main`, lines 119–127 | Restate briefly that \(W\) is a finite independent Bernoulli sum and \(f\) its pmf. | A slightly more self-contained theorem is easier to quote and consult. The current standing-notation formulation is acceptable. |
| **11** | Lines 670–672 | Say simply that the compact range requires the unsymmetrized weights; move the \(H=7/2\) countercheck until after \(S,T\) have been defined. | This removes a small forward dependency without changing the explanation. |

# What to leave alone

- **Keep the mathematical voice.** “We,” “put,” “hence,” and “it remains to prove” are conventional. There is no pervasive phrasing problem, and nothing in the main mathematical prose warrants the verdict “a mathematician would never write this.”
- **Keep the curvature motivation and the distinction between local deficit and global variance.** They give the question a clear mathematical purpose.
- **Keep the sharpness family, the ULC counterexample, the balanced Bernoulli family, and the Haar-unitary application.** They respectively address the attainable constant, the necessary class structure, the weakness of prior bounds within the target class, and probabilistic applicability. Their functions are distinct.
- **Keep Table 1.** It makes the comparison easier to inspect.
- **Keep the proof heuristic and explicit mass calculations.** Accessibility does not require suppressing explanations.
- **Keep precise source locators and qualified novelty claims.**
- **Keep the visible positivity certificates.** Their large coefficients are not a stylistic defect.
- **Keep the open constant question beside its bracket.**
- **Do not add a conclusion, more elementary examples, generic applications rhetoric, or simulations merely to satisfy an imagined ECP convention.**
- **Do not conceal computer assistance or substantive AI contributions.** Improve their organization and terminology.

# What you could not judge

I did not assess correctness, optimality, novelty against the full literature, certificate validity, or agreement between the Lean statements and the manuscript. Those fall outside the requested editorial review.

The seven comparators do not establish ECP’s current requirements for AI disclosure, supplementary archives, or formal verification. My suggestions on those subjects concern clarity and proportion, not a verified journal policy.

I assessed content and exposition from the LaTeX and extracted PDF text. The PDF text establishes that the introduction extends through page 6; I did not evaluate typography or visual layout.

Finally, the comparison supports a limited editorial prediction: **the manuscript’s register is already credible; readers are more likely to ask for a better-organized opening and a clearer presentation of verification materials than for more motivation or more examples.**
