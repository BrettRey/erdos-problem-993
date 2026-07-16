# Venue decision for the Poisson--binomial article

> Preliminary review artifact. The approved canonical decision is
> `submission/venue-decision-2026-07-16.md`; that record supersedes this note.

## Recommendation

Target **Electronic Communications in Probability (ECP)** first, subject to an early page-fit gate. Draft directly in the ECP class and require the complete paper to fit in at most 12 pages, all inclusive, without moving any load-bearing mathematical explanation outside the article.

Use **Combinatorics, Probability and Computing (CPC)** as the planned alternative if the honest proof and certificate boundary do not fit, or if the theorem grows into a broader package with sharpness, signed corollaries, and applications.

This is a journal-first plan. No public preprint is required for either venue, so lack of mathematics arXiv access is not an obstacle.

## Why ECP is the first target

- ECP publishes short peer-reviewed papers in all areas of probability and explicitly directs papers under 12 pages to ECP rather than its full-length sister journal, EJP.
- The proposed article is currently one intrinsic probability theorem about sums of independent Bernoulli variables, with a signed corollary and an exact finite verification. That is the natural shape of a communication, not yet a full-length paper.
- ECP is immediate open access under CC BY 4.0 and requires no author or publication fee.
- The paper can remain completely independent of Erdős Problem 993. The tree problem belongs in a brief origin note or acknowledgment, not in the title, abstract, or main motivation.

The working title should be close to:

> A variance-normalized curvature bound at the first descent of Poisson--binomial laws

## The early page-fit gate

Before polishing prose, build a complete ECP-template skeleton containing:

1. theorem, notation, and the first-descent interpretation;
2. exact provenance for the Hillion--Johnson and Bobkov--Marsiglietti--Melbourne inputs;
3. deterministic-parameter stripping, support and mode facts;
4. curvature propagation and both modal mass windows;
5. variance synthesis and scalar reduction;
6. the analytic range and exact finite-cell lemma;
7. a mathematically explicit certificate/trust-boundary statement;
8. the signed corollary, constant-status discussion, AI-use disclosure, data/code availability, and references.

The target is 11 pages in ECP style, leaving one page of safety. If this skeleton exceeds 12 pages after ordinary compression, pivot to CPC immediately. Do not omit endpoint cases, theorem-to-certificate interfaces, or enough exact-verification detail for a referee to understand why the computation is rigorous.

## Load-bearing assumptions and falsification conditions

| Assumption | What would make it false? | Response |
|---|---|---|
| The result is naturally a short communication. | The complete proof boundary needs more than 12 ECP pages. | Move to CPC; do not compress away rigor. |
| The exact computation is acceptable as part of a probability theorem. | A referee must trust opaque software or a bespoke environment to understand the proof. | State the finite rational/Bernstein reduction in the article and provide a small independent verifier plus immutable certificate. |
| The constant has enough mathematical meaning. | The paper cannot explain why variance is the right scale or why $1/4$ is non-arbitrary. | Add the finite Poisson-limit obstruction toward the $1/3$ ceiling, or weaken the venue ambition. |
| The theorem stands without the graph project. | The introduction can motivate it only as an auxiliary lemma for Erdős 993. | Reframe around Bernoulli sums, modal curvature, ULC nonimplication, and finite anti-concentration. |
| Novelty is secure enough for submission. | MathSciNet or an expert identifies a prior theorem implying the bound. | Stop and reassess the standalone contribution before submitting. |

An unstated requirement is that a probability referee should be able to audit the mathematical reduction without running the code. The verifier confirms a finite exact claim; it must not substitute for stating that claim precisely.

## Why CPC is the planned alternative

CPC has unusually direct subject precedent: it published Baillon--Cominetti--Vaisman's sharp bound for nonhomogeneous Bernoulli sums and Bobkov--Marsiglietti--Melbourne's discrete log-concave concentration theorem, which supplies one of the present proof's inputs. It also accepts computer-assisted proofs, permits supplementary material and preprints, requires a data-availability statement, and has an explicit policy allowing disclosed substantive AI assistance under human responsibility.

CPC has no public hard page limit. Ordinary subscription publication carries no APC; optional gold open access is expensive and unnecessary. It is therefore the right destination if the natural paper is roughly 12--18 pages or if the constant/sharpness story expands materially.

## Venues not selected

- **Electronic Journal of Probability:** wrong for an 8--12 page paper; reconsider only for a genuinely expanded full-length theorem package.
- **Statistics & Probability Letters:** its six-journal-page cap would likely force essential proof and certificate reasoning out of the article. Its optional SSRN route is useful, but not enough to justify that compression.
- **ALEA:** mathematically suitable, but its current AI policy explicitly mentions AI only for language, grammar, and clarity. This project used AI more substantively, so submission would require prior editorial clarification.
- **Discrete Mathematics:** credible fallback if the work is ultimately framed as a coefficient or discrete-probability theorem, but it reaches the intended probability audience less directly than ECP or CPC.

## Submission gates

Before submission:

1. complete the authenticated MathSciNet pass and a submission-date delta search;
2. obtain one focused expert novelty/proof-interface check;
3. freeze the verifier, exact certificates, replay instructions, dependency record, and hashes in a durable DOI archive;
4. include the project's standard page-one AI-use disclosure and clear human-responsibility statement;
5. run an independent proof/certificate audit;
6. verify the final page count in the actual ECP class.

## Decision rule

Start in ECP format. Stay with ECP only if the article is complete and readable at 12 pages or fewer. Otherwise move once, early, to CPC and use the extra room well.

---
comments:
  c1:
    body: I have no expertise here, but these seem reasonable. I have a skill for
      venue selection or something like that, don't I?
    by: user
    at: 2026-07-16T17:40:40.831Z
