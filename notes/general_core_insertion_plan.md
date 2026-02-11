# General core + s leaves: insertion plan for paper

Goal: add a new theorem (general fixed core with s leaves at a vertex)
without disrupting the existing broom asymptotics narrative.

## Suggested placement

Section: `paper/main.tex`, inside `\section{Broom asymptotics}`.

1) Immediately after the paragraph introducing the broom question (before
   Theorem~\ref{thm:broom-sgeqp}):
   - Insert a short paragraph: “The broom analysis extends to any fixed core
     H with s leaves attached at a vertex v.” This sets the scope.

2) Insert a new theorem environment:
   - Title: “Asymptotic near-miss for leaf attachments”.
   - Hypotheses: fixed tree H with distinguished vertex v; define
     H_s by attaching s leaves to v.
   - Definitions:
     A(x) = I(H − v; x), B(x) = I(H − N[v]; x),
     μ = A'(1)/A(1), m = ceil(μ + 1/2),
     C = (4m+2) − 4μ ∈ [4,8).
   - Conclusion: nm(H_s) = 1 − C/s + O(1/s^2); hence H_s is unimodal for
     all sufficiently large s.

3) Proof sketch block (or full proof if space allows):
   - Note that I(H_s;x) = (1+x)^s A(x) + x B(x).
   - The proof is identical to the broom case, replacing path-specific A,B
     by the fixed-core polynomials; the binomial-ratio expansion is the
     only asymptotic step.
   - Emphasize: the “fixed core” assumption only enters by bounding the
     degree of A and B and treating xB as O(1).

4) After the broom-specific theorem, add a sentence:
   - “The broom asymptotic is the special case H = P_p with v an endpoint.”

## Minimal edits (line-level guidance)

Look for the line:
  “The targeted results raised a natural question: does the broom near-miss ratio …”

Insert the new paragraph and theorem right after that line.

## Narrative framing

Keep the broom result as the concrete worked example, but position the
general-core theorem as the conceptual contribution. This avoids relying
on the spider theorem for novelty and makes the asymptotic method the
paper’s “owned” proof contribution.
