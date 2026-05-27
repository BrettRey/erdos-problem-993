# AxiomProver Scan for Erdős #993

Date: 2026-05-27

## Purpose

Record the scan prompted by the recent AxiomProver/Axiom Math papers and
extract proof-search lessons for the tree independence-polynomial unimodality
project.

## Sources Checked

- Axios article, 2026-05-26:
  `https://www.axios.com/2026/05/26/axiom-ai-math-journal`
- Axiom public papers endpoint, checked 2026-05-27:
  `https://axiom-backend-0w7x.onrender.com/api/papers?sort%5B0%5D=Date:desc`
- Axiom Fel conjecture Lean repository:
  `https://github.com/AxiomMath/fel-polynomial`
- Relevant arXiv entries listed by the Axiom papers endpoint:
  - `https://arxiv.org/abs/2605.21718`
  - `https://arxiv.org/abs/2604.25246`
  - `https://arxiv.org/abs/2604.13238`
  - `https://arxiv.org/abs/2603.29970`
  - `https://arxiv.org/abs/2603.23928`
  - `https://arxiv.org/abs/2602.05090`
  - `https://arxiv.org/abs/2602.03716`
  - `https://arxiv.org/abs/2602.03722`

## What Transfers

The useful pattern is not broad theorem discovery.  It is:

```text
narrow natural-language theorem
-> exact algebraic/combinatorial formal target
-> Lean/mathlib certificate or counterexample
-> human proof note around the certificate
```

This matches the Axiom Fel repository structure, where a task statement and
natural-language theorem are separated from Lean proof files.  That is a better
fit for this project than asking a prover to attack Erdős #993 globally.

## Project Mapping

The best target is the fixed-`r` Route-2 certificate machinery:

- The theorem is already decomposed into finite exact checks, hub-off mode
  margins, hub-off reserve, hub-on mode domination, and hub-on perturbation.
- The proof obligations are algebraic and certificate-shaped.
- The route does not require a new global tree induction.
- The current notes already identify proof-writing gaps:
  - adjacent `F` margins imply global `F` margin;
  - hub-on perturbation changes the tie fugacity and bridge mean by a bounded
    amount;
  - positive shifted-coefficient certificates need a theorem-level statement.

Second-best targets:

- The tie-fugacity bridge in `gpt_attack/SG3_ROUTE2_PACKET.md`.
- Restricted STP2 closure in `Formal/STP2Closure.lean`, but only with stronger
  tree-realizable or gap-free hypotheses; the abstract closure is known false
  under the current weak hypotheses.

Avoid:

- Asking for a direct proof of the full tree unimodality conjecture.
- Asking for generic log-concavity machinery.  Tree log-concavity is false.
- Reopening obsolete SCC or broad compression arguments.

## Concrete Follow-Up

Create an AxiomProver-style handoff packet for the fixed-`r` certificate lemma.
It should include:

- a concise README;
- a natural-language theorem target;
- a task prompt for a proof agent;
- a Lean-facing `problem.lean` with small formal lemmas and the intended
  certificate bridge statements.

This is now added at:

```text
gpt_attack/axiom_fixed_r_certificate/
```

The first pass also proves the adjacent-margin and fugacity-perturbation
arithmetic in Lean; `problem.lean` checks with no `sorry`s.  The remaining
formalization target is the finite-support Gibbs mean-shift lemma and the full
certificate criterion.

## Parking Decision

After creating the packet, the main project should be mothballed again.
Reopen only if:

- E-JC reviews arrive;
- AxiomProver or another strong Lean agent is available to run against
  `gpt_attack/axiom_fixed_r_certificate/`;
- a deliberately bounded proof session targets only the finite-support Gibbs
  mean-shift lemma.

Do not reopen broad ECMS, Conjecture A, STP2 closure, or large-compute threads
without a new reason.

## Verification Note

I ran `lake build` during the scan.  It replayed through
`Formal/STP2Closure.lean`, confirmed the two expected `sorry` declarations
there, and then hung during `Formal/P3.lean` replay.  I terminated the build
process.  No result verification should be inferred from this run beyond the
observed warning set.
