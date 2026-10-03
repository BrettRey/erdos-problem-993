# Re-verification after the paper began citing the full Lean proof of Theorem 1.1
<!-- SUMMARY: pass re-checks after the abstract, outline, supplement and acks were updated for the Theorem 1.1 formalization · status: diff-scoped, no new findings · updated: 2026-10-03 -->

The diff covers the abstract's last sentence, one outline sentence, the
supplement description's paragraph on the formalizations, one word in the
acknowledgement ("both" became "the"), and one keyword ("formal verification").
Mathematics and numbers are unchanged.

| Pass | Check | Result |
|---|---|---|
| negative-claims-audit, source-reread | Every claim about the formalization matches `formalization/pb_deduction_aristotle_result/LOCAL_REPLAY.md`. That covers Theorem 1.1 for $0<p_i<1$ from the pgf, the HJ cubics, strict log-concavity, the discrete proof of the maximal-mass bound, and standard axioms only. "The other results of Section 1 are not formalized" is true: only (1.6) from (1.3) is covered, by the conditional project, as the text says. | Supported. |
| contribution-alignment, adversarial-cold-read | The abstract's method sentence now says the main inequality and its cited inputs are verified in Lean. That strengthens the contribution statement, and the problem and theorem sentences are unchanged. Portal abstract re-synced: exact match. | Round-9 verdict stands. |
| terminological-hygiene, level-category-audit | "Lean proof assistant" is glossed. "Formalized" and "verified" are used with Lean as the agent and the proof as the object, which is correctly attributed. | No issue. |
| build-integrity, proofread | 0 undefined references, 0 overfull boxes, 13 pages. Read the rendered supplement and abstract. | Clean. |
| others (reader, coherence, de-ai, scar tissue, prose passes) | Read each changed sentence in context. | No issue. |
