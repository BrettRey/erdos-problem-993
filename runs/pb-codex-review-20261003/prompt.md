You are an independent referee and verifier from a different model family than the author's assistants. Be adversarial. Report errors, gaps and overclaims, not praise. You have a read-only shell: you may read files and run read-only commands, including `python3` for exact arithmetic with `fractions`, but do not modify anything. Give your full report on stdout as your final answer.

Repository root (your working directory): the Erdős #993 research repository. The paper under review:

- Manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex` ("Variance and local log-concavity of Poisson–binomial laws", for Electronic Communications in Probability). Main result, Theorem 1.1: V·δ_D ≥ 1/4 for Poisson–binomial laws with V ≥ 1.
- Supplement audit manifest: `paper/poisson_binomial/CERTIFICATE.md`
- Check program for §4.1: `scripts/verify_pb_compact_monotone.py`; for §4.2: `scripts/verify_pb_large_h_range.py`; for Example 1.6: `scripts/verify_pb_cue_threshold.py`
- Lean proof of Proposition 3.1: `formalization/pb_scalar_inequality_aristotle_result/` (definitions `PBScalar/Defs.lean`, statements `PBScalar/Statement.lean`, compact range `PBScalar/Compact.lean`, large range `PBScalar/Large.lean`)
- Lean proof of Theorem 1.1: `formalization/pb_deduction_aristotle_result/` (definitions `PBDeduction/Defs.lean`, statements `PBDeduction/Statement.lean`, proofs in the other `PBDeduction/*.lean` files)
- The replay records `formalization/*/LOCAL_REPLAY.md` were written by the author's assistant. Treat them as claims, not evidence.

Recent changes you should scrutinize most closely:
(a) §4.1 (label `sec:compact`) now proves the scalar inequality on 3 < H ≤ 16 by a monotonicity argument (A_K(δ) nonincreasing in δ for fixed K, the target decreasing) plus 32 exact endpoint checks, with rounded-down margins printed. §4.2 (label `sec:large`) now begins with the symmetrization.
(b) The paper now claims that Theorem 1.1, including the Hillion–Johnson cubic inequalities, strict log-concavity and the maximal-mass bound, is formalized in Lean (abstract last sentence; proof outline; `\sdescription{...}` in the supplement box; acknowledgement).

Tasks:
1. **Mathematics.** Check §4 completely (preamble, §4.1, §4.2), recomputing what you can in exact arithmetic (the 32 endpoint ratios, the printed J = 5 Bernstein coefficients, a spot check of the J ≥ 6 identity). Then check §§2–3: the recurrence (2.4), the reciprocal bound (2.5), Lemma 2.1, the mass bounds, the pairwise variance bound, the maximal-mass bound and its proof, and the deduction of Theorem 1.1 from Proposition 3.1. Report any error, gap, unjustified step, off-by-one or boundary case.
2. **Formal statement fidelity.** Read `PBDeduction/Defs.lean` and `PBDeduction/Statement.lean`. Decide whether `theorem_1_1` states exactly Theorem 1.1 of the paper, under the paper's standing assumption 0 < p_i < 1. Check the definitions of `pbPmf`, `pbVar`, `IsFirstDescent` and `deficit`, the index conventions, and the zero extension. Do the same for `PBScalar.scalar_inequality` against Proposition 3.1 and the definition of A(δ) (eq. `A-def`, with K from eq. `K-a`). Look for any way the formal statements could be weaker than the paper's: vacuous hypotheses, a definition that trivializes something, coercion or `toNat` pitfalls, or division by zero conventions. Grep the proof files for `sorry`, `admit`, `axiom`, `native_decide`, `implemented_by` and `ofReduceBool`. You cannot run `lake build`; say so.
3. **Claims.** Check every claim the manuscript, CERTIFICATE.md and the supplement description make about the computations and formalizations against the files. Flag any overclaim or anything stale.
4. **Anything else** a referee at ECP would flag as a mathematical error or a seriously misleading statement anywhere in the paper. Keep pure style comments to a short list at the end, or omit them.

Report format:
- **Verdict:** one paragraph.
- **Findings:** each with location (line number or label), severity (error / gap / imprecision / minor), a concrete description, and a proposed fix.
- **What you checked and how:** list the computations you ran, with their results.
- **What you could not check.**

For any negative claim ("no other errors", "statement is faithful"), say what you read or ran to support it.
