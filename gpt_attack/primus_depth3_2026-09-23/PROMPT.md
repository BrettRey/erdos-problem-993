# Assignment: the residual depth-three window for tree independence sequences

Work on one precise open target, with the attached package as your starting
point. The motivation is Erdős Problem #993, but this assignment is not a
request to claim that the full problem follows from a single remaining lemma.

For every finite simple tree T with n vertices, independence number alpha,
33 <= n <= 38, alpha in {17,18,19}, and 2 alpha - n <= 5, prove or refute

    s_3^2 >= s_2 s_4,

where s_d is the number of independent sets of size alpha-d.

Read `MATHEMATICS.md` before starting. It gives intrinsic extendable and
blocked counts e_d and b_d, with s_d=e_d+b_d, and the exact remaining regime.
Established subcases are supplied in Lean source: b_1=0 and low b_4/e_4 density.
The remaining branch has b_1>0, high density, and a negative combined correction
D. The Pascal reserve and blocked-shadow bound are already available. Do not
spend the run re-proving those subcases unless you find a specific defect.

First run `replay.py` and the regression tests to validate the supplied
counterexamples to false auxiliary claims. Treat a proposed new inequality as
falsifiable: check its hypotheses and adversarial examples before promoting it.
Use exact integer/rational arithmetic for signs. No float64 root-finding.

You may use a different proof route. A sufficient bound that the Pascal margin
absorbs D would be useful, but it is not the primary target. Refuting that bound
does not refute the target. Refuting the target would not by itself refute
unimodality or solve #993. Do not assume blocked-profile log-concavity,
componentwise nonadversity, real-rootedness, mode–mean localization, or an
unproved Holant representation. See `FAILED_ROUTES.md`.

Return one of the following grades, supported by the corresponding artifacts:

- **PROVED:** a complete proof of the unchanged universal window target, with
  every reduction justified. Stronger strict positivity is welcome but optional.
- **REFUTED:** an exact graph6/edge-list witness in the stated window, its full
  independence polynomial, and a standalone exact verification. Explicitly
  distinguish target refutation from any auxiliary refutation.
- **PARTIAL:** at least one new proved lemma or precisely delimited structural
  reduction, with its usefulness and remaining gap stated. Do not label a larger
  finite sample as a proof.
- **NO ADVANCE:** report failed approaches and exact obstructions honestly.

The deliverable should contain a concise written argument, a status/dependency
table, scripts and certificates, exact reproduction commands and versions,
and any Lean source/build/axiom audit. Separate machine-checked proofs, written
proofs, exact finite checks, and conjectural evidence. The supplied formal target
definitions must remain unchanged; no new hypothesis may be silently added.
Do not replace a graph theorem with a scalar implication that assumes its hard
counting inequality.

Prefer one verifiable advance over a broad survey or a speculative paper. Work
within the run budget selected by the user, and return partial results rather
than starting additional paid runs, purchasing resources, or publishing work.
