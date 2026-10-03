INDEPENDENT CHECKS OF THE POISSON-BINOMIAL TURAN PAPER
Prepared: 2026-10-03

Manuscript checked:
  variance-scaled-turan-first-descent(1).pdf
  SHA-256: 0d17c6cc01e99f111599a045757605a5798c90f4beeb2b9df4920dae8b2b8781

These programs independently reconstruct selected calculations in the attached
revision. The author's certificate supplement was not attached and was not used.
The manuscript itself is not included in this archive.

REQUIREMENTS AND USE
Python 3.10 or later. Only the Python standard library is used.
Extract this archive and run:
  python3 audit_certificates.py
  python3 audit_cue.py
Each program writes JSON and a human-readable log next to the program.
The included JSON and log files record the completed independent checks.

SECTION 4 CERTIFICATES
- Rebuilds the unsymmetrized K=3 numerator in equation (4.6) from the
  defining mass bounds and pairwise variance form.
- Rebuilds the twelve symmetric cell polynomials for m=4,...,15.
- Converts all finite polynomials into Bernstein form with exact rational
  arithmetic and verifies the conversions by symbolic reconstruction.
- Checks that all 275 finite-range Bernstein coefficients are positive.
- Reconstructs the J=5 quartic coefficients and all five symbolic
  factorizations for J>=6; verifies positivity of the stated polynomial
  coefficients.
- Checks the four finite-sum identities by polynomial differences and their
  zero base cases.
- The recorded run passed all 74 named exact algebraic assertions, as well
  as the reconstruction assertions. Full coefficients are in the JSON.

CUE HALF-CIRCLE NUMERICS
- Evaluates the finite rational sums in the manuscript's exact variance
  formula, using rigorous rational upper and lower bounds for pi from
  Machin's identity and alternating arctangent series.
- Certifies V_1996 < 1 < V_1998.
- Certifies the displayed bounds at N=10000.
- Decimal intervals in the output are rounded outward.

SCOPE AND LIMITS
The written probabilistic arguments, the determinantal representation,
the variance formula, its monotonicity and asymptotic expansion were checked
analytically in the accompanying review. These scripts do not formalize those
arguments. They are independent arithmetic checks, not a machine-checked
formalization of the whole theorem. They do not inspect or validate the
author's original supplement, file manifest, repository, or Lean development.
