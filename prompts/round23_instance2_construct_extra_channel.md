# Round 23 (Instance 2): Construct the Minimal Extra Channel for Defect

Use only your prior outputs plus this prompt. No new scan claims.

## Locked context

- Fixed-family reserve `R_shift = C10+C01+C11` at `lambda0` is insufficient at `n=19`.
- Defect framework is in place:
  `sum_err <= D + lambda0*R_shift + Defect`.

Goal now is to absorb `Defect` with one explicit extra channel.

## Task

Build a minimal additional channel family from shifted `Lambda_old` moments, e.g.

`C_ab(k) = sum_i Lambda_old(i) P(k-i-a) Q(k-i-b)`

and find a parsimonious correction term

`Extra(k) = sum_j theta_j C_{a_j,b_j}(k)` with `theta_j >= 0`

such that a plausible theorem is

`Defect <= c * Extra(k)`

with small `c`.

Requirements:
1. Keep number of new shifts/channels minimal.
2. Explain why chosen shifts target the observed `n=19` failures.
3. Provide one explicit candidate `(Extra, c)` and one fallback candidate.
4. Show resulting upgraded global inequality explicitly.

## Output format

1. `Primary Extra channel formula`
2. `Fallback Extra channel formula`
3. `Resulting global inequalities`
4. `Single missing lemma`
5. `Targeted finite test plan`
