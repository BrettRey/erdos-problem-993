# Round 30 (Instance 1): Drift Envelope Theorem (alpha vs lambda)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Lower-face calibration has drifted to:
  - `alpha_19 = 0.2437206585182262`
  - `alpha_* (n<=21) = 0.21034113597068071`
- Previously locked upper-face constant:
  - `lambda19 = 0.08144365672607116`
- Current gap at n<=21 alpha level:
  - `alpha_* - lambda19 = 0.12889747924460955`
- Lower face and upper face both use `sum_all` (not odd-only).

## Task

Write a theorem-level upgrade from fixed constants to a drift envelope:

1. State closure with `alpha(n)` and `lambda(n)` as functions (or bounds) of size cutoff.
2. Give a one-line criterion for when closure remains valid (`alpha(n) > lambda(n)`).
3. Derive an explicit break condition and a quantitative warning threshold based on drift speed.
4. Instantiate with current locked values (`alpha_*`, `lambda19`) and report certified margin.
5. Make the theorem robust to future updates where only one face drifts first.

## Output format

1. `Drift-envelope theorem statement`
2. `Closure criterion and quantitative margin`
3. `Break condition and warning thresholds`
4. `Current instantiation (n<=21 alpha, n<=19 lambda)`
5. `Single unresolved lemma interface`

