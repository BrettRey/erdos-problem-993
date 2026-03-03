# Round 34 (Instance 1): Regime-Adaptive Bucket Theorem After n=24 Cross-Line Shift

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- n=23 frontiers:
  - `alpha_front(23)=0.18243252946998884` at class `(2,19)`
  - `lambda_front(23)=0.24039868912946966` at class `(2,18)`
  - gap `=-0.05796615965948082`
- n=24 frontiers (mod=256 complete):
  - `alpha_front(24)=0.16161242603550297` at class `(2,20)`
  - `lambda_front(24)=0.280781720999777` at class `(4,18)`
  - gap `=-0.11916929496427403`
- Existing theorem stack:
  - Round 33 three-bucket split (`S18`, `Srest`, `L`)
  - bucket gaps `g_B = alpha_B - lambda_B`
  - witness obstruction format: same-bucket pair with `r_upper > r_lower` forces `g_B<0`.

## Task

Build a final regime-adaptive theorem statement that handles both geometries:

1. internal same-line conflict (n=23 style: `(2,18)` vs `(2,19)`), and
2. cross-line conflict (n=24 style: line `a=2` vs line `a=4`).

Specifically:

1. State a partition-selection rule `Pi(n)` from witness diagnostics (not ad hoc per round).
2. State the bucketwise closure theorem under `Pi(n)` with exact ratio notation.
3. Give a minimality argument for why `Pi(23)` must separate `(2,18)` from `(2,19)`, and what `Pi(24)` must at least separate given the locked witnesses.
4. Give a coarsening obstruction lemma (witness-pair style) that certifies when a proposed partition cannot close.
5. Provide a finite computation/falsification checklist that outputs bucket gaps and the explicit failing witness pair in each bad bucket.

## Output format

1. `Regime-adaptive repaired theorem (final form)`
2. `Minimal partition requirements at n=23 and n=24`
3. `Coarsening obstruction lemma`
4. `Finite partition-evaluation protocol`
5. `Failure localization checklist`
