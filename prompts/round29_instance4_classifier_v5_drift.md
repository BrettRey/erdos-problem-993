# Round 29 (Instance 4): Classifier v5 with Drift-Aware Proof Margin

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Existing classifier scores:
- `Gamma_G` (global reserve hardness)
- `Gamma_L` (local neighbour-compensation stress)
- `Gamma_M` (alpha proof-margin slack)

New drift signal from `n<=21` alpha scan:
- `alpha_min` decreased from `0.2437206585` (n<=19) to
  `alpha_* = 0.21034113597068071` (n<=21).
- odd-only alpha remains strongly negative.

## Task

Upgrade classifier to explicitly track drift risk.

Requirements:
1. Define one drift score (e.g., `Gamma_D`) that measures alpha erosion by class/size.
2. Couple `Gamma_D` with `Gamma_G`, `Gamma_L`, `Gamma_M` into one risk index for next-jump triage.
3. Produce a refined ranking for next checks that separates:
   - likely global-lambda jump classes,
   - likely alpha-route break classes.
4. Give explicit danger thresholds in terms of `alpha - lambda` gap and margin scores.

## Output format

1. `Drift-aware score definitions`
2. `Coupled risk interpretation`
3. `Refined next-jump / next-break ranking`
4. `Danger thresholds`
5. `Minimal validation pipeline`
