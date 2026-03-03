# Round 30 (Instance 4): Classifier v6 (Frontier Tracking + Route Risk)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Current score set:
- `Gamma_G`: global reserve hardness (lambda demand),
- `Gamma_L`: local neighbour-compensation stress,
- `Gamma_M*`: margin to drifted alpha calibration,
- `Gamma_D`: alpha-erosion drift score.

Current constants:
- `alpha_19 = 0.2437206585182262`
- `alpha_*  = 0.21034113597068071`
- `lambda19 = 0.08144365672607116`
- `Delta_alpha = alpha_19 - alpha_* = 0.033379522547545476`.

## Task

Upgrade to a frontier-tracking classifier that is directly operational:

1. Define one compact frontier state vector at cutoff `n` containing:
   - current lambda frontier,
   - current alpha frontier,
   - classwise drift drivers.
2. Define one scalar risk index that combines `Gamma_G, Gamma_L, Gamma_M*, Gamma_D`.
3. Give a ranked “next checks” schedule split into:
   - global-lambda jump queue,
   - alpha-erosion/break queue.
4. Add explicit trigger rules for `warning`, `critical`, and `break`.
5. Provide a minimal data pipeline to compute this from stored witness logs.

## Output format

1. `Frontier state definition`
2. `Unified risk index`
3. `Ranked next-check queues`
4. `Trigger rules (warning/critical/break)`
5. `Minimal computation pipeline`

