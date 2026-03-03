# Round 30 Prompt Prep (2026-03-03)

Goal: push the drift-aware framework from Round 29 into a proof-ready envelope and an operational frontier tracker.

## Locked carry-forward facts

- All-diagonal alpha calibration drift:
  - `alpha_19 = 0.2437206585182262`
  - `alpha_* (n<=21) = 0.21034113597068071`
  - erosion `Delta_alpha = 0.033379522547545476`
- Upper reserve constant currently locked through prior regime:
  - `lambda19 = 0.08144365672607116`
- Current alpha-lambda gap at drifted alpha:
  - `alpha_* - lambda19 = 0.12889747924460955`
- Odd-only lower-face route remains unusable (`min_alpha_odd < 0`).
- Lower-face mechanism candidate is packet-pair local inequality summing to
  `G_k >= alpha_* R_shift`.

## Round 30 objectives

1. Convert fixed-constant theorem to drift-envelope form `alpha(n), lambda(n)`.
2. Convert packet-pair mechanism into a finite inequality family (`a in {2,3,4}`).
3. Produce induction pack v8 with explicit two-bottleneck interface and drift threshold failure mode.
4. Upgrade classifier to frontier-tracking mode with trigger thresholds and two prioritized queues.

## Prompt design decisions

- Keep all statements in all-diagonal form (`sum_all`) to avoid odd-only dead ends.
- Keep the broken bookkeeping inequality excluded.
- Make each instance produce one clear artifact:
  - theorem envelope,
  - local bottleneck reduction,
  - induction interface,
  - operational classifier.

