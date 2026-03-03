# Round 33 Prompt Prep (2026-03-03)

Goal: convert Round 32 theory outputs into certifiable finite interfaces that can be checked directly from witness/channel data.

## Inputs received (Round 32, 4 outputs)

1. Internal-break theorem repair:
   - Three-way bucket split at n=23:
     - `S18`: `(a,b)=(2,18)`
     - `Srest`: `a<=3` and not `(2,18)`
     - `L`: `a>=4`
   - closure criterion via bucket gaps `(g_S18, g_Srest, g_L)`.

2. a=2 local repair family:
   - two a=2 signatures (center mismatch vs edge-band deficit)
   - correction mass `M = rho*Psi + lambda*Phi`
   - transfer from `a=2` slot to `a=3` sink, sum-preserving.

3. Induction v10:
   - hard-step `t=2` requires class-stratified channel-vector sandwich
   - bridge `t>=3` can stay weaker
   - non-circular interface retained via `Lambda = D + X_pos - sum_all`.

4. Classifier v8:
   - internal-deficit diagnostics (line-level and within-line mismatch)
   - priority moves now focused on `a=2` line at n=23
   - refined queue and stop conditions.

## Locked frontier constants (for prompt grounding)

- `alpha_front(23) = 0.18243252946998884`
- `lambda_front(23) = 0.24039868912946966`
- `gap = -0.05796615965948082`
- alpha witness class `(2,19)`, lambda witness class `(2,18)`.

## Round 33 objectives

1. Turn three-bucket theorem from structural proposal into a minimality/necessity claim with explicit falsification geometry.
2. Turn `(rho, lambda)` local repair into a concrete feasibility region (linear constraints + calibration workflow).
3. Compress v10 hard-step obligations to smallest finite class interface that still handles n=23 internal break.
4. Turn v8 diagnostics into an implementable computation protocol (exact aggregates + queue decision logic).

## Constraints

- Keep all-diagonal framework (`sum_all`), boundary-correct indexing, and prefix filter fixed.
- No odd-only reopening.
- No fresh scan claims; use locked data and clearly mark any conditional statements.
- Every prompt must demand a finite falsification/certification checklist.
