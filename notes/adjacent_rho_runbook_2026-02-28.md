# Adjacent `(m,lambda,rho)` Runbook (2026-02-28)

## Live run
- Process:
  - `python scripts/adjacent_rho_split_scan_minu.py --min-n 3 --max-n 27 --m-min 4 --progress-every 0 --out results/adjacent_rho_split_scan_minu_mge4_n27_exact.json`
- Artifact target:
  - `results/adjacent_rho_split_scan_minu_mge4_n27_exact.json`

## Immediate post-finish commands
1. Validate witness constraints (`|delta_N|=1`, exact key equality):
   - `python3 scripts/report_adjacent_rho_result.py --in results/adjacent_rho_split_scan_minu_mge4_n27_exact.json`
2. Summarize all key artifacts in one pass:
   - `python3 scripts/summarize_mlambda_rho_outputs.py results/canonical_projection_battery_minu_mlambda_rho_sigma_mge4_n26_exact.json results/c1c2_coverage_mlambda_rho_mge4_n24_exact.json results/adjacent_rho_split_scan_minu_mge4_n27_exact.json`

## Interpretation matrix
- If `adjacent_split_found=false`:
  - No adjacent split through `n<=27` for projection `(m,lambda,rho)`.
  - Combine with existing `n<=26` no-split battery to report:
    - no `(m,lambda,rho)` split witness through `n<=27`.
- If `adjacent_split_found=true`:
  - Record witness `g6` pair, `N_A,N_B`, `|delta_N|`.
  - Confirm exact equality of `m`, `lambda`, `rho` and include rebuilt records.

## Current theory position (A/B consolidation)
- Proven from allowed inequalities:
  - adjacent overlap is impossible for **even** `N` in pair `(N,N+1)`.
- Remaining unresolved window:
  - odd->even pairs `(2k-1,2k)`.
- Sharpened blocker:
  - any odd->even overlap must satisfy the extremal even remainder regime
    (`deg P_even = N_even/2`), i.e. balanced bipartite remainder.
- Existing c1/c2 coverage artifact:
  - `results/c1c2_coverage_mlambda_rho_mge4_n24_exact.json` shows zero coverage
    (`c1_true=0`, `c2_true=0`, `c1c2_true=0`), so that criterion did not activate.
