# Orchestrator v13 Artifact Run (n<=25) — 2026-03-03

This run uses the new deterministic `orchestrator_v13` and `pi_n` runtimes committed in `5c4f8ac` + `b1fa4af`.

## Inputs

- Input file: `results/orchestrator_v13_input_from_modal_n24_n25.json`
- Source material:
  - `results/alpha_bookkeeping_modal_n24_n24_w256.json`
  - `results/lambda_frontier_modal_n24_n24_w256.json`
  - `results/n25_modal_frontier_authoritative_2026-03-03.json`

### Normalization note

The source modal frontier files do not natively contain full `WitnessRaw` rows. The input file is a deterministic heuristic normalization used for executable pipeline testing:

- `C10 := R`, `C01 := 0`, `C11 := 0`
- `sum_err := sum_all`
- `Lambda := X + D` (or `Lambda := X` when `D` missing)
- `m := k + 1`
- `Gamma_L := 0`

This makes the run reproducible but should be treated as a runtime integration test artifact, not a theorem-level calibrated dataset.

## Commands

```bash
PYTHONPATH=. python3 -m orchestrator_v13.cli \
  --input results/orchestrator_v13_input_from_modal_n24_n25.json \
  --out-dir results/orchestrator_v13_run_n24_n25 \
  --verify
```

Output: `verify: OK`

## Produced artifacts

- `results/orchestrator_v13_run_n24_n25/frontier_state.json`
- `results/orchestrator_v13_run_n24_n25/partition_plan.json`
- `results/orchestrator_v13_run_n24_n25/repair_queue.json`
- `results/orchestrator_v13_run_n24_n25/v13_obligations.json`

Observed selection/closure:

- selected partition: `PI2`
- global closure status: `OPEN`
- reason: A2 bucket deficit (`lambda_hat > lambda_max`), `Phi > 0`

## Pi(n) certificate from same normalized dataset

Input rows were cast to integer `InstanceRecord` fields and evaluated with `pi_n.compute_pi`.

- output: `results/pi_n_cert_from_modal_n24_n25.json`
- status: `FAIL_CLASS`
- replay verification: passed (`verify_fail_class`)
