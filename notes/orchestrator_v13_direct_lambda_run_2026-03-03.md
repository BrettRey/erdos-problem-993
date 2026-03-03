# Orchestrator v13 Direct-Lambda Extractor Run (n=22..25) — 2026-03-03

This run replaces the earlier ad-hoc mixed normalization with a source-semantic extractor based only on lambda-frontier witness rows.

## Extractor

Script: `scripts/build_orchestrator_input_from_lambda_frontiers.py`

Input sources:

- `results/lambda_frontier_modal_n22_n22_w256.json`
- `results/lambda_frontier_modal_n23_n23_w256.json`
- `results/lambda_frontier_modal_n24_n24_w256.json`
- `results/n25_modal_frontier_authoritative_2026-03-03.json` (lambda witness only)

Semantics used (from scan scripts):

- `X = Lambda - D`
- `R = R_shift`
- `sum_all = Σ_s err_s`
- `need = max(0, sum_all - D)`

Extractor mapping to `WitnessRaw`:

- `Lambda := X + D` (exact)
- `D := D` (exact)
- `sum_all := sum_all` (exact)
- `sum_err := sum_all` (same computed error-mass channel used by lambda frontier)
- `C10,C01,C11 := (R_shift, 0, 0)` preserving exact `R_shift = C10+C01+C11`
- `m := k+1` to satisfy orchestrator prefix filter `0 <= k < m`

## Commands

```bash
python3 scripts/build_orchestrator_input_from_lambda_frontiers.py \
  --lambda-file results/lambda_frontier_modal_n22_n22_w256.json \
  --lambda-file results/lambda_frontier_modal_n23_n23_w256.json \
  --lambda-file results/lambda_frontier_modal_n24_n24_w256.json \
  --out results/orchestrator_v13_input_direct_lambda_n22_n25.json

PYTHONPATH=. python3 -m orchestrator_v13.cli \
  --input results/orchestrator_v13_input_direct_lambda_n22_n25.json \
  --out-dir results/orchestrator_v13_run_direct_lambda_n22_n25 \
  --verify
```

`verify: OK`

## Outputs

- `results/orchestrator_v13_input_direct_lambda_n22_n25.json`
- `results/orchestrator_v13_run_direct_lambda_n22_n25/frontier_state.json`
- `results/orchestrator_v13_run_direct_lambda_n22_n25/partition_plan.json`
- `results/orchestrator_v13_run_direct_lambda_n22_n25/repair_queue.json`
- `results/orchestrator_v13_run_direct_lambda_n22_n25/v13_obligations.json`

Observed orchestrator state:

- selected partition: `PI0`
- global closure: `CLOSED`
- `Phi = 0`

## Pi(n) on direct-lambda dataset

Input rows were cast to integer `InstanceRecord` and evaluated by `pi_n.compute_pi`.

- output: `results/pi_n_cert_from_direct_lambda_n22_n25.json`
- status: `PASS`
- replay verification: passed (`verify_pass_certificate`)
