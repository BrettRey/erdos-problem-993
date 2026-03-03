# Modal Jobs Launched (2026-03-03)

Launched from repo root with profile `brettrey`.

## Active jobs

1. n=28 exhaustive unimodality
- script: `search_modal_exhaustive.py`
- command: `modal run --detach search_modal_exhaustive.py --n 28 --workers 1024`
- app id: `ap-CniND5AhuvPE22I4Rj3BlR`

2. n=27 LC + near-miss analysis
- script: `analyze_modal_lc_nm.py`
- command: `modal run --detach analyze_modal_lc_nm.py --n 27 --workers 1024 --top-k 200 --lc-top-k 200`
- app id: `ap-I9ythYbuBdJJAZxIperBod`

3. alpha-bookkeeping diagnostics (`n=3..21`)
- script: `scan_modal_alpha_bookkeeping.py`
- command: `modal run --detach scan_modal_alpha_bookkeeping.py --min-n 3 --max-n 21 --workers 512`
- app id: `ap-2EaZ4rTTZlEpjJoOJyXa3j`

## Additional launches (same session)

4. alpha-bookkeeping diagnostics (`n=3..22`)
- script: `scan_modal_alpha_bookkeeping.py`
- command: `modal run scan_modal_alpha_bookkeeping.py::dispatch --min-n 3 --max-n 22 --workers 512`
- app id: `ap-UTckD3DCnlRW8b4enlD9JP`

5. n=29 exhaustive unimodality (dedicated app-name wrapper)
- script: `search_modal_exhaustive_n29.py`
- command: `modal run --detach search_modal_exhaustive_n29.py::main --n 29 --workers 1024`
- app id: `ap-McOXFdv4JRAHV0awz0nR6I`

6. n=28 LC + near-miss (dedicated app-name wrapper)
- script: `analyze_modal_lc_nm_n28.py`
- command: `modal run --detach analyze_modal_lc_nm_n28.py::main --n 28 --workers 1024 --top-k 200 --lc-top-k 200`
- app id: `ap-HBh4OE85E9jViar6JA8Sye`

7. alpha-bookkeeping drift extension (`n=22` only)
- script: `scan_modal_alpha_bookkeeping.py`
- command: `modal run scan_modal_alpha_bookkeeping.py::dispatch --min-n 22 --max-n 22 --workers 1024`
- app id: `ap-vR5lPbMsEPxaAYZdiQ76Bm`
- dict: `erdos-993-alpha-n22-n22`
- immediate snapshot after dispatch: `1 / 1024` partition rows present in dict (key `22/0/1024`)
- first shard witness snapshot (`22/0/1024`):
  - `min_alpha_all = 0.1875868239504603`
  - witness `(a,b)=(2,18)`, `step=2`, `k=9`
  - `min_alpha_odd = -2.3333333333333335`

8. lambda-frontier drift scan (`n=22` only)
- script: `scan_modal_lambda_frontier.py`
- command: `modal run scan_modal_lambda_frontier.py::dispatch --min-n 22 --max-n 22 --workers 1024`
- app id: `ap-sGXfBUv6k5ng3r7pBqelx9`
- dict: `erdos-993-lambda-frontier-n22-n22`
- immediate snapshot after dispatch: `20 / 1024` partition rows present in dict
- first shard witness snapshot (`22/0/1024`):
  - shard `lambda_max = 0.1968360500404575`
  - witness pair `(a,b)=(4,16)`, `step=2`, `k=5`
  - witness values: `need=3584476`, `R_shift=18210465`, `sum_all=23143747`, `D=19559271`

## Consolidated n=22 snapshot (later run, mod=256 complete)

- alpha dict key coverage:
  - `22/*/256`: 256 keys (complete)
  - `22/*/1024`: 1 key (partial legacy dispatch)
- lambda dict key coverage:
  - `22/*/256`: 256 keys (complete)
  - `22/*/1024`: 28 keys (partial legacy dispatch)
- `min_alpha_all` over `22/*/256`: `0.1875868239504603` (key `22/0/256`)
- `max_lambda_needed` over `22/*/256`: `0.1968360500404575` (key `22/0/256`)
- derived snapshot gap: `alpha_min - lambda_max = -0.0092492260899972`
- artifact saved:
  - `results/modal_n22_frontier_snapshot_2026-03-03.json`

## Completed full mod=256 output artifacts

- `results/alpha_bookkeeping_modal_n22_n22_w256.json`
  - `xneg_total = 11889449`
  - `min_alpha_all = 0.1875868239504603`
  - `min_alpha_odd = -2.3333333333333335`
- `results/lambda_frontier_modal_n22_n22_w256.json`
  - `xneg_total = 11889449`
  - `xneg_step2_total = 11652945`
  - `lambda_max = 0.1968360500404575`
  - `impossible_lambda = false`

## n=23 dispatch attempts (status)

- attempted alpha dispatch:
  - `modal run scan_modal_alpha_bookkeeping.py::dispatch --min-n 23 --max-n 23 --workers 256`
  - app id: `ap-jBufELMOxb9pxBXZdfJnri`
- attempted lambda dispatch:
  - `modal run scan_modal_lambda_frontier.py::dispatch --min-n 23 --max-n 23 --workers 256`
  - app id: `ap-z3cAxKbsHJeX2qPyoMQ1qO`
- current status at logging time:
  - both apps are `stopped` with no active tasks
  - no named dict rows materialized yet for `n=23`

## Launcher reliability patch + n=23 active run

To avoid local-dispatch fragility, both scanners were patched with a server-side
`launch_partitions` function:

- `scan_modal_alpha_bookkeeping.py::launch_partitions`
- `scan_modal_lambda_frontier.py::launch_partitions`

Then full `main` scans were started directly for `n=23, workers=256`:

- alpha app: `ap-TKd3pIfia7KZyf66oZBeyJ`
- lambda app: `ap-KoS2UgijCoIjt72pItuYyc`

In-progress snapshot (while running):

- alpha dict `erdos-993-alpha-n23-n23`: `45/256` shard keys present
- lambda dict `erdos-993-lambda-frontier-n23-n23`: `46/256` shard keys present

## Monitoring

- list apps: `modal app list`
- stream logs: `modal app logs <APP_ID>`
- stop app: `modal app stop <APP_ID>`
- inspect dict rows: `modal dict items <DICT_NAME>`
- count completed partitions (example):
  `modal dict items erdos-993-alpha-n22-n22 | rg -o "22/[0-9]+/1024" | wc -l`

## Notes

- All three scripts are single-local-entrypoint (`main`) Modal launchers.
- Detached local-entrypoint runs keep worker tasks active in the app; use app IDs above for monitoring.

## Consolidated n=23 snapshot (mod=256 complete)

- alpha dict key coverage:
  - `23/*/256`: 256 keys (complete)
- lambda dict key coverage:
  - `23/*/256`: 256 keys (complete)
- `min_alpha_all` over `23/*/256`: `0.18243252946998884`
  - witness: `(a,b)=(2,19)`, `step=2`, `k=9`, `root=0`
- `max_lambda_needed` over `23/*/256`: `0.24039868912946966`
  - witness: `(a,b)=(2,18)`, `step=2`, `k=5`, `root=20`
- derived snapshot gap: `alpha_min - lambda_max = -0.05796615965948082`
- drift vs n=22:
  - `delta_alpha = -0.005154294480471472`
  - `delta_lambda = +0.043562639089012145`
  - `delta_gap = -0.04871693356948362`
- artifact saved:
  - `results/modal_n23_frontier_snapshot_2026-03-03.json`

## Completed full mod=256 output artifacts (n=23)

- `results/alpha_bookkeeping_modal_n23_n23_w256.json`
  - `xneg_total = 36238085`
  - `min_alpha_all = 0.18243252946998884`
  - `min_alpha_odd = -2.5789473684210527`
  - `impossible_all = false`
- `results/lambda_frontier_modal_n23_n23_w256.json`
  - `xneg_total = 36238085`
  - `xneg_step2_total = 35474970`
  - `lambda_max = 0.24039868912946966`
  - `impossible_lambda = false`

## n=23 geometry note for prompt design

- At n=23, both global frontiers sit in the small-first line (`a=2`):
  - alpha witness class: `(2,19)`
  - lambda witness class: `(2,18)`
- This means the Round 31 two-bucket split (`a<=3` vs `a>=4`) is no longer sufficient by itself; the break is now internal to the small-first bucket.

## n=24 launch attempt (diagnostic)

Attempted to advance alpha/lambda frontier to `n=24` using the same scripts and workers=`256`.

Launches executed:

- detached main (alpha): `ap-dLkFkOpDXAmcyXziqml0Ur`
- detached main (lambda): `ap-tFjNqn32pRcrwf1uN9nJZQ`
- server-launch attempts via `launch_partitions` / `dispatch`:
  - alpha apps: `ap-6uhg4eD4FqA8W6NlO6mw3E`, `ap-IIU9IvwK73Pa4PCW6zL4VP`, `ap-8n0c3fPCedzUty12SUlXEx`
  - lambda apps: `ap-igHQc7QELf2DGX2wJ79ruW`, `ap-3AAiONRn3BwcZsSloZZ1Xg`, `ap-3UBDhYNRdU3C9vyRpWEuNu`, `ap-Mbriv0j2yFmk7d8iOIQgNw`

Observed behavior:

- shard materialization for both dicts stalled at:
  - `erdos-993-alpha-n24-n24`: `20/256`
  - `erdos-993-lambda-frontier-n24-n24`: `20/256`
- log samples from active apps showed repeated heartbeat failures and worker drops, e.g.:
  - `Modal Worker Heartbeat attempt failed ... ConnectionError('Deadline exceeded')`
  - `Runner failed with exception: Worker disappeared, in-progress inputs will be re-scheduled.`

Current status:

- n=24 data is partial and not suitable for frontier claims.
- most n=24 attempt apps have been stopped/are stopping to avoid churn.
- rerun should use a more stable launch context before treating n=24 as authoritative.

### n=24 partial frontier read (non-authoritative)

Using full dict dumps (`modal dict items --json`) after stopping stalled apps:

- alpha dict rows present: `127/256`
- lambda dict rows present: `119/256`

Provisional extrema over currently materialized rows only:

- provisional `min_alpha_all = 0.16161242603550297`
  - shard key `24/0/256`
  - witness class `(a,b)=(2,20)`, step 2, `k=10`
- provisional `max_lambda_needed = 0.280781720999777`
  - shard key `24/0/256`
  - witness class `(a,b)=(4,18)`, step 2, `k=5`

These are incomplete-frontier diagnostics only and must not be used as final n=24 claims until `256/256` rows are present in both dicts.

## n=24 completion (after missing-shard relaunch)

Final dict coverage:

- alpha dict `erdos-993-alpha-n24-n24`: `256/256`
- lambda dict `erdos-993-lambda-frontier-n24-n24`: `256/256`

Final frontiers from complete mod=256 rows:

- `alpha_front(24) = 0.16161242603550297`
  - witness class `(a,b)=(2,20)`, step 2, `k=10`
- `lambda_front(24) = 0.280781720999777`
  - witness class `(a,b)=(4,18)`, step 2, `k=5`
- derived gap:
  - `alpha_front(24) - lambda_front(24) = -0.11916929496427403`

Drift vs n=23:

- `delta_alpha = -0.02082010343448587`
- `delta_lambda = +0.040383031870307355`
- `delta_gap = -0.061203135304793214`

Artifacts written:

- `results/alpha_bookkeeping_modal_n24_n24_w256.json`
- `results/lambda_frontier_modal_n24_n24_w256.json`
- `results/modal_n24_frontier_snapshot_2026-03-03.json`

Launcher note:

- Added `launch_missing_partitions` + `dispatch_missing` to both scanner scripts to recover incomplete dicts by spawning only absent shard keys.

## n=25 launch (in progress)

Started full frontier scans at `workers=256` using missing-shard capable dispatch:

- alpha dispatch app: `ap-wJsDYjc7gGLb4nFtjnT7kf`
  - command: `modal run --detach scan_modal_alpha_bookkeeping.py::dispatch_missing --min-n 25 --max-n 25 --workers 256`
- lambda dispatch app: `ap-S908csHXcJMrnlGjUNUI3L`
  - command: `modal run --detach scan_modal_lambda_frontier.py::dispatch_missing --min-n 25 --max-n 25 --workers 256`

Initial dict materialization snapshot right after launch:

- `erdos-993-alpha-n25-n25`: `1/256`
- `erdos-993-lambda-frontier-n25-n25`: `1/256`

Use `dispatch_missing` repeatedly as needed to backfill any stalled residues until both dicts reach `256/256`.
