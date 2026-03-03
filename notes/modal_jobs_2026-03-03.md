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
