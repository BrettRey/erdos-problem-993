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

## Monitoring

- list apps: `modal app list`
- stream logs: `modal app logs <APP_ID>`
- stop app: `modal app stop <APP_ID>`

## Notes

- All three scripts are single-local-entrypoint (`main`) Modal launchers.
- Detached local-entrypoint runs keep worker tasks active in the app; use app IDs above for monitoring.
