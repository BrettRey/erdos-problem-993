# Modal Jobs Launched (2026-03-05)

Launched from repo root with profile `brettrey`.

## Summary

- Performed doc cleanup for the committed `n=27` LC + near-miss artifact:
  - `results/analysis_n27_modal_lc_nm.json`
- Raised Modal worker timeouts from `2h` to `12h` for the exhaustive and LC/NM
  workers so `n>=28` partitions are not artificially capped.
- Added a local collector:
  - `scripts/collect_modal_results.py`
  - this merges dict-backed shard outputs into the usual JSON artifacts without
    rerunning the jobs.

## Scripts touched for reliability

- `search_modal_exhaustive.py`
- `search_modal_exhaustive_n29.py`
- `analyze_modal_lc_nm.py`
- `analyze_modal_lc_nm_n28.py`

Changes:

- worker timeout `7200 -> 43200`
- fixed stale default `n` values in the `n=28` / `n=29` wrappers
- added `launch_partitions` to the `n=28` / `n=29` wrappers to match the base scripts

## Runs started

### 1. n=28 exhaustive unimodality

- dict: `erdos-993-n28-unimodality-results`
- expected tree count: `2,023,443,032`

Launch paths used in-session:

```bash
modal run --detach search_modal_exhaustive.py::main \
  --n 28 --workers 1024 \
  --expected 2023443032 \
  --out-json results/analysis_n28_modal_unimodality.json

modal run search_modal_exhaustive.py::dispatch --n 28 --workers 1024
```

### 2. n=28 LC + near-miss

- dict: `erdos-993-n28-lc-nm-results`

Launch paths used in-session:

```bash
modal run --detach analyze_modal_lc_nm_n28.py::main \
  --n 28 --workers 1024 \
  --top-k 200 --lc-top-k 200 \
  --out-json results/analysis_n28_modal_lc_nm.json

modal run analyze_modal_lc_nm_n28.py::dispatch \
  --n 28 --workers 1024 \
  --top-k 200 --lc-top-k 200
```

### 3. n=29 exhaustive unimodality

- dict: `erdos-993-n29-unimodality-results`
- expected tree count: `5,469,566,585`

Launch paths used in-session:

```bash
modal run --detach search_modal_exhaustive_n29.py::main \
  --n 29 --workers 1024 \
  --expected 5469566585 \
  --out-json results/analysis_n29_modal_unimodality.json

modal run search_modal_exhaustive_n29.py::dispatch --n 29 --workers 1024
```

## Live snapshot

Snapshot collected with:

```bash
python3 scripts/collect_modal_results.py status \
  --kind unimodality --n 28 --workers 1024 --expected 2023443032

python3 scripts/collect_modal_results.py status \
  --kind lc_nm --n 28 --workers 1024 --top-k 200 --lc-top-k 200

python3 scripts/collect_modal_results.py status \
  --kind unimodality --n 29 --workers 1024 --expected 5469566585
```

Observed snapshot:

- `n=28` unimodality:
  - completed partitions: `396 / 1024`
  - trees checked so far: `614,576,414`
  - counterexamples so far: `0`
- `n=28` LC + near-miss:
  - completed partitions: `75 / 1024`
  - trees checked so far: `97,539,110`
  - non-unimodal so far: `0`
  - LC failures so far: `2`
  - best near-miss so far: `0.8565665724120973`
- `n=29` unimodality:
  - completed partitions: `36 / 1024`
  - trees checked so far: `111,044,351`
  - counterexamples so far: `0`

These are in-progress counts only.

## Collector usage

Status-only snapshot:

```bash
python3 scripts/collect_modal_results.py status \
  --kind unimodality --n 28 --workers 1024 --expected 2023443032

python3 scripts/collect_modal_results.py status \
  --kind lc_nm --n 28 --workers 1024 --top-k 200 --lc-top-k 200

python3 scripts/collect_modal_results.py status \
  --kind unimodality --n 29 --workers 1024 --expected 5469566585
```

Merge finished dicts into local artifacts:

```bash
python3 scripts/collect_modal_results.py collect \
  --kind unimodality --n 28 --workers 1024 --expected 2023443032 \
  --out results/analysis_n28_modal_unimodality.json

python3 scripts/collect_modal_results.py collect \
  --kind lc_nm --n 28 --workers 1024 --top-k 200 --lc-top-k 200 \
  --out results/analysis_n28_modal_lc_nm.json

python3 scripts/collect_modal_results.py collect \
  --kind unimodality --n 29 --workers 1024 --expected 5469566585 \
  --out results/analysis_n29_modal_unimodality.json
```

## Recommendation logged this session

- `n=28`: yes
- `n=29` unimodality only: yes
- `n=30`: no

Reason:

- extrapolating from the committed `n=27` Modal artifacts gives rough `1024`-worker
  wall times of about:
  - `n=28` search: `3.5h`
  - `n=28` LC/NM: `7.2h`
  - `n=29` search: `9.5h`
  - `n=29` LC/NM: `19.5h`
  - `n=30` search: `25.6h`
  - `n=30` LC/NM: `52.9h`

So `n=30` is not the next sensible burn under the current setup.
