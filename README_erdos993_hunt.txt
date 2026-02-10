# Erdos #993 Hunt (Parallel Mixed Search)

## Status snapshot (as of 2026-02-09)
- Erdos Problem #993 is marked open on the Erdos Problems site; the problem page shows a last edit on 01 Feb 2026.
- The discussion thread reports exhaustive verification for trees up to n <= 29 (comment dated 07 Jan 2026).
- Recent work:
  - "Independent set sequence of some linear hypertrees" (arXiv:2409.15555; submitted 23 Sep 2024, last revised 30 Jan 2026) notes the tree unimodality question remains open.
  - "Trees with non log-concave independent set sequences" (arXiv:2502.10654; submitted 15 Feb 2025, last revised 23 Jan 2026) constructs non-log-concave tree sequences; this script biases toward those families.

Sources:
- https://www.erdosproblems.com/993
- https://www.erdosproblems.com/forum/thread/993
- https://arxiv.org/abs/2409.15555
- https://arxiv.org/abs/2502.10654

## What this script does
- Builds (or loads) a cached library of ~2000 trees (random + structured + "wiggly" families).
- Runs two tracks in parallel:
  1) Forest search: pick k components (default k=2..5) and convolve their polynomials.
  2) Star-of-subtrees search: attach a new root to k library roots (default k=3..10).
- Also includes a double-star search (two hubs connected by an edge, each with its own branches).
- Every ~15 minutes, worker 0 runs a guided beam search in polynomial space over the top of the library.
  Disable with: --guided-every 0
  The guided pass ranks candidates using log-concavity breaks and total-variation "wiggliness."

## Run
- python erdos993_hunt.py
  Mixed mode, runs until found or stopped. Default workers = CPU count minus one.
- python erdos993_hunt.py --time-limit 36000
  Stop after ~10 hours.
- python erdos993_hunt.py --workers 2 --library-size 800 --max-n 180
  Lower resource usage.
- python erdos993_hunt.py --mode double --double-left-min 2 --double-left-max 6 --double-right-min 2 --double-right-max 6
  Search only the double-star construction.
- python erdos993_hunt.py --galvin-sweep --sweep-max-n 260 --sweep-d-max 80 --sweep-e-max 30 --sweep-t-max 10 --sweep-top-k 200
  Deterministic Galvin-style sweep; writes `out_erdos993/galvin_sweep.json` and exits.
- python erdos993_hunt.py --sweep-report out_erdos993/galvin_sweep.json --sweep-include-top 100
  Seeds the main library with the top sweep entries before the mixed search.
- python erdos993_hunt.py --sweep-report out_erdos993/galvin_sweep.json --sweep-include-top 200 --sweep-only
  Run the main search using only sweep-ranked candidates (no random library).
- python erdos993_hunt.py --sweep-report out_erdos993/galvin_sweep.json --sweep-include-top 200 --sweep-only --guided-once
  Run a single guided beam-search pass and exit (fast, deterministic).

Requirements:
- Python 3.9+
- Optional: numpy for faster polynomial convolution.

## Outputs
- out_erdos993/library.pkl
  Cache of library trees.
- out_erdos993/progress_worker*.json
  Progress checkpoints (every ~5 minutes by default).
- out_erdos993/galvin_sweep.json
  Report of top candidates from the deterministic Galvin-style sweep.
- Certificates:
  - out_erdos993/certificate_forest_YYYYMMDD_HHMMSS.json
  - out_erdos993/certificate_tree_star_YYYYMMDD_HHMMSS.json
  - out_erdos993/certificate_tree_YYYYMMDD_HHMMSS.json (if found during library build)
  - out_erdos993/certificate_tree_sweep_YYYYMMDD_HHMMSS.json (if found during Galvin sweep)

## If a certificate is found
1) Confirm forest_unimodal or unimodal is false.
2) Post the JSON and the exact command/parameters to the Erdos Problems #993 forum thread:
   https://www.erdosproblems.com/forum/thread/993
