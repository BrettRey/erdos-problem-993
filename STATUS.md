# Erdős Problem #993 — Independent Set Sequence Unimodality for Trees

**Problem:** Is the independent set sequence of every tree unimodal?

**Source:** Alavi, Malde, Schwenk, and Erdős (1987)
**Reference:** [erdosproblems.com #993](https://erdosproblems.com/993)

## Source of truth

This status mirrors the current manuscript in `paper/main.tex`. Numeric snapshots live in `results/*.json` where available. This file is not automatically recomputed.

## Results (as reported in `paper/main.tex`)

### Exhaustive verification through n = 26

The manuscript reports no unimodality violations among all 447,672,596 non-isomorphic trees on `n <= 26` vertices. Tree counts match OEIS A000055.

| n | Trees | Time |
|---|------:|------:|
| 1--15 | 13,188 | <1s |
| 16 | 19,320 | 1s |
| 17 | 48,629 | 3s |
| 18 | 123,867 | 9s |
| 19 | 317,955 | 23s |
| 20 | 823,065 | 68s |
| 21 | 2,144,505 | 55s |
| 22 | 5,623,756 | 1m 44s |
| 23 | 14,828,074 | 4m 41s |
| 24 | 39,299,897 | 12m 5s |
| 25 | 104,636,890 | 38m 33s |
| 26 | 279,793,450 | 4h 51m |
| **Total** | **447,672,596** | |

Additional n=26 details from the manuscript:
- Exactly 2 log-concavity failures (both at k = 13).
- Best near-miss ratio `nm = 0.845` (well below the violation threshold of 1).

### Targeted search on structured families (n up to 500)

The manuscript reports 145,362 tested trees across five families, with no unimodality violations.

| Family | Trees | LC failures | Best nm |
|--------|------:|----------:|---------|
| Galvin SST T_{m,t,1} | 571 | 108 | 0.936 |
| Generalized SST T_{m,t,d} | 680 | 268 | 0.981 |
| Caterpillars | 5,196 | 0 | -- |
| Spiders and brooms | 133,915 | 0 | 0.992 |
| Random Ramos-Sun-style | 5,000 | 2 | 0.804 |
| **Total** | **145,362** | **378** | |

Patterns highlighted in the manuscript:
- Subdivided stars dominate log-concavity failures.
- Caterpillars are log-concave through `n = 500` in the tested set.
- Brooms produce the highest near-miss ratios without log-concavity failure.

### Broom asymptotics

For broom `B(13, s)` the manuscript reports:

| s | n | nm | 1 - nm | s * (1 - nm) |
|---:|---:|:---------|:-------|:-------------|
| 1,000 | 1,013 | 0.995911 | 0.00409 | 4.089 |
| 2,000 | 2,013 | 0.997947 | 0.00205 | 4.105 |
| 5,000 | 5,013 | 0.999177 | 0.00082 | 4.115 |
| 10,000 | 10,013 | 0.999588 | 0.00041 | 4.119 |
| 20,000 | 20,013 | 0.999794 | 0.00021 | 4.120 |

The data are consistent with `nm(s) = 1 - C/s + O(1/s^2)` and `C ~ 4.12`.

## Artifacts in this repo

- `results/analysis_n26.json` contains the n=26 exhaustive LC and near-miss summary.
- `results/analysis_n15.json` contains a smaller-scope exhaustive summary.
- `results/targeted_n500.json` contains the targeted search totals and top near-misses.
- `results/targeted_families.json` contains the per-family summary for the targeted search.
- `results/two_branch_lc_n24.json` contains the exhaustive C2 scan (`#{v:deg(v)>=3} <= 2`) through `n=24`: 196,635 C2 trees, 0 LC failures, worst LC ratio `0.846153846...`.
- `results/squeeze_n12_networkx.json` and `results/squeeze_n20_geng.json` contain squeeze scans (first descent vs. Levit–Mandrescu tail). Through n=20 the strict condition “first descent >= tail start” fails, with worst margin `-2`.
- `out_erdos993/*.json` contains heuristic/parallel search progress logs and sweeps.

The broom asymptotic table above is not saved to disk; reproduce it with `python broom_asymptotic.py`.

## Method (from manuscript)

- Tree enumeration via nauty/geng (`geng -c n n-1:n-1`), partitioned with `res/mod`.
- Independence polynomial DP on a rooted tree, using polynomial convolution.
- Unimodality check by forbidding any rise after the first strict descent.
- Log-concavity check using integer arithmetic (`i_k^2 >= i_{k-1} i_{k+1}`).
- Parallelism via `multiprocessing` when using geng.

## Reproduction

```bash
pip install networkx numpy
brew install nauty

python test_all.py
python search.py --max-n 26 --workers 8
python analyze.py 26 --workers 8
python targeted.py --max-n 500 --random-count 5000
python broom_asymptotic.py
```

## Open directions (from manuscript)

1. Prove an elementary broom-unimodality argument independent of the spider theorem, or extend log-concavity proofs to larger structural classes (e.g., trees with at most two branch vertices) via transfer-matrix or rooted-pair invariants.
2. Extend exhaustive search to n = 27 (about 751M trees; current Python code estimated at 12–15 hours on 8 cores).
3. Adapt PatternBoost-style search to optimize near-miss ratio rather than log-concavity failure.

## Literature notes (cited in the manuscript)

- Kadrawi & Levit (2023) found two log-concavity failures at n = 26.
- Galvin (2025) constructed subdivided-star families with log-concavity failures deep in the sequence.
- Ramos & Sun (2025) used PatternBoost to find many log-concavity failures up to n = 101 with no unimodality failures.
- Levit & Mandrescu (2006) proved a monotone tail bound for independence sequences of trees.
- Li, Li, Yang, and Zhang (2025) proved all spiders are strongly log-concave, implying broom log-concavity.
