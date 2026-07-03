# Ramos-Sun / Li Literature Map for Erdos 993
Date: 2026-07-03

Scope: source-grounded map of the Ramos-Sun PatternBoost paper, the Li 2026 two-family paper, and the Li--Li--Yang--Zhang 2025 symmetric-function paper, with implications for the current Erdős 993 project strategy.
## Bottom Line
The Fable 5 terrain note is broadly right that log-concavity is no longer a plausible global route and that direct unimodality is the real target. The main correction is that the new counterexample literature does **not** support a simple claim that log-concavity failures are intrinsically localized. Ramos--Sun's corpus is heavily shaped by a single-index reward function, and Bautista-Ramos/Galvin-style results show that multiple or far-tail breaks are real.

For this repository, the useful lesson is not "mine LC failures for a localized invariant." It is sharper:

1. Treat LC defects as a stress-test corpus for any candidate direct-unimodality invariant.

2. Do not use global 2-Schur positivity as an intermediate invariant, because it is equivalent to log-concavity.

3. Consider a partial 2-Schur/tail-monotonicity hybrid only if it can be stated below log-concavity and checked against the vertex-deletion/subdivision recurrences.

4. Keep the existing PNP / `d_leaf <= 1` / mode-mean lane as the primary proof architecture unless a new invariant survives computational closure tests.

## Status Corrections
- The Alavi--Erdős--Malde--Schwenk forest/tree unimodality conjecture remains open as of 2026-07-03.

- Exhaustive tree unimodality verification through `n <= 29` is a local/project result and is also reported on the Erdős Problems thread. The current manuscript reports `8,691,747,673` trees through `n <= 29`.

- The statement "the known non-log-concave families were recently shown unimodal" is too broad. Li 2026 proves unimodality for the two Kadrawi--Levit families `T_{3,m,n}` and `T^*_{3,m,n}`, not for all Galvin, Ramos--Sun, or Bautista-Ramos families.

- The statement "LC violations are localized" is false as a theorem and misleading as a heuristic. Ramos--Sun observed localization within their search setup, but Galvin constructs arbitrary-distance breaks and Bautista-Ramos constructs arbitrary-many breaks.

## Ramos--Sun 2025
Source: Eric Ramos and Sunny Sun, "An AI enhanced approach to the tree unimodality conjecture," arXiv:2510.18826v2, 2025.
### What They Establish
- They use PatternBoost with Prüfer-code encodings of labelled trees.

- The scoring function targets a single log-concavity inequality at a user-chosen index: `score_i = a_{i-1} a_{i+1} - a_i^2`.

- Positive score means log-concavity fails at that index.

- They report tens of thousands of new non-log-concave trees, with vertex counts from 27 to 101.

- Their public repository contains the code and a 60-vertex output folder; the paper says this includes around 35,000 non-log-concave trees on 60 vertices.

- Their experiments most readily found breaks at `N/2`; with path-punishing local search, they also found some breaks at `N/2 - 1`.

- They report strong search biases: convergence toward paths, a tendency to shrink the independence number toward the minimum, and difficulty finding odd-`N` or further-from-middle breaks.

### What Is Not Proved
- No theorem toward tree unimodality.

- No structural classification of non-log-concave trees.

- No proof that LC breaks are localized.

- No proof that a non-unimodal tree, if it exists, must be huge.

- No proof that the PatternBoost corpus is distributionally representative of all non-log-concave trees.

### Useful Conjectural/Empirical Leads
- Ramos--Sun's Conjecture 4.2: if a tree on `N` vertices is non-log-concave at index `N/2 - beta`, then `alpha(T) >= N/2 + beta + 1`.

- Their consistent motif across found examples is excess degree-2 structure, especially degree-2 vertices adjacent to leaves.

- Their method's failure to find global-shape disruptions is weak evidence only, because the reward function did not optimize global shape.

### Local Audit of the Public 60-Vertex Corpus
I cloned `https://github.com/ericgramos/TreeUnimodalityPatternBoost` outside the repo at `/tmp/erdos993-ramos-sun` and analyzed the final top-performer file:

```bash
python3 scripts/analyze_prufer_corpus.py \
  /tmp/erdos993-ramos-sun/60_vertex_output/search_output_11.txt \
  --scores /tmp/erdos993-ramos-sun/60_vertex_output/training_distribution.txt \
  --limit 0 \
  --target-k 30 \
  --canonicalize \
  --out results/ramos_sun_60_epoch11_summary.json
```

The corpus file contains 50,000 labelled Prüfer codes of length 58, so these are 60-vertex trees. The exact recomputation of the final epoch found:

- 39,766 non-log-concave trees, all with the single LC defect at `k = 30`.
- 0 non-unimodal trees.
- 43,595 unique nauty-canonical graph6 trees among the 50,000 labelled rows.
- `alpha = 31` for 39,979 rows, `alpha = 30` for 10,019 rows, and `alpha = 32` for 2 rows.
- `d_leaf <= 1` for 49,998 rows; only 2 rows have a multi-leaf hub.
- maximum LC ratio `3.3358165094122856`; maximum near-miss ratio `0.8493940107104357`.

Interpretation: this final-epoch file is a useful adversarial test set for the `d_leaf <= 1` lane, but it is not representative of all non-LC trees. It is explicitly selected for the `k = 30` LC score under the no-path local objective.

I then ran the same exact recomputation on all eleven `search_output_i.txt` top-50,000 epoch files and aggregated the trend in `results/ramos_sun_60_epoch_trends.json`. The number of non-LC rows by epoch is:

```text
1, 282, 6222, 20033, 24163, 27281, 30345, 32708, 34953, 37688, 39766
```

Every analyzed row in every epoch remained unimodal. The maximum LC ratio increased from `1.29827176778651` in epoch 1 to `3.3358165094122856` by epoch 5 and then plateaued. The maximum near-miss ratio did **not** move toward 1 as LC failures became dense: it was `0.8764440260284707` in epoch 1 and then stabilized near `0.8493940107104357` from epoch 9 onward. This is a real empirical signal for the current project: optimizing a severe single-index LC defect is not, in this corpus, optimizing proximity to a unimodality valley.

### Local Audit of Structured Non-LC Families
I added Galvin and Bautista--Ramos generators to the same exact-analysis harness:

```bash
python3 scripts/analyze_prufer_corpus.py \
  --family galvin \
  --m-values 1-50 \
  --t-values 1-50 \
  --max-n 500 \
  --limit 0 \
  --top 10 \
  --out results/galvin_family_grid_m1-50_t1-50_nle500_summary.json

python3 scripts/analyze_prufer_corpus.py \
  --family bautista-ramos \
  --m-values 1-8 \
  --t-values 1-8 \
  --limit 0 \
  --top 10 \
  --out results/bautista_ramos_family_grid_m1-8_t1-8_summary.json

python3 scripts/analyze_prufer_corpus.py \
  --family li \
  --m-values 1-50 \
  --t-values 1-50 \
  --limit 0 \
  --top 10 \
  --out results/li_family_grid_m1-50_n1-50_summary.json

python3 scripts/analyze_prufer_corpus.py \
  --family li-star \
  --m-values 1-50 \
  --t-values 1-50 \
  --limit 0 \
  --top 10 \
  --out results/li_star_family_grid_m1-50_n1-50_summary.json
```

The Galvin grid processed 759 trees with `n <= 500`: 390 were non-log-concave, all expected single break positions were hit for those 390, and 0 were non-unimodal. The maximum near-miss ratio was `0.9808181929611438` at `(m,t) = (21,11)`.

The Bautista--Ramos grid processed 64 trees: 28 were non-log-concave, 23 had all asymptotic expected break positions realized in this small-`t` window, and 0 were non-unimodal. The maximum near-miss ratio was `0.972953314235538` at `(m,t) = (8,8)`, where the polynomial has 8 LC defects.

The Li/Kadrawi--Levit grids processed 2,500 trees in each of the two families. For `T_{3,m,n}`, 229 were non-log-concave; for `T^*_{3,m,n}`, 321 were non-log-concave. Both grids had 0 non-unimodal rows. The original 26-vertex witnesses `T_{3,4,4}` and `T^*_{3,3,4}` are regression tests: both have an LC break at `k = 13`, are `d_leaf <= 1`, and remain unimodal.

The small paper example `TG_{2,5}` is now a regression test: the generated tree has 70 vertices, independence number 37, LC breaks at indices 34 and 36, and a unimodal independence sequence.

After extending the harness with Conjecture-A stress metrics, all three current sources satisfy the same low-mode/mean-bound pattern:

| Source | Rows | Non-LC | Low-mode rows | `mu < n/3` rows | `mode <= ceil(mu)` rows | Min low-mode surplus | Min `n/3 - mu` | Max near-miss |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Ramos--Sun epoch 11 | 50,000 | 39,766 | 50,000 | 50,000 | 50,000 | 3 | 1.9606941767 | 0.8493940107 |
| Galvin grid | 759 | 390 | 759 | 759 | 759 | 0 | 0.0833333333 | 0.9808181930 |
| Bautista--Ramos grid | 64 | 28 | 64 | 64 | 64 | 1 | 0.3543123543 | 0.9729533142 |
| Li `T_{3,m,n}` grid | 2,500 | 229 | 2,500 | 2,500 | 2,500 | 1 | 0.5047704234 | 0.9507426453 |
| Li `T^*_{3,m,n}` grid | 2,500 | 321 | 2,500 | 2,500 | 2,500 | 1 | 0.6060606061 | 0.9467941207 |
| Exhaustive `n=28` LC failures | 19 | 19 | 19 | 19 | 19 | 1 | 0.8503646480 | 0.7059671963 |

The tight cases split into two different stress modes:

- Closest to `mu = n/3`: Galvin `(m,t) = (1,50)`, with `n = 102`, mode `34`, `mu = 33.8333333272`, low-mode surplus `1`, and no LC defect.
- Closest to tail-rise/non-unimodality: Galvin `(21,11)`, with near-miss `0.9808181930`, but its single LC defect is at `k = 233`, 79 positions after the mode and 78 positions after the near-miss transition.
- Multiple-break stress: Bautista--Ramos `(8,8)`, with 8 LC defects and near-miss `0.9729533142`; the first LC defect is 70 positions after the mode and 69 positions after the near-miss transition.
- Li-family stress: `T_{3,50,50}` has near-miss `0.9507426453`, LC defect at `k = 105`, and mode `69`; `T^*_{3,50,50}` has LC ratio about `3.25`, LC defect at `k = 106`, and mode `70`.
- Ramos--Sun top near-miss row has mode `17`, near-miss at transition `18 -> 19`, and its LC defect at `k = 30`, 13 positions after the mode.

Ratio-profile interpretation: LC breaks are bumps in the consecutive coefficient ratios `r_k = a_k/a_{k-1}`. In the current stress sources, the maximum post-mode ratio occurs immediately after the mode, while the LC bumps occur much deeper in the tail where the absolute ratios are tiny. For example, Galvin `(21,11)` has near-miss ratio `0.9808181930` near the mode, but at its LC defect the ratios jump only from about `7.27e-5` to `0.002808`. Bautista--Ramos `(8,8)` similarly jumps from about `4.07e-4` to `0.001292` at its first LC defect. This is a better invariant lead than "LC defects are localized": the observed invariant is **separation between tail-rise pressure and LC-defect pressure**.

The exhaustive `n=28` LC failures show the same qualitative separation in a tighter small-order form: every one of the 19 LC failures is `d_leaf <= 1`, low-mode, mean-bounded, and unimodal. Their LC defects occur only 5 or 6 positions after the mode, but the right-hand ratio at the LC bump is at most `0.0210526316`, far below the near-mode tail ratios. The top `n=28` near-miss tree is different: it is not `d_leaf <= 1`, fails low-mode, and has no LC defect, which is consistent with the PNP framework routing such cases through the multi-leaf-hub transfer rather than Conjecture A.

I then ran a bounded forest/product stress audit over single, pair, and triple products of nine representative hard components: Ramos--Sun top near-miss/top-LC, Galvin top near-miss/top-mean, Bautista--Ramos top near-miss, Li top rows, and the `n=28` top near-miss/top-LC witnesses. This produced 219 forest products, 188 of which were non-log-concave, with 0 non-unimodal rows. All 188 non-LC product rows still had LC defects after the mode and after the near-miss transition. The maximum near-miss ratio was `0.9923801835`, attained by the triple product of Bautista--Ramos `TG_{8,8}`. The `d_leaf <= 1` product subfamily had 164 rows; all 164 were low-mode and mean-bounded, 148 were non-LC, and 0 were non-unimodal. Low-mode/mean failures occurred only in products containing the non-`d_leaf <= 1` `n=28` top near-miss witness.

I then added a ratio-profile power audit, which is the more direct stress test for the separation invariant. It takes hard components and repeated products, with single-component powers up to `20` and mixed-component powers up to `4`. Across 184 product powers, 52 were non-log-concave and 0 were non-unimodal. The near-mode tail pressure became very tight: the maximum post-descent transition ratio reached `0.9988814019`. But the right-hand ratio at the strongest LC-defect bump never exceeded `0.0166666667`, every LC row had its LC-bump right ratio below the maximum post-descent tail ratio, and no row crossed the audit's `0.95` unsafe threshold. In other words, product powers can push the sequence arbitrarily close to a post-mode plateau without making LC defects coincide with the ratio pressure that would create a valley.

Finally, I added a direct beam search over hard-component products. This search retained states that maximized post-mode tail pressure, LC-bump pressure, and bump-to-tail pressure, rather than only following preset powers. With `max_factors = 8` and beam width `80`, it examined 2,866 products; 2,628 were non-log-concave, with 0 non-unimodal rows and 0 unsafe LC-bump rows. The search found mixed products with tail ratios up to `0.9971493049`, but it did not beat the simple product-power scan, and the maximum LC-bump right ratio again stayed at `0.0166666667`. Mixing hard components therefore made LC failure dense but did not couple LC defects to valley pressure.

I then tested the first proposed proof mechanism for the separation lemma: at each post-descent LC bump, set `lambda = a_k/a_{k+1}` so sizes `k` and `k+1` tie, then compute `mu(lambda)-k`. This falsified the broad tie-fugacity mechanism. Across the same 3,050 product rows, there were 8,777 post-descent LC-bump events, 6,590 had negative `mu(lambda)-k`, and the minimum was `-3.8776208619` at Bautista--Ramos `TG_{8,8}`. There were still 0 post-descent upward-transition events. So the ratio separation survives, but the mean/tie route cannot be applied to all LC bumps. It can only be considered for actual upward crossings, which remain absent in the stress set.

Compact audit artifacts: `results/stress_invariant_audit_2026-07-03.json`, `results/forest_product_stress_audit_2026-07-03.json`, `results/ratio_profile_power_audit_2026-07-03.json`, `results/product_ratio_beam_search_2026-07-03.json`, and `results/tie_fugacity_bump_audit_2026-07-03.json`.

## Li--Li--Yang--Zhang 2025
Source: Ethan Y. H. Li, Grace M. X. Li, Arthur L. B. Yang, Zhong-Xue Zhang, "A symmetric function approach to log-concavity of independence polynomials," arXiv:2501.04245v1, 2025.
### What They Prove
- They define 2-Schur positivity for symmetric functions by nonnegativity of Schur coefficients indexed by partitions of length at most 2.

- Their Theorem 1.3 gives an equivalence: for a positive-coefficient polynomial `P`, 2-Schur positivity of the associated `F_P` is equivalent to log-concavity and strong log-concavity of `P`.

- For a graph `G`, this translates to 2-Schur positivity of the associated `Y_G` being equivalent to log-concavity of `I_G(t)`.

- They prove all spider independence polynomials are strongly log-concave.

- They also prove pineapple graphs have log-concave independence polynomials.

### Consequence for This Project
Global 2-Schur positivity cannot be the missing invariant between log-concavity and unimodality. It is log-concavity in richer language. The reusable part is the proof technique: pair negative 2-Schur terms with positive partners, or try to prove only the subset of Schur coefficients needed before the known decreasing tail takes over.
## Li 2026
Source: Grace M. X. Li, "Unimodality of independence polynomials of two family of trees," arXiv:2603.03025v1, 2026.
### What She Proves
- For every `m,n >= 1`, `I(T_{3,m,n}; t)` is unimodal.

- For every `m,n >= 1`, `I(T^*_{3,m,n}; t)` is unimodal.

- These families include the original Kadrawi--Levit 26-vertex log-concavity failures and their natural infinite extensions.

### Proof Shape
- Use the 2-Schur framework to show prefix log-concavity up to just before the dangerous top index.

- Use Levit--Mandrescu tail monotonicity to finish the decreasing tail.

- The method is direct-unimodality in effect, but the machinery still works by proving enough log-concavity on a prefix and then importing a known tail theorem.

### What Is Not Proved
- No general theorem for arbitrary trees.

- No extension to all non-log-concave tree families.

- No preserved invariant for the standard tree recursion.

## Bautista-Ramos 2025
Source: César Bautista-Ramos, "Multiple breaks of log-concavity in the independence polynomials of trees," arXiv:2511.00334v1, 2025.
### Relevance
This paper constructs infinite families of trees whose independence polynomials violate log-concavity at an arbitrary number of indices. It kills the naive "one localized LC defect" version of Route 1. It does not, by itself, threaten tree unimodality.
## Strategy Implications
### Route 1: Corpus Mining
Proceed, but redefine the target. Do not look for "LC except within bounded distance of a canonical position." That target is already too fragile. Instead:

- Compute valley margins, mode position, LC-defect locations, `alpha`, degree-2 abundance, leaf-support patterns, `d_leaf`, and mean/mode gaps across the Ramos--Sun corpus.

- Ask whether LC defects, even when multiple, remain separated from actual descent-then-rise patterns by a stable reserve inequality.

- Compare Ramos--Sun examples with Galvin and Bautista-Ramos structured families to avoid overfitting to PatternBoost's `N/2` bias.

- Test any candidate invariant under `I_T = I_{T-v} + x I_{T-N[v]}` and under the subdivision identity `I(T_e) = I(T) + x I(T/e)`.

### Route 2: Symmetric-Function Lifting
Treat as background machinery, not the main project lane. A plausible subproject is:

- Formalize "partial 2-Schur positivity through index `q` plus Levit--Mandrescu tail monotonicity implies unimodality."

- Search for a graph/tree condition that implies that partial 2-Schur positivity without implying full log-concavity.

- Test whether the PNP or `d_leaf <= 1` reductions naturally imply such a partial condition.


This should not displace the current PNP/mode-mean architecture unless it produces a concrete, checkable invariant.
### Route 3: Probabilistic / Hard-Core
Still relevant to the local project. The current `mu(T) < n/3` and mode-mean direction is closer to this route than the external papers are. The priority remains proving or refuting the mode-mean inequality in the class needed by Conjecture A.
### Counterexample Search
Ramos--Sun shows modern search can break stronger hypotheses, not that it has seriously attacked non-unimodality. A non-unimodality search needs a global reward:

- penalize first descent followed by later rise;

- track valley depth and recovery height;

- include products/forests, since forest polynomials can amplify shape interactions;

- seed with non-LC families but avoid only optimizing single LC defects.

## Recommended Next Actions
1. Treat the Ramos--Sun GitHub corpus as normalized through the current summary artifacts, and expand only when a candidate invariant needs row-level witnesses.

2. Draft the separation invariant suggested by the current audit: after the first descent, the maximum ratio `a_{j+1}/a_j` occurs near the mode, while LC-defect ratio bumps occur far enough right that both adjacent ratios remain well below 1.

3. Use the added Galvin and Bautista--Ramos generators to stress-test any proposed invariant against single-break, multi-break, and PatternBoost-discovered non-LC examples.

4. Do not pursue a tie-fugacity proof for all LC bumps. The next proof target should be crossing-only: show that after the first descent, `r_{k+1} >= 1` cannot occur, using mean reserve, tail lock, or structural transfer.

5. Test crossing-only sufficient conditions against coefficient-level artifacts before trying to prove them. The current compact JSONs are enough for ratio conclusions, but not enough for all proposed mean/tie mechanisms.

## Project Decision
The external literature supports keeping the public posture conservative:

- no proof of Erdős 993;

- unimodality verified through `n <= 29`;

- log-concavity is decisively false;

- symmetric-function methods prove important families but are not a general bridge;

- the most useful AI-heavy work is invariant mining plus formal/certificate stress testing, not broad proof search.
