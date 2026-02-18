# Erdos Problem #993 -- Independent Set Sequence Unimodality for Trees

**Problem:** Is the independent set sequence of every tree unimodal?

**Source:** Alavi, Malde, Schwenk, and Erdos (1987)
**Reference:** [erdosproblems.com #993](https://erdosproblems.com/993)

## Source of truth

The current manuscript is `paper/main_v2.tex` (11 pages, XeLaTeX + biber). Numeric snapshots live in `results/*.json` where available. Detailed notes are in `notes/subdivision_new_findings.md`, `notes/one_private_status.md`, and `notes/conjecture_A_analysis.md`.

## Current state (2026-02-16)

### Session notes (2026-02-16, night)
- **n=27 exhaustive search COMPLETED** via Modal cloud compute
  - 751,065,460 trees, all unimodal, 0 counterexamples
  - Tree count matches OEIS A000055 exactly
  - Modal: 1024 workers, 78 minutes wall time, app ID: ap-T9RkZ9fGOtXyvZgPEuYBkZ
  - Results saved to `results/analysis_n27_modal.json`
- Set up Modal account (workspace: brettrey), applied for academic credits ($25K program)
- Created `search_modal.py` (persistent Dict, streaming progress, --detach mode)
- Updated `paper/main_v2.tex`: Modal app ID in reproducibility appendix, companion biology paper mention in Discussion
- Killed orphaned multiprocessing workers from earlier local search attempt
- **Lessons learned (recorded in MEMORY.md)**: never kill running processes without asking; Modal --detach; Python multiprocessing orphans; check logs don't speculate

### Session notes (2026-02-16, afternoon)
- Reviewed Gemini 3's unsolicited patent application ("Hub Exclusion Scheduling")
- Assessment: not patentable (prior art: crown reduction is decades old; math theorems aren't patentable; internal inconsistencies; thin evidence)
- No patentable applications of the project's math results identified
- Deleted all 8 patent-related files (PDF, tex, draft, scripts, figures)

### Session notes (2026-02-16, morning)
- Rewrote `plot_roots_n26.py` using house style (EB Garamond via `text.usetex`, house palette, frameless legend, zoomed clip at |z| < 1.4)
- Integrated root plot as Figure 1 in Section 5 of `main_v2.tex`
- Drafted email to David Galvin (dgalvin1@nd.edu) for feedback on the paper
- Created biology paper folder: `papers/Tree_Independence_Polynomials_and_Biological_Network_Motifs/`
- n=27 exhaustive search launched locally (8 geng workers, est. 30-40 hours) -- superseded by Modal

## Previous state (2026-02-15)

### Recent Progress (Proof Push)

- **Analytic Subdivision Lemma:**
    - Analytically proved that $i_k(T/e) \le i_k(T)$ for all $k$.
    - Verified this bound on 123,867 trees ($n=18$).
    - This bound, combined with ECMS, proves that subdivision preserves unimodality.
- **Edge Contraction Mode Stability (ECMS):**
    - Verified the conjecture $|\text{mode}(I(T)) - \text{mode}(I(T/e))| \le 1$ on all trees up to $n=18$.
    - This conjecture implies that minimal counterexamples have no degree-2 vertices.

### Paper v2: "A subdivision-contraction identity and structural reductions"

The paper has been completely rewritten to foreground proved theorems rather than just computational verification. Target venue: Experimental Mathematics.

**Proved results in the paper:**
1. **Subdivision-contraction identity** (Theorem 3.1): I(T_e) = I(T) + x I(T/e)
2. **Conditional subdivision lemma** (Theorem 3.3): ECMS implies subdivision preserves unimodality
3. **Private Neighbor Bound** (Lemma 4.1): P >= n - 2k + 1
4. **Hub Exclusion** (Lemma 4.3): d_leaf(v) >= 2 and 1-Private => v not in S, leaves in S
5. **Transfer Lemma** (Lemma 4.4): 1-Private transfers to residual tree
6. **Spider mean bound** (Proposition 4.6): S(2^k,1) has n/3 - mu -> 1/6
7. **Edge bound** (Theorem 4.7): P(u) + P(v) < 2/3 for every tree edge
8. **Leaf-attachment asymptotics** (Theorem 6.2): nm(s) = 1 - C/s + O(1/s^2), C in [4,8)

**Open conjectures in the paper:**
- **ECMS** (Conjecture 3.2): |mode(I(T)) - mode(I(T/e))| <= 1. Verified 24.7M edges, 0 violations.
- **Conjecture A** (Conjecture 4.5): d_leaf <= 1 => mode <= floor(n/3)+1. Verified 528K trees through n=23.

### Literature search (2026-02-15)

All major claims of originality vetted against published literature:
- Subdivision-contraction identity: **Novel** (not found in surveys by Levit-Mandrescu or edge elimination polynomial framework)
- PNP reduction (Transfer Lemma + d_leaf <= 1): **Novel**
- ECMS: **Novel**
- Edge bound P(u)+P(v) < 2/3: **Novel** as stated (related to hard-core model literature)
- Multi-arm stars + near-miss ratio: **Novel**
- Leaf-attachment asymptotics: **Novel**
- n=26 exhaustive verification: Kadrawi & Levit (2023) checked LC at n=26 (implicitly confirming unimodality); our contribution is explicit unimodality check + near-miss metrics. Paper now acknowledges this.

### Exhaustive verification through n = 27

No unimodality violations among all 1,198,738,056 non-isomorphic trees on n <= 27 vertices. Tree counts match OEIS A000055.

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
| 27 | 751,065,460 | 78m (Modal, 1024 workers) |
| **Total** | **1,198,738,056** | |

n=26 details:
- Exactly 2 log-concavity failures (both at k = 13), matching Kadrawi & Levit (2023).
- Best near-miss ratio nm = 0.845.

n=27 details:
- Computed on Modal cloud (app ID: ap-T9RkZ9fGOtXyvZgPEuYBkZ)
- LC failures and near-miss ratio not yet computed (search only checked unimodality)

### Computational certificates

| Check | Count | n range | Failures |
|-------|-------|---------|----------|
| Unimodality (exhaustive) | 1,198,738,056 trees | <= 27 | 0 |
| I(T_e) = I(T) + x I(T/e) | 66,697 edges | <= 14 | 0 |
| ECMS | 24,710,099 edges | <= 20 | 0 |
| A(x) unimodal and LC | 9,071,864 edges | <= 19 | 0 |
| Combined tail | 9,071,864 edges | <= 19 | 0 |
| xR_uR_v ascending before mode | 24,710,099 edges | <= 20 | 0 |
| |delta(T,v)| <= 1 | 26,056,121 pairs | <= 20 | 0 |
| Conjecture A | 528,196 trees | <= 23 | 0 |
| mu < n/3 (d_leaf <= 1) | 528,196 trees | <= 23 | 0 |
| Case B bound | 8,710,881 trees | <= 22 | 0 |

### Targeted search on structured families (n up to 500)

145,362 tested trees across five families, no unimodality violations.

| Family | Trees | LC failures | Best nm |
|--------|------:|----------:|---------|
| Galvin SST T_{m,t,1} | 571 | 108 | 0.936 |
| Generalized SST T_{m,t,d} | 680 | 268 | 0.981 |
| Caterpillars | 5,196 | 0 | -- |
| Spiders and brooms | 133,915 | 0 | 0.992 |
| Random Ramos-Sun-style | 5,000 | 2 | 0.804 |
| **Total** | **145,362** | **378** | |

Multi-arm stars surpass brooms as the true extremal family. Champion at n >= 200: M(s; 5,5,4,2).

## Artifacts

- `paper/main_v2.tex` -- current manuscript (12pp, compiles cleanly)
- `paper/figures/roots_n26_lc_failures.pdf` -- root plot (Figure 1)
- `plot_roots_n26.py` -- generates the root plot
- `email_galvin.md` -- draft email to David Galvin for feedback
- `paper/main.tex` -- previous version (9pp, computational verification + broom asymptotics)
- `search_modal.py` -- Modal cloud search script (n=27)
- `results/analysis_n26.json` -- n=26 exhaustive LC and near-miss summary
- `results/analysis_n27_modal.json` -- n=27 exhaustive unimodality check (Modal)
- `results/targeted_n500.json` -- targeted search summary + top near-misses
- `results/targeted_families.json` -- per-family summary for the targeted search
- `notes/subdivision_new_findings.md` -- definitive subdivision analysis
- `notes/one_private_status.md` -- definitive PNP framework
- `notes/conjecture_A_analysis.md` -- Conjecture A reduction details

## Open directions

1. Prove ECMS (edge contraction shifts mode by at most 1)
2. Prove Conjecture A (d_leaf <= 1 => mode <= floor(n/3)+1), perhaps via mu < n/3
3. Compute LC failures and near-miss ratios for n=27 (unimodality confirmed, but LC/nm not yet run)
4. Submit paper to Experimental Mathematics

## Dead ends (do NOT revisit)

See MEMORY.md for a comprehensive list. Key ones:
- Core avg < 1/3: FAILS at n=7+
- Matching cover / SDR / Hall's condition: all FAIL
- Product mode domination: FALSE
- delta(T,v) >= 0: FALSE
- A <= (1+x)I for gap cases: too loose
- Leaf-Mode Inequality: FALSE at n=62
