# Decisions Log

Append-only record of project decisions. Agents: add an entry whenever a non-trivial decision is made during a session (structural changes, venue choices, theoretical commitments, scope changes, reviewer feedback acted on). Keep entries short.

Format: `## YYYY-MM-DD` then bullet points with **bold topic** and brief rationale.

---

## 2026-03-11

- **Prospect-aware search scoring** Added archive-plus-prospect scoring to `nm_optimizer.py` and `scripts/lc_breaker_optimizer.py` rather than trying to reconstruct unavailable AlphaEvolve code; the transferable idea is to reward lineages that empirically produce stronger children while retaining best distinct exact trees.
- **Ramsey paper citation placement** Added the Nagda--Raghavan--Thakurta Ramsey/AlphaEvolve paper to `paper/references.bib` and cited it in `paper/main_v2.tex` as broader AI-assisted extremal-combinatorics context, not as direct tree-unimodality prior art.
- **No `n=29` LC/NM burn** Estimated Modal cost for a full `n=29` LC + near-miss run was on the order of $1k before credits, which is outside budget; an initial dispatch was stopped immediately and no sustained large compute will be pursued.
- **Project shelved** With exhaustive unimodality already closed through `n=29` and no affordable high-value computation left, the project is being put to bed; only low-cost paper/admin cleanup would justify reopening.
