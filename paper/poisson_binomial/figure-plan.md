# Figure plan: Variance and local log-concavity of Poisson–binomial laws
<!-- SUMMARY: plan-figures menu for the ECP paper (12-page cap, PDF at 12); answer to "should it have figures?": not required, one candidate earns its place if space is made · status: menu, awaiting Brett's trim · updated: 2026-10-03 -->

**Brett's question: should this paper have any figures at all?** No figure is required. ECP proof papers often have none, and every claim in this paper is exact. One candidate earns its place: the comparison of lower bounds in the balanced variance-one family, #1 below. The ChatGPT Pro referee named that comparison as the strongest case for the paper. Two cheap ways to show it:
- **(a)** As the small table #8. That costs about six lines and probably keeps the paper at 12 pages.
- **(b)** As figure #1. That costs a third to half a page. It would push the paper to 13 pages, which the 2025 IMS editorial report tolerates, unless something is cut. The obvious cut is the 12-line smoothing proof of (3.3), replaced by the citation to BMM Corollary 3.2.

Nothing below has been built.

**Sources.** Every candidate plots a quantity defined in the paper and computed exactly or to certified precision. None plots empirical data. Where no saved output exists, the entry says `COMPUTE`: a short script must produce the numbers from the stated formulas before the figure is drawn, and no number goes on a plot from memory. Existing outputs:
- `scripts/verify_pb_corollaries_20261003.py`: random laws, ULC family values, CUE values;
- `runs/pb-revision2-20261003/chatgpt-pro-independent-checks/audit_cue.json`: certified CUE values.

| # | Fig | Kind | Makes clearer | Type | Source | Keep |
|---|-----|------|---------------|------|--------|------|
| 1 | Lower bounds on δ_D in the balanced V=1 family, against m | data (computed) | That every earlier bound degenerates at fixed variance while Theorem 1.1 holds. Compared to what: the four earlier bounds plotted against 1/4 and the true δ_D | log–log lines with direct labels: ULC (∼2/m), 1/(D+1) (∼1/m), Johnson (O(m⁻²)), Pitman (20) (∼√2/m), true δ_D, horizontal line at 1/4 | `COMPUTE`. The formulas are in §1; the Pitman values were checked in session at m = 5, 10, 40, 200, 1000, and the true δ_D needs an exact pmf | **must, if any figure** |
| 2 | Bound at D for the CUE half-circle count, against N | data (computed) | That variance and degree give different scales: 1/log N against 1/N | log-x lines with direct labels: 1/(4V_N) and (N+1)/((N/2+2)N/2), with the even threshold N = 1998 marked | `COMPUTE` from the exact V_N formula; spot values in `audit_cue.json` | nice |
| 3 | 1/δ_k against k for one Bernoulli sum, with the cone 1/δ_D + \|k−D\| | data (computed) | What (2.5) and Corollary 1.3 say: the reciprocal deficit moves at most one per step, so the bound at D spreads | points for 1/δ_k; two slope-±1 lines from (D, 1/δ_D); a horizontal guide at 4V | `COMPUTE` (one explicit law, exact arithmetic) | nice |
| 4 | Same axes for the ULC law g^(m) next to a Bernoulli sum of equal variance | data (computed) | What ULC permits and Bernoulli sums forbid: g^(m) jumps from 1/δ_m ≤ 4/3 to 1/δ_{m+1} ≈ m/2 | two small multiples on one scale | `COMPUTE` (the closed forms of Proposition 1.5) | nice, best paired with #3 |
| 5 | Proof architecture | conceptual DAG | Which step uses which, and which steps are analytic, computer-checked, or conditionally formalized. Arrow meaning: "is used to prove" | TikZ DAG: HJ cubic inequalities → (2.5) → Lemma 2.1 → mass bounds → V ≥ M²A(δ); BMM → M² ≥ (1+12V)⁻¹; both → Prop. 3.1 ← {Bernstein certificates, symbolic J ≥ 6}; Prop. 3.1 → Thm 1.1 → Cor. 1.3, Prop. 1.4 | conceptual | stretch |
| 6 | Slack in the scalar inequality, ST/Q(H) and A(δ)/Q against H | data (computed) | Where the proof is tight (small H, near δ = 1/4) and where it has room (large H, ratio → 8π/3). Bears on "where does 1/4 come from?" | single line, log-x, reference line at 1 | `COMPUTE` from (4.4)–(4.6); the large-H limit is heuristic and needs checking | stretch |
| 7 | V·δ_D for random Bernoulli sums, with lines at 1/4 and 1/3 | data (computed) | The bracket (1.4) in context | scatter, range-framed | Partly in `verify_pb_corollaries_20261003.py`; the output would need saving | stretch. **Caution:** random laws never come near the extremal families (sampled minimum ≈ 0.50), so the plot could mislead about where κ⋆ lies |
| 8 | Table, not figure: the bounds of #1 at D, with their rates as m → ∞ | table | The same comparison as #1, in about six lines | 5-row table (ULC ∼2/m, 1/(D+1) ∼1/m, Johnson O(m⁻²), Pitman (20) ∼√2/m, Theorem 1.1 ≥ 1/4) | Rates from §1 and the ChatGPT Pro report. The paper states only the limits, and the Pitman rate is numerical, so the table would say "tends to 0" or give the rate with its basis | **must-alternative to #1** |

Dropped after the house-doctrine check:
- a Bernstein-coefficient table (the lab-notebook material removed in the first revision);
- any "pmf with mode and first descent marked" schematic, which is decoration.

**To trim:** pick #8 (table, no page cost) or #1 (figure, page cost); say whether #3 and #4 are worth a panel pair; the rest can go. Every survivor needs only a `COMPUTE` script, no data collection.
