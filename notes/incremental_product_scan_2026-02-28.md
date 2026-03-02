> **OBSOLETE (2026-03-01):** SCC / Condition C fails at n=28. The SCC invariant propagation reported here is valid only for n≤22. See `notes/scc_false_n28_2026-03-01.md`.

# Incremental Product Scan Results (2026-02-28)

## Purpose

Validates the "incremental factor induction" approach (GPT 5.2 Round 5, Prompt 1/Instance 1).
At each support vertex r of each tree, we build E^{(k)} and J^{(k)} incrementally
by multiplying in one child factor at a time:

```
E^{(0)} = 1,  J^{(0)} = 1
E^{(k+1)} = E^{(k)} · I_c,   J^{(k+1)} = J^{(k)} · E_c
```

At every intermediate stage we check:
1. **Strong Condition C**: b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0
2. **LC of E^{(k)}**: b_k² ≥ b_{k-1}·b_{k+1}
3. **J^{(k)} ≤ E^{(k)}**: coefficientwise dominance

## Results (n ≤ 22, ALL trees)

| Metric | Value |
|--------|-------|
| Total trees | 9,114,283 |
| Total intermediate stages | 89,888,083 |
| Condition C failures | **0** |
| LC failures (E^{(k)}) | **0** |
| J^{(k)} > E^{(k)} failures | **0** |
| Min nontrivial margin | **6** |
| Achieved at | n=6, g6=E?qo, stage 1/1 |
| Runtime | 1119.2s |

## Nontrivial margin definition

The raw minimum margin is 0 (achieved everywhere from n=4 onward) because at
indices k beyond deg(E), all coefficients are zero and the SCC identity holds
as 0 ≥ 0 trivially. The "nontrivial margin" excludes indices where b_k = 0,
measuring the tightest genuine constraint.

## Minimum margin tree

The graph6 string `E?qo` at n=6 achieves min nontrivial margin = 6. This tree
has a support vertex with a single child (stage 1/1), and the SCC margin at
the tightest interior index is exactly 6.

## Significance

- **All three invariants propagate through every intermediate product stage.**
- The margin of 6 is substantial (not approaching 0), suggesting significant slack.
- The incremental framework E^{(k+1)} = E^{(k)}·I_c, J^{(k+1)} = J^{(k)}·E_c
  preserves all the structure needed for an inductive proof.
- J^{(k)} ≤ E^{(k)} holding at ALL stages (not just the final product) is a
  key new finding — it means the dominance constraint is available as a hypothesis
  at every induction step.

## Leaf-factor algebra (same session)

See `notes/leaf_factor_SCC_closure_2026-02-28.md` for the full decomposition.

Key identities for attaching a leaf (I_c = 1+x, E_c = 1):
- d'_k = d_k + e_k  where e_k = b_{k-1}·j_k - b_k·j_{k-1}
- c'_k = c_k + c_{k-1} + γ_k  where γ_k = b_k·b_{k-1} - b_{k-2}·b_{k+1} ≥ 0 (LC)
- S'_k = S_k + (α+β') + β + γ₁ + γ₂

Profile (523K checks, n≤15): β < 0 in 95.3%, α+β' < 0 in 0.01%, min ratio nonneg/|neg| = 3.0.

## Intermediate stage margin analysis

The global minimum nontrivial SCC across ALL intermediate stages (n≤17, 776K stages,
48.6K trees) is **exactly 2**, achieved at the very first stage (one leaf child only):
E = [1,1], J = [1], I = [1,2]. For stars K_{1,s}: min SCC = s at k=s (grows with s).

The worst-case trees at each n are always stars K_{1,n-1}.

## Script

`scan_incremental_product.py`
