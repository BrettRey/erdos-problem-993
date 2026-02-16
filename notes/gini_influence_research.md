# Influence Gini Coefficient - Research Notes

## Summary

We explored whether the "influence concentration" of vertices in a tree correlates with unimodality near-miss behavior. The idea was to apply Gini coefficient concepts from inequality measurement to understand why certain trees get close to violating unimodality.

## Approach

1. **Definition**: For each vertex v, compute S(v) = {k : v can be in an independent set of size k}. Let a(v) = average of S(v). Then G(T) = Gini({a(v)}).

2. **Hypothesis**: Higher G → higher near-miss ratio (closer to unimodality violation).

## Key Data

| Tree Type | G (Gini) | nm (near-miss) |
|-----------|-----------|-----------------|
| Broom(10,500) | 0.89 | 0.985 |
| Star(200) | 0.50 | 0.97 |
| Path(100) | 0.33 | 0.84 |

## Findings

- **Brooms**: High G → high nm (confirms hypothesis)
- **Stars**: Moderate G, high nm (breaks simple monotonic relationship)
- **Paths**: Low G → lower nm

The correlation is real but not clean enough for a simple theorem.

## Files

Created in branch `causal-influence-analysis`:
- `timing_analysis.py` - timing metrics
- `causal_influence.py` - influence cone computation
- `support_analysis.py`, `quantify_gini.py` - correlation analysis  
- `exact_gini.py` - verified exact values
- `THEOREM_GINI_UNIMODALITY.txt` - theorem statement draft
- `formal_proof.py` - proof sketch

## Verdict

Worth publishing as **empirical observation/conjecture** but not as a formal theorem. The relationship is too messy for a clean mathematical result.

## Conclusion

The structural intuition is valuable (concentrated influence → near-miss), but the quantitative relationship doesn't yield a clean theorem. Could be mentioned in paper as "observed correlation between influence concentration and near-miss behavior."
