import json
from graph6 import parse_graph6
from timing_analysis import compute_timing_imbalance
import sys
sys.path.insert(0, '.')
from indpoly import independence_poly, is_log_concave, log_concavity_ratio

with open('results/analysis_n26.json') as f:
    data = json.load(f)

for failure in data['lc_failures']:
    g6 = failure['graph6']
    print(f"Graph6: {g6}")
    print(f"  Poly: {failure['poly']}")
    
    # Try to decode (may not work with all formats)
    try:
        G = parse_graph6(g6)
        n = len(G.nodes())
        adj = [list(G.neighbors(v)) for v in range(n)]
        timing = compute_timing_imbalance(n, adj)
        print(f"  n={n}, gini={timing['gini']:.3f}, min_spread={timing['min_spread']}")
    except Exception as e:
        print(f"  Could not decode: {e}")
    print()
