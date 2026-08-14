#!/bin/bash
# LC-failure census ladder, n = 27..32 (2026-08-14).
# Validation gates: n=27 must yield 0 failures; n=28 must reproduce the 19
# stored failures by exact polynomial multiset. Any gate failure or ALARM
# aborts the ladder. Restartable: completed shards are skipped by the driver.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

for n in 27 28 29 30 31 32; do
  echo "=== ladder: n=$n start $(date) ==="
  python3 scripts/run_lc_census_20260814.py --n "$n" --mod 64 --workers 9
  python3 - "$n" <<'EOF'
import json, sys
sys.path.insert(0, '.')
n = int(sys.argv[1])
s = json.load(open(f'results/lc_census_20260814/summary_n{n}.json'))
assert s['complete'], f'n={n} INCOMPLETE: {s["shards_missing"]}, trees={s["trees"]} vs {s["expected_trees"]}'
assert s['alarms'] == 0, f'n={n}: ALARM lines present — inspect summary NOW'
if n == 27:
    assert s['failures'] == 0, f'n=27 gate: expected 0 failures, got {s["failures"]}'
    print('gate n=27 passed: 0 failures, complete')
elif n == 28:
    import networkx as nx
    from verify_partial_sync_obstruction_20260811 import tree_ipoly
    stored = json.load(open('results/analysis_n28_modal_lc_nm.json'))
    want = []
    for f in stored['top_lc_failures']:
        p = f['poly']
        if p is None:
            G = nx.from_graph6_bytes(f['graph6'].encode())
            adj = [list(G.neighbors(v)) for v in range(G.number_of_nodes())]
            p = tree_ipoly(len(adj), adj, 0)
        want.append(tuple(p))
    got = []
    for line in s['fail_lines']:
        d = dict(x.split('=', 1) for x in line.split()[1:])
        got.append(tuple(int(c) for c in d['poly'].split(',')))
    assert sorted(got) == sorted(want), f'n=28 gate: polynomial multiset mismatch ({len(got)} vs {len(want)})'
    print('gate n=28 passed: 19 failures, exact polynomial multiset match')
else:
    print(f'n={n} complete: failures={s["failures"]}, min dist in fail lines to inspect at merge')
EOF
done
echo "=== ladder complete $(date) ==="
