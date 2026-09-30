"""Summarise jsonl result files: global max lower end, min upper end, per-family bests,
S(2/3) violations and empty rows. Floats for display only; exact ends in the records."""
import json, glob, sys, collections
files = sys.argv[1:] or glob.glob('data/fam_p*.jsonl')
R = [json.loads(l) for f in files for l in open(f)]
print(len(R), 'trees; S(2/3)-violating trees:', sum(1 for r in R if r['n_viol23']), '; empty-row trees:', sum(1 for r in R if r['empty_row']))
fam = collections.defaultdict(list)
for r in R: fam[r['family']].append(r)
for f, L in fam.items():
    a = max(L, key=lambda r: r['low_f']); b = min(L, key=lambda r: r['up_f']); c = max(L, key=lambda r: r['maxL23'])
    print(f"{f:10s} N={len(L):5d} maxlow={a['low_f']:.5f} {a['params']} n={a['n']} at{a['low_at']} | minup={b['up_f']:.5f} {b['params']} n={b['n']} at{b['up_at']} | maxL23={c['maxL23']:.3e} {c['params']}")
a = max(R, key=lambda r: r['low_f']); b = min(R, key=lambda r: r['up_f'])
print('GLOBAL interval [%.6f, %.6f]' % (a['low_f'], b['up_f']), a['family'], a['params'], '|', b['family'], b['params'])
