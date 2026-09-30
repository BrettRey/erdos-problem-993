"""Collect exact data for the reported R1 counterexamples (from the independent rechecks) plus r1lib diagnostics."""
import json
from fractions import Fraction
from r1lib import adj_from_edges, evaluate
out = []
for f in ['data/recheck_mixed_4x9_4x10.json', 'data/recheck_S_8_10_2.json', 'data/recheck_S_9_9_2.json']:
    d = json.load(open(f))
    adj = adj_from_edges(d['n'], [tuple(e) for e in d['edges']])
    r = evaluate(adj)
    lv = []
    for l in d['levels']:
        L = Fraction(l['L_centre'])
        lv.append({"k": l['k'], "L_centre_exact": l['L_centre'], "L_centre_float": float(L),
                   "lc_margin_exact": l['lc_margin_ik2_minus_prod'], "lc_margin_rel": l['lc_margin_rel'],
                   "normA_L_over_k_ik1_ik": l['L_over_k_ik1_ik']})
    out.append({"source": f, "n": d['n'], "alpha": d['alpha'], "window": d['window'],
                "hubs_m": d.get('hubs_m', [d.get('m')] * d.get('h', 0)),
                "levels": lv, "r1lib_viol": [(v['u'], v['k']) for v in r['viol']],
                "r1lib_F": r['F'], "r1lib_F_at": r['F_at'], "r1lib_C": r['C'], "r1lib_C_at": r['C_at'], "r1lib_lc_fail": r['lc_fail']})
json.dump(out, open('data/counterexamples_summary.json', 'w'), indent=1)
for o in out:
    print(o['n'], o['alpha'], o['window'], o['hubs_m'], o['r1lib_viol'], 'F', round(o['r1lib_F'], 4), 'C', round(o['r1lib_C'], 4), o['r1lib_lc_fail'])
    for l in o['levels']:
        print('   k', l['k'], 'L', l['L_centre_float'], 'LCmargin', l['lc_margin_exact'], 'rel', l['lc_margin_rel'], 'normA', l['normA_L_over_k_ik1_ik'])
