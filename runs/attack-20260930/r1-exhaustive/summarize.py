"""Print a compact per-n summary table from results.json (diagnostic floats labelled)."""
import json
import sys

d = json.load(open(sys.argv[1]))
tot = {"trees": 0, "window_levels": 0, "window_vk": 0, "r1_viol": 0, "pv_viol": 0, "r1_zero": 0}
for n, r in d["per_n"].items():
    if "stats" not in r:
        print(n, r)
        continue
    s = r["stats"]
    for k in tot:
        tot[k] += int(s[k])
    pv = r["pv_max_ratio_lhs_over_rhs"]
    t1, t2 = r["TOP_NORM1"][0], r["TOP_NORM2"][0]
    print(f"n={n} trees={s['trees']} oeis_ok={r['count_ok']} src={r['source_dir']} same={r['r1_census_vs_r1_ext_inwindow_identical']} "
          f"levels={s['window_levels']} (u,k)={s['window_vk']} R1viol={s['r1_viol']} R1zero={s['r1_zero']} PVviol={s['pv_viol']} PVtrees={s['pv_viol_trees']}")
    print(f"   PVmax={pv['exact']} (float {pv['float_diag']:.5f}) k={pv['k']} deg={pv['deg']} par={pv['par']}")
    print(f"   norm1 max={t1['norm1_exact']} (float {t1['norm1_float_diag']:.5f}) k={t1['k']} W={t1['window']} alpha={t1['alpha']} deg_u={t1['deg_u']} par={t1['par']}")
    print(f"   norm2 max={t2['norm2_exact']} (float {t2['norm2_float_diag']:.5f}) k={t2['k']} W={t2['window']} alpha={t2['alpha']} deg_u={t2['deg_u']} par={t2['par']}")
    if "outside_window_ext" in r:
        e = r["outside_window_ext"]
        print(f"   EXT {e['stats']} hist={e['hist']} min={e['r1_min_offset_above_q']}")
print("TOTAL", tot)
print("violations rechecked:", len(d["violations"]))
