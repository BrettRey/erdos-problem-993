import hashlib, json, subprocess, glob, os, sys
sys.set_int_max_str_digits(0)
def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
files = sorted(glob.glob('*.py') + glob.glob('data/*'))
key = dict(
  verdict_S23='REFUTED (exact, three independent implementations agree on H(75,5))',
  S23_counterexample=dict(tree='H(75,5): hub 0 -> 75 supports -> 5 leaves each', n=451, edges='data/H_75_5_edges.json',
      upper_end=json.load(open('data/recheck_H_75_5.json'))['up'][:2], upper_end_float=0.665344, binding='hub u=0 (deg 75), k=222',
      violating_rows='u=0, k=220..224 (5 rows); LC holds at these k'),
  constant_rule_empty=dict(tree='MSH(38;11,2): centre 0 -> 38 hubs -> 11 supports -> 2 leaves', n=1293, edges='data/MSH_38_11_2_edges.json',
      lower=json.load(open('data/recheck_MSH_38_11_2.json'))['rows'][0], upper=json.load(open('data/recheck_MSH_38_11_2.json'))['rows'][1],
      certified_by=['orbit.py', 'slib.py (explicit rerooting)', 'recheck_rows.py (fresh forest DP, Fractions)']),
  cross_tree_empty=dict(lower_tree='MSH(25;11,2) n=851 lower 0.669288 (centre, k=383; rechecked)', upper_tree='H(75,5) n=451 upper 0.665344'))
json.dump(key, open('data/key_results.json', 'w'), indent=1, default=str)
files.append('data/key_results.json')
with open('manifest.yaml', 'w') as f:
    f.write('# lp-discharge adversary run, 2026-09-30\nmodel: claude-opus-5-5 (Claude Code subagent)\n')
    f.write(f'git_head: {head}\npython: venv/bin/python (python-flint 0.9.0)\ncores_used: <=3\n')
    f.write('arithmetic: exact integers/Fractions for every violation and interval claim; floats only in *_f ranking fields\n')
    f.write('inputs_reused:\n')
    for p in ['../../r1-adversarial/r1lib.py', '../../r1-adversarial/sa_r1.py', '../../r1_third_check_20260930.py', '../dlib.py', '../../route-freecount/pv_lib.py']:
        f.write(f'  - {{path: {p}, sha256: {sha(p)}}}\n')
    f.write('files:\n')
    for p in files:
        f.write(f'  - {{path: {p}, bytes: {os.path.getsize(p)}, sha256: {sha(p)}}}\n')
print(open('manifest.yaml').read()[:1500])
