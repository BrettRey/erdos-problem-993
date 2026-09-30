"""Private-mass sufficient condition over the whole layered scan (depth 2..6)."""
import json, time
exec(open('private_mass.py').read().split("for name, b in")[0])
exec(open('scan_layered.py').read().split("grids = {")[0].split("t0 = time.time()")[0].replace("from layered import layered, coeffs", "").replace("def analyse", "def _unused"))
src = open('scan_layered.py').read(); g = src[src.index("grids = {"):src.index("for depth, grid in grids.items():")]
exec(g)
t0 = time.time(); tot = dict(trees=0, cases=0, fails=0); worst = None
for depth, grid in grids.items():
    for b in grid:
        n, t, f, w, dg = check(list(b)); tot['trees'] += 1; tot['cases'] += t; tot['fails'] += f
        if w is not None and (worst is None or w < worst[0]): worst = (w, list(b), n)
    print(f'depth {depth} done ({time.time()-t0:.0f}s): {tot}', flush=True)
print('min floor(private negative / positive) over all positive cases:', worst)
json.dump(dict(totals=tot, worst=worst), open('private_mass_scan.json', 'w'))
