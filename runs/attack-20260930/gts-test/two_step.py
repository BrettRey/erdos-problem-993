"""Two-step repair of the weak form: for each tree whose every proper GTS image
has a larger central margin, is there a proper GTS image T' and a proper GTS
image T'' of T' with m(T'') <= m(T)? (Chaining two-step moves up the poset would
still reach the star.) Exact Fractions."""
import sys, json, networkx as nx
sys.argv = [sys.argv[0], '0']
exec(open('gts_margin_test.py').read().split("N = int(sys.argv[1])")[0])
res = {}
for n in range(7, 17):
    stuck = two_ok = 0; still = []
    for T in nx.nonisomorphic_trees(n):
        if max(d for _, d in T.degree()) == n - 1: continue
        mT = margin(T, 4); imgs = list(proper_gts_images(T))
        if any(margin(U, 4) <= mT for _, U in imgs): continue
        stuck += 1; ok = False
        for _, U in imgs:
            if max(d for _, d in U.degree()) == n - 1:   # U is the star: m(star) <= m(T) by census for n<=27
                ok = margin(U, 4) <= mT
            else:
                ok = any(margin(V, 4) <= mT for _, V in proper_gts_images(U))
            if ok: break
        if ok: two_ok += 1
        else: still.append(nx.to_graph6_bytes(T, header=False).decode().strip())
    res[n] = dict(stuck_one_step=stuck, repaired_two_step=two_ok, still_stuck=still)
    print(f'n={n} stuck(1-step)={stuck} repaired by 2 steps={two_ok} still stuck={len(still)} {still[:3]}', flush=True)
json.dump(res, open('two_step_n7_16.json', 'w'), indent=1)
