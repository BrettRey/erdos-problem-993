"""No-dedup recheck of SoS(4,3), n=17: no chain of <= 2 proper GTS moves reaches
margin <= m(T); some chain of 3 does."""
import sys, time, networkx as nx
t0 = time.time()
sys.argv = [sys.argv[0], '0']
exec(open('gts_margin_test.py').read().split("N = int(sys.argv[1])")[0])
exec(open('star_of_stars_steps.py').read().split("for m in range(3, 7):")[0].split('exec(open')[0])
def sos(m, s):
    T = nx.Graph(); T.add_node(0); c = 1
    for _ in range(m):
        h = c; T.add_edge(0, h); c += 1
        for _ in range(s): T.add_edge(h, c); c += 1
    return T
T = sos(4, 3); mT = margin(T, 4)
L1 = [U for _, U in proper_gts_images(T)]
L2 = [V for U in L1 for _, V in proper_gts_images(U)]
print('images: 1-step', len(L1), '2-step', len(L2))
print('any 1-step <= m(T):', any(margin(U, 4) <= mT for U in L1))
print('any 2-step <= m(T):', any(margin(V, 4) <= mT for V in L2))
hit = next((W for V in L2 for _, W in proper_gts_images(V) if margin(W, 4) <= mT), None)
print('3-step hit:', hit is not None, 'leaves', sum(1 for v in hit if hit.degree(v) == 1) if hit else None, 'of max', T.number_of_nodes() - 1, f'({time.time()-t0:.1f}s)')
