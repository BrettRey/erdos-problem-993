import sys, time
sys.path.insert(0,'scripts')
import vatter_lib_20261003 as V
fams = [
 ("star K_{1,40}", V.hubstar(40,0)),
 ("path P_60", V.path(60)),
 ("caterpillar 20x2", V.caterpillar(20,2)),
 ("spider S(3^15)", V.spider([3]*15)),
 ("H(9,2) n=28 [PV killer]", V.hubstar(9,2)),
 ("H(20,3)", V.hubstar(20,3)),
 ("MSH(4x9,4x10) n=237 [R1 killer]", V.msh([9]*4+[10]*4,2)),
 ("S(8,10,2) n=249", V.msh([10]*8,2)),
 ("S(9,9,2) n=253", V.msh([9]*9,2)),
 ("H(75,5) n=451 [2/3-rule killer]", V.hubstar(75,5)),
]
for name, adj in fams:
    t=time.time(); n=len(adj); a=V.alpha(adj); W=V.window(n,a)
    f=V.first_vertex_filter(adj, W)
    dead=[k for k in W if not f[k]]
    print(f"{name}: n={n} alpha={a} window=[{W.start},{W.stop-1}] k with no first vertex: {len(dead)} {dead[:8]}{'...' if len(dead)>8 else ''}  min #first-vertex candidates={min(len(f[k]) for k in W)}  ({time.time()-t:.1f}s)", flush=True)
