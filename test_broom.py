from indpoly import independence_poly, is_log_concave, log_concavity_ratio
from targeted import make_broom

n, adj = make_broom(10, 20)
poly = independence_poly(n, adj)
lc = is_log_concave(poly)
viol, pos = log_concavity_ratio(poly)
print(f'Broom(10,20) n={n}: LC={lc}, max_viol={viol:.4f} at pos {pos}')
