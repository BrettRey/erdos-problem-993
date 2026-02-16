#!/usr/bin/env python3
"""Test: mode(∏R_c) <= mode(∏I_c) for all trees T and vertices v with children c.

This is the inductive step for proving δ(T,v) <= 1.

R_c = IS(T_c - {c}) and I_c = IS(T_c), so I_c >= R_c coefficientwise.
The claim is that the product of the "smaller" polynomials has mode
at most equal to the product of the "larger" ones.

Also test: mode(∏R_c) <= mode(∏I_c) - Σδ_c + 1 (tighter bound needed for proof).
"""
import subprocess
import sys
import time
from collections import Counter, deque

from indpoly import _polyadd, _polymul, independence_poly

GENG = "/opt/homebrew/bin/geng"


def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode("ascii")]
    n = data[0]
    idx = 1
    adj = [[] for _ in range(n)]
    bit = 5
    word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5
                idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj


def first_descent(seq):
    for k in range(len(seq) - 1):
        if seq[k] > seq[k + 1]:
            return k
    return len(seq) - 1


def tree_polys(adj, vertices, root):
    """Compute I(T) and R_root(T) = IS(T - {root}) for subtree."""
    vset = set(vertices)
    parent = {root: -1}
    order = [root]
    queue = deque([root])
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y in vset and y not in parent:
                parent[y] = x
                queue.append(y)
                order.append(y)
    dp_in = {}
    dp_out = {}
    for v in reversed(order):
        children = [y for y in adj[v] if y in vset and parent.get(y) == v]
        if not children:
            dp_in[v] = [0, 1]
            dp_out[v] = [1]
        else:
            prod_out = [1]
            for c in children:
                prod_out = _polymul(prod_out, dp_out[c])
            dp_in[v] = [0] + prod_out
            prod_both = [1]
            for c in children:
                prod_both = _polymul(prod_both, _polyadd(dp_in[c], dp_out[c]))
            dp_out[v] = prod_both
    I_T = _polyadd(dp_in[root], dp_out[root])
    R_root = dp_out[root]  # IS with root excluded
    return I_T, R_root, dp_in, dp_out, parent


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    print(f"Testing mode(∏R_c) <= mode(∏I_c) for tree vertices, n up to {max_n}", flush=True)
    print("=" * 90, flush=True)

    t0 = time.time()
    total = 0
    product_mode_fail = 0
    max_gap_product = 0  # mode(∏R_c) - mode(∏I_c)
    delta_ge0_fail = 0  # mode(I) < mode(R_v) (δ < 0)
    delta_dist = Counter()

    for n in range(3, max_n + 1):
        tn = time.time()
        n_total = 0
        n_fail = 0
        n_delta_neg = 0

        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)

            # For each vertex v, compute I(T_i) and R_{c_i} for children
            for v in range(nn):
                n_total += 1
                total += 1

                # Get tree structure rooted at v
                I_T, R_v, dp_in, dp_out, parent = tree_polys(adj, set(range(nn)), v)

                children = [y for y in adj[v] if parent.get(y) == v]
                if not children:
                    continue  # leaf, no children to check

                # Compute ∏I(T_i) and ∏R_{c_i}
                prod_I = [1]  # product of I(T_i)
                prod_R = [1]  # product of R_{c_i}
                sum_delta = 0

                for c in children:
                    # I(T_c) = dp_in[c] + dp_out[c]
                    I_c = _polyadd(dp_in[c], dp_out[c])
                    R_c = dp_out[c]  # IS with c excluded

                    m_Ic = first_descent(I_c)
                    m_Rc = first_descent(R_c)
                    delta_c = m_Ic - m_Rc
                    sum_delta += delta_c

                    prod_I = _polymul(prod_I, I_c)
                    prod_R = _polymul(prod_R, R_c)

                m_prodI = first_descent(prod_I)
                m_prodR = first_descent(prod_R)

                gap = m_prodR - m_prodI
                if gap > max_gap_product:
                    max_gap_product = gap
                if gap > 0:
                    n_fail += 1
                    product_mode_fail += 1
                    if product_mode_fail <= 5:
                        print(f"  FAIL: n={nn}, v={v}, deg={len(children)}, "
                              f"mode(∏R_c)={m_prodR}, mode(∏I_c)={m_prodI}, "
                              f"Σδ_c={sum_delta}", flush=True)

                # Also check δ(T, v)
                m_I = first_descent(I_T)
                m_Rv = first_descent(R_v)  # = m_prodI
                delta = m_I - m_Rv
                delta_dist[delta] += 1
                if delta < 0:
                    n_delta_neg += 1
                    delta_ge0_fail += 1

        proc.wait()
        elapsed = time.time() - tn
        print(f"n={n:2d}: verts={n_total:>10,}  prod_mode_fail={n_fail}  δ<0={n_delta_neg}  ({elapsed:.1f}s)",
              flush=True)

    total_time = time.time() - t0
    print("=" * 90, flush=True)
    print(f"Total: {total:,} vertex checks in {total_time:.1f}s", flush=True)
    print(flush=True)
    print(f"mode(∏R_c) > mode(∏I_c):  "
          f"{'NEVER' if product_mode_fail == 0 else f'{product_mode_fail}'}",
          flush=True)
    print(f"Max mode(∏R_c) - mode(∏I_c): {max_gap_product}", flush=True)
    print(f"δ(T,v) < 0: {'NEVER' if delta_ge0_fail == 0 else f'{delta_ge0_fail}'}", flush=True)
    print(flush=True)
    print("δ(T,v) distribution:", flush=True)
    for g in sorted(delta_dist.keys()):
        count = delta_dist[g]
        pct = 100 * count / total
        print(f"  {g:+d}: {count:>10,} ({pct:.2f}%)", flush=True)


if __name__ == "__main__":
    main()
