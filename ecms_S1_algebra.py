#!/usr/bin/env python3
"""Verify algebraic formula for S_1 and bound |remote| < 1/2.

KEY ALGEBRAIC RESULT (derived from cavity method):

For edge e = (u,v) with cavity messages A = R_{u→v}, B = R_{v→u}:
  For each w ∈ N(u)\{v}, with q_w = R_{w→u}:

    ΔP(w) = F_u × q_w/(1+q_w)

  where F_u = A(B²+B-1) / [(1+A+B)(1+AB)]

Similarly for w ∈ N(v)\{u}:
    ΔP(w) = F_v × R_{w→v}/(1+R_{w→v})

  where F_v = B(A²+A-1) / [(1+A+B)(1+AB)]

Therefore:
  S_1^u = F_u × Σ_{w ∈ N(u)\{v}} q_w/(1+q_w)
  S_1^v = F_v × Σ_{w ∈ N(v)\{u}} r_w/(1+r_w)

And: Σ q_w/(1+q_w) ≤ -log(A) (since ∏ 1/(1+q_w) = A and log(1+x) ≥ x/(1+x))

So: |S_1^u| ≤ |F_u| × (-log A), |S_1^v| ≤ |F_v| × (-log B)
And: |S_1| ≤ |S_1^u| + |S_1^v| ≤ |F_u|(-log A) + |F_v|(-log B)

We test whether this bound is always < 1/2.
"""

import json
import math
import os
import time
from collections import defaultdict

from trees import trees

MAX_N = 18


def cavity_messages(n, adj):
    if n <= 1:
        return {}
    parent = [-1] * n
    children = [[] for _ in range(n)]
    visited = [False] * n
    order = []
    visited[0] = True
    queue = [0]
    head = 0
    while head < len(queue):
        v = queue[head]; head += 1
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)
    msgs = {}
    for v in reversed(order):
        p = parent[v]
        if p == -1:
            continue
        prod = 1.0
        for c in children[v]:
            prod *= 1.0 / (1.0 + msgs[(c, v)])
        msgs[(v, p)] = prod
    for v in order:
        for c in children[v]:
            prod = 1.0
            if parent[v] != -1:
                prod *= 1.0 / (1.0 + msgs[(parent[v], v)])
            for c2 in children[v]:
                if c2 != c:
                    prod *= 1.0 / (1.0 + msgs[(c2, v)])
            msgs[(v, c)] = prod
    return msgs


def occupation_probs(n, adj, msgs=None):
    if msgs is None:
        msgs = cavity_messages(n, adj)
    P = [0.0] * n
    for v in range(n):
        Rv = 1.0
        for u in adj[v]:
            Rv *= 1.0 / (1.0 + msgs.get((u, v), 0.0))
        P[v] = Rv / (1.0 + Rv)
    return P


def contract_edge(n, adj, u, v):
    merged_neighbors = set()
    for w in adj[u]:
        if w != v:
            merged_neighbors.add(w)
    for w in adj[v]:
        if w != u:
            merged_neighbors.add(w)
    old_to_new = {}
    for i in range(n):
        if i == v:
            continue
        old_to_new[i] = i if i < v else i - 1
    n_new = n - 1
    adj_new = [[] for _ in range(n_new)]
    u_new = old_to_new[u]
    for w in sorted(merged_neighbors):
        w_new = old_to_new[w]
        adj_new[u_new].append(w_new)
        adj_new[w_new].append(u_new)
    for i in range(n):
        if i == u or i == v:
            continue
        i_new = old_to_new[i]
        for j in adj[i]:
            if j == u or j == v:
                continue
            j_new = old_to_new[j]
            if j_new not in adj_new[i_new]:
                adj_new[i_new].append(j_new)
    for i in range(n_new):
        adj_new[i].sort()
    return n_new, adj_new, old_to_new


def main():
    t0 = time.time()

    # TEST 1: Verify the algebraic formula for ΔP(w) at distance 1
    formula_errors = 0
    formula_max_err = 0.0

    # TEST 2: Verify S_1 = F_u × sum + F_v × sum
    S1_formula_errors = 0

    # TEST 3: Bound |S_1| ≤ |F_u|(-log A) + |F_v|(-log B) and check < 1/2
    max_S1_bound = 0.0
    max_S1_bound_info = {}

    # TEST 4: Actual max |S_1|
    max_actual_S1 = 0.0

    # TEST 5: Maximize g(A,B) = |F_u|(-log A) + |F_v|(-log B) over all edges
    max_g = 0.0
    max_g_info = {}

    # TEST 6: Check tighter bound using actual sum vs -log bound
    max_tighter_ratio = 0.0  # actual sum / (-log A)

    # TEST 7: remote as function of S_1 and higher-order terms
    max_remote_minus_S1 = 0.0  # |remote - S_1| -- the "tail" after dist 1

    total_edges = 0

    for n in range(3, MAX_N + 1):
        tn = time.time()
        n_edges = 0
        n_form_err = 0

        for _, adj in trees(n):
            msgs = cavity_messages(n, adj)
            P_T = occupation_probs(n, adj, msgs)

            seen = set()
            for u in range(n):
                for v in adj[u]:
                    e_key = (min(u, v), max(u, v))
                    if e_key in seen:
                        continue
                    seen.add(e_key)
                    total_edges += 1
                    n_edges += 1

                    A = msgs.get((u, v), 0.0)
                    B = msgs.get((v, u), 0.0)

                    # Factors
                    denom = (1 + A + B) * (1 + A * B)
                    F_u = A * (B * B + B - 1) / denom if denom > 0 else 0
                    F_v = B * (A * A + A - 1) / denom if denom > 0 else 0

                    # Contract and compute actual ΔP
                    nc, adjc, o2n = contract_edge(n, adj, u, v)
                    msgs_c = cavity_messages(nc, adjc)
                    P_Te = occupation_probs(nc, adjc, msgs_c)

                    # Verify formula for each w at distance 1
                    actual_S1_u = 0.0
                    formula_S1_u = 0.0
                    sum_Pw_u = 0.0

                    for w in adj[u]:
                        if w == v:
                            continue
                        q_w = msgs.get((w, u), 0.0)
                        Pw_cav = q_w / (1 + q_w)
                        predicted = F_u * Pw_cav
                        sum_Pw_u += Pw_cav

                        w_new = o2n[w]
                        actual = P_T[w] - P_Te[w_new]

                        err = abs(predicted - actual)
                        if err > 1e-10:
                            formula_errors += 1
                            n_form_err += 1
                        if err > formula_max_err:
                            formula_max_err = err

                        actual_S1_u += actual
                        formula_S1_u += predicted

                    actual_S1_v = 0.0
                    formula_S1_v = 0.0
                    sum_Pw_v = 0.0

                    for w in adj[v]:
                        if w == u:
                            continue
                        r_w = msgs.get((w, v), 0.0)
                        Pw_cav = r_w / (1 + r_w)
                        predicted = F_v * Pw_cav
                        sum_Pw_v += Pw_cav

                        w_new = o2n[w]
                        actual = P_T[w] - P_Te[w_new]

                        err = abs(predicted - actual)
                        if err > 1e-10:
                            formula_errors += 1
                            n_form_err += 1
                        if err > formula_max_err:
                            formula_max_err = err

                        actual_S1_v += actual
                        formula_S1_v += predicted

                    actual_S1 = actual_S1_u + actual_S1_v

                    # Verify S_1 formula
                    S1_err = abs((formula_S1_u + formula_S1_v) - actual_S1)
                    if S1_err > 1e-10:
                        S1_formula_errors += 1

                    if abs(actual_S1) > max_actual_S1:
                        max_actual_S1 = abs(actual_S1)

                    # TEST 3: Bound via -log
                    log_bound_u = abs(F_u) * (-math.log(A)) if A > 0 else 0
                    log_bound_v = abs(F_v) * (-math.log(B)) if B > 0 else 0
                    S1_bound = log_bound_u + log_bound_v

                    if S1_bound > max_S1_bound:
                        max_S1_bound = S1_bound
                        max_S1_bound_info = {
                            "n": n, "edge": (u, v),
                            "A": A, "B": B, "F_u": F_u, "F_v": F_v,
                            "sum_u": sum_Pw_u, "sum_v": sum_Pw_v,
                            "log_A": -math.log(A) if A > 0 else 0,
                            "log_B": -math.log(B) if B > 0 else 0,
                            "bound": S1_bound, "actual": actual_S1,
                        }

                    # TEST 5: g(A,B) -- the algebraic bound
                    g = S1_bound
                    if g > max_g:
                        max_g = g
                        max_g_info = {"A": A, "B": B, "g": g, "n": n}

                    # TEST 6: Tighter bound
                    if A > 0 and sum_Pw_u > 0:
                        r = sum_Pw_u / (-math.log(A))
                        if r > max_tighter_ratio:
                            max_tighter_ratio = r
                    if B > 0 and sum_Pw_v > 0:
                        r = sum_Pw_v / (-math.log(B))
                        if r > max_tighter_ratio:
                            max_tighter_ratio = r

                    # TEST 7: |remote - S_1|
                    dmu = sum(P_T) - sum(P_Te)
                    local = (A + B - A * B) / denom if denom > 0 else 0
                    remote = dmu - local
                    tail = abs(remote - actual_S1)
                    if tail > max_remote_minus_S1:
                        max_remote_minus_S1 = tail

        elapsed = time.time() - tn
        print(f"n={n}: {n_edges} edges, form_err={n_form_err}, "
              f"max_bound={max_S1_bound:.6f}, {elapsed:.1f}s", flush=True)

    total_time = time.time() - t0

    print(f"\n=== TEST 1: FORMULA VERIFICATION ===")
    print(f"ΔP(w) = F × q_w/(1+q_w) matches actual ΔP?")
    print(f"Mismatches (>1e-10): {formula_errors}")
    print(f"Max error: {formula_max_err:.2e}")

    print(f"\n=== TEST 2: S_1 FORMULA ===")
    print(f"S_1 = F_u×Σ + F_v×Σ mismatches: {S1_formula_errors}")

    print(f"\n=== TEST 3: |S_1| BOUND ===")
    print(f"|S_1| ≤ |F_u|(-log A) + |F_v|(-log B)")
    print(f"Max bound value: {max_S1_bound:.10f}")
    print(f"Max actual |S_1|: {max_actual_S1:.10f}")
    print(f"{'PASSES' if max_S1_bound < 0.5 else 'FAILS'}: bound < 1/2")
    print(f"Info: {max_S1_bound_info}")

    print(f"\n=== TEST 5: ALGEBRAIC BOUND g(A,B) ===")
    print(f"Max g: {max_g:.10f}")
    print(f"Info: {max_g_info}")

    print(f"\n=== TEST 6: SUM vs -LOG RATIO ===")
    print(f"Max (Σ P_w^cavity) / (-log A): {max_tighter_ratio:.6f}")
    print(f"(If always < 1, the -log bound is never achieved)")

    print(f"\n=== TEST 7: TAIL AFTER S_1 ===")
    print(f"Max |remote - S_1|: {max_remote_minus_S1:.8f}")
    print(f"(The higher-order contribution beyond distance 1)")

    # Numerical optimization of g(A,B) on (0,1]^2
    print(f"\n=== NUMERICAL MAX OF g(A,B) ===")
    best_g = 0
    best_AB = (0, 0)
    for i in range(1, 1000):
        A = i / 1000.0
        for j in range(1, 1000):
            B = j / 1000.0
            d = (1 + A + B) * (1 + A * B)
            Fu = abs(A * (B * B + B - 1)) / d
            Fv = abs(B * (A * A + A - 1)) / d
            g = Fu * (-math.log(A)) + Fv * (-math.log(B))
            if g > best_g:
                best_g = g
                best_AB = (A, B)
    print(f"Max g(A,B) on grid: {best_g:.10f}")
    print(f"Achieved at A={best_AB[0]:.3f}, B={best_AB[1]:.3f}")
    print(f"{'PASSES' if best_g < 0.5 else 'FAILS'}: max g < 1/2")

    results = {
        "max_n": MAX_N,
        "total_edges": total_edges,
        "formula_errors": formula_errors,
        "formula_max_err": formula_max_err,
        "S1_formula_errors": S1_formula_errors,
        "max_actual_S1": round(max_actual_S1, 10),
        "max_S1_bound": round(max_S1_bound, 10),
        "max_S1_bound_info": {k: (round(v, 8) if isinstance(v, float) else
                                   list(v) if isinstance(v, tuple) else v)
                               for k, v in max_S1_bound_info.items()},
        "max_g": round(max_g, 10),
        "max_tighter_ratio": round(max_tighter_ratio, 8),
        "max_remote_minus_S1": round(max_remote_minus_S1, 8),
        "numerical_max_g": round(best_g, 10),
        "numerical_max_AB": list(best_AB),
        "total_time_s": round(total_time, 2),
    }

    out_path = "results/ecms_S1_algebra.json"
    os.makedirs("results", exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
