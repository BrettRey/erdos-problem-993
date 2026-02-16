#!/usr/bin/env python3
"""Verify Hall's condition for heavy vertices in d_leaf ≤ 1 trees.

For each d_leaf ≤ 1 tree, compute H = {v : P(v) > 1/3} and check whether
H has an SDR (system of distinct representatives) in L = V \ H.

If Hall's condition holds for all such trees, then the edge bound
P(u)+P(v) < 2/3 immediately gives μ < n/3 (Conjecture A).
"""
import subprocess, sys, struct, itertools
from fractions import Fraction

def parse_graph6(s):
    s = s.strip()
    data = [c - 63 for c in s.encode('ascii')]
    n = data[0]; idx = 1
    adj = [[] for _ in range(n)]
    bit = 5; word = data[idx] if idx < len(data) else 0
    for j in range(n):
        for i in range(j):
            if word & (1 << bit):
                adj[i].append(j)
                adj[j].append(i)
            if bit == 0:
                bit = 5; idx += 1
                word = data[idx] if idx < len(data) else 0
            else:
                bit -= 1
    return n, adj

def is_dleaf_le1(n, adj):
    for v in range(n):
        leaf_children = sum(1 for u in adj[v] if len(adj[u]) == 1)
        if leaf_children > 1:
            return False
    return True

def compute_probs(n, adj):
    """Compute marginal P(v) using cavity method on tree."""
    if n == 1:
        return [Fraction(1, 2)]
    if n == 2:
        return [Fraction(1, 3), Fraction(1, 3)]

    # Root at vertex 0, BFS to get parent/children
    parent = [-1] * n
    children = [[] for _ in range(n)]
    order = []
    visited = [False] * n
    queue = [0]
    visited[0] = True
    while queue:
        v = queue.pop(0)
        order.append(v)
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                queue.append(u)

    # Bottom-up: compute R(v -> parent) = cavity message from v toward parent
    # R(v->parent) = 1 / prod_{c in children(v)} (1 + R(c->v))
    R_up = [Fraction(0)] * n  # R(v -> parent(v))
    for v in reversed(order):
        if not children[v]:  # leaf
            R_up[v] = Fraction(1)
        else:
            prod = Fraction(1)
            for c in children[v]:
                prod *= (1 + R_up[c])
            R_up[v] = Fraction(1, 1) / prod

    # Top-down: compute R(parent -> v) = cavity message from parent toward v
    # R(parent->v) = 1 / [(1 + R(parent->grandparent)) * prod_{sibling s} (1 + R(s->parent))]
    # But for the root, there's no parent message.
    R_down = [Fraction(0)] * n  # R(parent(v) -> v)

    for v in order:
        if v == 0:
            continue
        p = parent[v]
        # R(p -> v) = 1 / [(factor from p's parent) * prod_{sibling s of v} (1 + R(s->p))]
        # factor from p's parent: if p is root, no parent factor; else (1 + R_down[p])
        prod = Fraction(1)
        if parent[p] != -1:
            prod *= (1 + R_down[p])
        for s in children[p]:
            if s != v:
                prod *= (1 + R_up[s])
        R_down[v] = Fraction(1, 1) / prod if prod != 0 else Fraction(0)

    # Marginal probability P(v) = R_full(v) / (1 + R_full(v))
    # R_full(v) = 1 / [prod over all neighbors u of (1 + R(u->v))]
    # For root: R_full = 1 / prod_{c} (1 + R_up[c])
    # For non-root: R_full = 1 / [(1 + R_down[v]) * prod_{c in children(v)} (1 + R_up[c])]
    P = [Fraction(0)] * n
    for v in range(n):
        prod = Fraction(1)
        if v != 0 and parent[v] != -1:
            prod *= (1 + R_down[v])
        for c in children[v]:
            prod *= (1 + R_up[c])
        R_full = Fraction(1, 1) / prod if prod != 0 else Fraction(0)
        P[v] = R_full / (1 + R_full)

    return P

def check_hall(H, adj_H_to_L):
    """Check Hall's condition: for every S ⊆ H, |N(S)| ≥ |S|.

    adj_H_to_L: dict mapping h -> set of L-neighbors.
    Returns (True, None) or (False, violating_set).
    """
    H_list = list(H)
    # Check all subsets (exponential, but H is small for n ≤ 20)
    for size in range(1, len(H_list) + 1):
        for S in itertools.combinations(H_list, size):
            neighbors = set()
            for h in S:
                neighbors |= adj_H_to_L[h]
            if len(neighbors) < len(S):
                return False, set(S)
    return True, None

def find_max_matching_bipartite(H, L, adj):
    """Find maximum matching in bipartite graph H-L using augmenting paths."""
    match_H = {}  # h -> l
    match_L = {}  # l -> h

    def augment(h, visited):
        for l in adj[h]:
            if l in visited:
                continue
            visited.add(l)
            if l not in match_L or augment(match_L[l], visited):
                match_H[h] = l
                match_L[l] = h
                return True
        return False

    for h in H:
        augment(h, set())

    return match_H

def main():
    third = Fraction(1, 3)
    geng = "/opt/homebrew/bin/geng"

    print("Hall's condition check for H = {v : P(v) > 1/3} in d_leaf ≤ 1 trees")
    print("=" * 70)

    for n in range(3, 21):
        cmd = [geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        total_dleaf = 0
        total_with_H = 0
        hall_violations = 0
        max_H_size = 0
        sdr_failures = 0
        worst_tree = None

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            if not is_dleaf_le1(nn, adj):
                continue
            total_dleaf += 1

            P = compute_probs(nn, adj)
            H = {v for v in range(nn) if P[v] > third}

            if not H:
                continue
            total_with_H += 1
            max_H_size = max(max_H_size, len(H))

            L = set(range(nn)) - H
            adj_H_to_L = {h: {u for u in adj[h] if u in L} for h in H}

            # Quick check: maximum matching size
            match = find_max_matching_bipartite(H, L, adj_H_to_L)
            if len(match) < len(H):
                sdr_failures += 1
                worst_tree = (g6, nn, adj, H, L, P)

            # Full Hall check (only if H is small enough)
            if len(H) <= 15:
                ok, viol = check_hall(H, adj_H_to_L)
                if not ok:
                    hall_violations += 1

        proc.wait()
        print(f"n={n:2d}: d_leaf≤1={total_dleaf:7d}  with_H={total_with_H:6d}  "
              f"max|H|={max_H_size:2d}  SDR_fail={sdr_failures:4d}  "
              f"Hall_fail={hall_violations:4d}", flush=True)

        if worst_tree:
            g6, nn, adj, H, L, P = worst_tree
            print(f"  FAILURE: {g6}")
            print(f"  H={H}, |H|={len(H)}")
            print(f"  P = {[float(P[v]) for v in range(nn)]}")
            deg = [len(adj[v]) for v in range(nn)]
            print(f"  deg = {deg}")

if __name__ == "__main__":
    main()
