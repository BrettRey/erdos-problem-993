#!/usr/bin/env python3
"""Experiment 2: Sensitivity analysis around the champion multi-arm star.

Alon's idea: perturb the champion tree M(183; 5,5,4,2) by a single edge
(prune-and-regraft one vertex) and measure how nm changes.
Is the champion on a smooth ridge or an isolated peak?
"""

from indpoly import independence_poly, near_miss_ratio as _nmr
import random
import copy

def nm(poly):
    """Get just the ratio from near_miss_ratio."""
    return _nmr(poly)[0]

def build_multiarm_star(s, arms):
    """Build adjacency list for M(s; a1, ..., ak)."""
    n = 1 + s + sum(arms)
    adj = [[] for _ in range(n)]
    hub = 0
    v = 1
    for _ in range(s):
        adj[hub].append(v)
        adj[v].append(hub)
        v += 1
    for a in arms:
        prev = hub
        for _ in range(a):
            adj[prev].append(v)
            adj[v].append(prev)
            prev = v
            v += 1
    return n, adj

def spr_perturbations(n, adj):
    """Generate all single-edge SPR (subtree prune-and-regraft) perturbations.
    
    For each non-bridge leaf, detach it and reattach to every other vertex.
    Yields (new_adj, description) pairs.
    """
    leaves = [v for v in range(n) if len(adj[v]) == 1]
    
    for leaf in leaves:
        parent = adj[leaf][0]
        # Detach leaf from parent, reattach to every other non-leaf vertex
        for target in range(n):
            if target == leaf or target == parent:
                continue
            # Build new adjacency
            new_adj = [list(nbrs) for nbrs in adj]
            new_adj[parent].remove(leaf)
            new_adj[leaf] = [target]
            new_adj[target].append(leaf)
            
            # Check it's still connected (simple check: leaf is connected to target)
            # Need to verify the whole tree is connected after removing leaf-parent edge
            # Since we started with a tree and moved one edge, it's still a tree iff
            # removing leaf from parent doesn't disconnect anything.
            # A leaf removal never disconnects, so this is always a valid tree.
            yield new_adj, f"move leaf {leaf} from {parent} to {target}"

def main():
    print("=" * 70)
    print("EXPERIMENT: Sensitivity of Champion Multi-Arm Star")
    print("=" * 70)
    print()
    
    # Build champion: M(s; 5,5,4,2), n=200
    arms = [5, 5, 4, 2]
    n = 200
    s = n - 1 - sum(arms)  # 183
    
    n_tree, adj = build_multiarm_star(s, arms)
    assert n_tree == n
    
    poly = independence_poly(n, adj)
    base_nm_200 = nm(poly)
    print(f"Champion: M({s}; {','.join(map(str,arms))}), n={n}")
    print(f"  Base nm = {base_nm_200:.6f}")
    print()
    
    # Too many perturbations for n=200. Sample.
    # Actually, let's do a smaller version first: M(s; 5,5,4,2) at n=50
    arms_small = [5, 5, 4, 2]
    n_small = 50
    s_small = n_small - 1 - sum(arms_small)  # 33
    
    n_tree_s, adj_s = build_multiarm_star(s_small, arms_small)
    poly_s = independence_poly(n_tree_s, adj_s)
    base_nm_s = nm(poly_s)
    
    print(f"Small champion: M({s_small}; 5,5,4,2), n={n_small}")
    print(f"  Base nm = {base_nm_s:.6f}")
    print()
    
    # Enumerate all leaf SPR perturbations
    better = []
    worse = []
    same = []
    count = 0
    
    for new_adj, desc in spr_perturbations(n_small, adj_s):
        try:
            new_poly = independence_poly(n_small, new_adj)
            new_nm = nm(new_poly)
            delta = new_nm - base_nm_s
            count += 1
            
            if delta > 0.0001:
                better.append((delta, desc, new_nm))
            elif delta < -0.0001:
                worse.append((delta, desc, new_nm))
            else:
                same.append((delta, desc, new_nm))
        except Exception as e:
            pass  # skip invalid trees
    
    print(f"  Total perturbations tested: {count}")
    print(f"  Better (Δnm > 0.0001): {len(better)}")
    print(f"  Worse  (Δnm < -0.0001): {len(worse)}")
    print(f"  Same   (|Δnm| ≤ 0.0001): {len(same)}")
    print()
    
    if better:
        better.sort(reverse=True)
        print("  Top 5 improvements:")
        for delta, desc, nm_val in better[:5]:
            print(f"    Δ = +{delta:.6f} (nm={nm_val:.6f}): {desc}")
        print()
    
    if worse:
        worse.sort()
        print("  Top 5 degradations:")
        for delta, desc, nm_val in worse[:5]:
            print(f"    Δ = {delta:.6f} (nm={nm_val:.6f}): {desc}")
        print()
    
    # Distribution of deltas
    all_deltas = [d for d, _, _ in better] + [d for d, _, _ in worse] + [d for d, _, _ in same]
    if all_deltas:
        print(f"  Δnm range: [{min(all_deltas):.6f}, {max(all_deltas):.6f}]")
        mean_delta = sum(all_deltas) / len(all_deltas)
        print(f"  Mean Δnm: {mean_delta:.6f}")
        print()
        
        if len(better) == 0:
            print("  ⟹ Champion is a LOCAL MAXIMUM: no single-leaf move improves nm.")
        elif len(better) < len(worse):
            print("  ⟹ Champion is near a peak: few improvements, many degradations.")
        else:
            print("  ⟹ Champion is NOT a local max: many improvements exist.")

    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
