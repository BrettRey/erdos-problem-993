"""
Test ECMS mean shift |μ(T) - μ(T/e)| for parameterized tree families at large n.

For any tree T and edge e, ECMS conjectures |mode(I(T)) - mode(I(T/e))| <= 1.
Since both I(T) and I(T/e) are empirically log-concave, and for LC sequences
mode in {floor(mu), ceil(mu)}, we have |delta_mu| < 1 => ECMS.

This script tests specific families up to n=200, tracking max |delta_mu| over
all edges for each family member.
"""

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

# Add project dir to path for indpoly import
sys.path.insert(0, str(Path(__file__).parent))
from indpoly import independence_poly


def compute_mu_exact(coeffs):
    """Compute mu = I'(1)/I(1) = sum(k*c_k)/sum(c_k) using exact Fraction."""
    total = Fraction(0)
    weighted = Fraction(0)
    for k, c in enumerate(coeffs):
        total += c
        weighted += k * c
    if total == 0:
        return Fraction(0)
    return weighted / total


def contract_edge(n, adj, u, v):
    """Contract edge {u,v}: merge v into u, remove v, renumber vertices.

    Returns (n_new, adj_new, merged_vertex_new_index).
    """
    # Build new adjacency: merge v's neighbors (except u) into u
    # Then remove vertex v and renumber
    new_neighbors_u = set(adj[u]) - {v}
    for w in adj[v]:
        if w != u:
            new_neighbors_u.add(w)

    # Build adjacency with v removed
    # Renumber: vertices 0..n-1 minus v -> 0..n-2
    # Map old index to new index
    old_to_new = {}
    new_idx = 0
    for i in range(n):
        if i == v:
            continue
        old_to_new[i] = new_idx
        new_idx += 1

    n_new = n - 1
    adj_new = [[] for _ in range(n_new)]

    # For vertex u (now old_to_new[u]), its neighbors are new_neighbors_u
    u_new = old_to_new[u]
    for w in new_neighbors_u:
        w_new = old_to_new[w]
        adj_new[u_new].append(w_new)
        # Also update w's adjacency: remove old v, remove old u, add new u
        # We'll build all adjacency from scratch below

    # Actually, let's rebuild all adjacency cleanly
    adj_new = [[] for _ in range(n_new)]

    for i in range(n):
        if i == v:
            continue
        i_new = old_to_new[i]
        for j in adj[i]:
            if j == v:
                continue
            j_new = old_to_new[j]
            if j_new not in adj_new[i_new]:
                adj_new[i_new].append(j_new)

    # Now add edges from u to v's non-u neighbors
    for w in adj[v]:
        if w == u:
            continue
        w_new = old_to_new[w]
        if w_new not in adj_new[u_new]:
            adj_new[u_new].append(w_new)
        if u_new not in adj_new[w_new]:
            adj_new[w_new].append(u_new)

    return n_new, adj_new, u_new


def get_edges(n, adj):
    """Return list of unique edges (u, v) with u < v."""
    edges = []
    for u in range(n):
        for v in adj[u]:
            if u < v:
                edges.append((u, v))
    return edges


def max_delta_mu(n, adj):
    """Compute max |mu(T) - mu(T/e)| over all edges. Return details."""
    coeffs_T = independence_poly(n, adj)
    mu_T = compute_mu_exact(coeffs_T)

    edges = get_edges(n, adj)
    best_delta = Fraction(0)
    best_edge = None
    best_mu_Te = None

    for u, v in edges:
        n_c, adj_c, _ = contract_edge(n, adj, u, v)
        coeffs_Te = independence_poly(n_c, adj_c)
        mu_Te = compute_mu_exact(coeffs_Te)
        delta = abs(mu_T - mu_Te)
        if delta > best_delta:
            best_delta = delta
            best_edge = (u, v)
            best_mu_Te = mu_Te

    return {
        "n": n,
        "mu_T": mu_T,
        "max_delta_mu": best_delta,
        "best_edge": best_edge,
        "mu_Te": best_mu_Te,
    }


# ---- Tree Family Builders ----

def make_star(s):
    """Star K_{1,s}: hub=0 connected to s leaves 1..s. n = s+1."""
    n = s + 1
    adj = [[] for _ in range(n)]
    for i in range(1, n):
        adj[0].append(i)
        adj[i].append(0)
    return n, adj


def make_double_star(a, b):
    """Double star D(a,b): hub1=0 with a leaves, hub2=1 with b leaves.
    Edge between hub1 and hub2. n = a + b + 2."""
    n = a + b + 2
    adj = [[] for _ in range(n)]
    # Edge between hubs
    adj[0].append(1)
    adj[1].append(0)
    # hub1 leaves: vertices 2..a+1
    for i in range(2, a + 2):
        adj[0].append(i)
        adj[i].append(0)
    # hub2 leaves: vertices a+2..a+b+1
    for i in range(a + 2, n):
        adj[1].append(i)
        adj[i].append(1)
    return n, adj


def make_extended_star(s, k):
    """Extended star ES(s,k): star K_{1,s} where one leaf is replaced by
    a path of length k. Hub=0, leaves 1..s-1, then path 0-s-s+1-...-s+k-1.
    n = s + k."""
    n = s + k
    adj = [[] for _ in range(n)]
    # Hub=0 connected to s-1 leaves (vertices 1..s-1) plus start of path (vertex s)
    for i in range(1, s):
        adj[0].append(i)
        adj[i].append(0)
    # Path from hub: 0 - s - s+1 - ... - s+k-1
    if k >= 1:
        adj[0].append(s)
        adj[s].append(0)
    for i in range(s, s + k - 1):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    return n, adj


def make_subdivided_star(s, k):
    """Subdivided star SS(s,k): hub=0 with s arms, each arm is a path of
    length k. n = 1 + s*k."""
    n = 1 + s * k
    adj = [[] for _ in range(n)]
    for arm in range(s):
        # Arm vertices: 1 + arm*k, 1 + arm*k + 1, ..., 1 + arm*k + k - 1
        first = 1 + arm * k
        adj[0].append(first)
        adj[first].append(0)
        for j in range(k - 1):
            u = first + j
            v = first + j + 1
            adj[u].append(v)
            adj[v].append(u)
    return n, adj


def make_tripod_star(a, b, c):
    """Tripod star T(a,b,c): central vertex 0, connected to three bridge
    vertices 1,2,3. Bridge 1 has a pendant leaves, bridge 2 has b, bridge 3
    has c. n = 1 + 3 + a + b + c = 4 + a + b + c."""
    n = 4 + a + b + c
    adj = [[] for _ in range(n)]
    # Center = 0, bridges = 1, 2, 3
    for br in [1, 2, 3]:
        adj[0].append(br)
        adj[br].append(0)
    # Leaves for bridge 1: vertices 4..4+a-1
    idx = 4
    for i in range(a):
        adj[1].append(idx)
        adj[idx].append(1)
        idx += 1
    # Leaves for bridge 2
    for i in range(b):
        adj[2].append(idx)
        adj[idx].append(2)
        idx += 1
    # Leaves for bridge 3
    for i in range(c):
        adj[3].append(idx)
        adj[idx].append(3)
        idx += 1
    return n, adj


def label_edge(n, adj, u, v):
    """Describe an edge by vertex degrees for readability."""
    du = len(adj[u])
    dv = len(adj[v])
    return f"({u}[d={du}],{v}[d={dv}])"


def run_family(name, generator, params_list, results_all):
    """Run a family of trees, print and collect results."""
    print(f"\n{'='*70}")
    print(f"Family: {name}")
    print(f"{'='*70}")

    family_results = []
    family_max_delta = Fraction(0)
    family_max_info = None

    for params in params_list:
        n, adj = generator(*params)
        if n < 3:
            continue

        res = max_delta_mu(n, adj)
        delta_f = float(res["max_delta_mu"])

        edge_label = ""
        if res["best_edge"] is not None:
            edge_label = label_edge(n, adj, *res["best_edge"])

        entry = {
            "params": list(params),
            "n": res["n"],
            "mu_T": float(res["mu_T"]),
            "mu_Te": float(res["mu_Te"]) if res["mu_Te"] is not None else None,
            "max_delta_mu": delta_f,
            "best_edge": list(res["best_edge"]) if res["best_edge"] else None,
            "edge_label": edge_label,
        }
        family_results.append(entry)

        if res["max_delta_mu"] > family_max_delta:
            family_max_delta = res["max_delta_mu"]
            family_max_info = entry

        # Print progress for selected sizes
        if n <= 20 or n % 10 == 0 or delta_f > 0.5:
            print(
                f"  params={params}, n={n}: max|dmu|={delta_f:.6f}  "
                f"edge={edge_label}  mu_T={float(res['mu_T']):.6f}  "
                f"mu_Te={float(res['mu_Te']):.6f}" if res["mu_Te"] is not None else
                f"  params={params}, n={n}: max|dmu|={delta_f:.6f}"
            )

    print(f"\n  Family max |delta_mu| = {float(family_max_delta):.8f}")
    if family_max_info:
        print(f"  Achieved at params={family_max_info['params']}, n={family_max_info['n']}")
        print(f"  Edge: {family_max_info['edge_label']}")

    results_all[name] = {
        "family_max_delta_mu": float(family_max_delta),
        "family_max_info": family_max_info,
        "entries": family_results,
    }

    return family_max_delta


def main():
    t0 = time.time()
    results_all = {}
    overall_max = Fraction(0)
    overall_family = ""

    # --- Family a: Stars K_{1,s} for s=2..100 ---
    star_params = [(s,) for s in range(2, 101)]
    delta = run_family("Stars K_{1,s}", make_star, star_params, results_all)
    if delta > overall_max:
        overall_max = delta
        overall_family = "Stars K_{1,s}"

    # --- Family b: Double stars D(a,b) ---
    # Sweep: a=b from 1..49 (balanced), plus some unbalanced
    ds_params = []
    # Balanced: a=b
    for a in range(1, 50):
        ds_params.append((a, a))
    # Unbalanced: a=1, b up to 97
    for b in range(2, 98):
        ds_params.append((1, b))
    # Unbalanced: a=2, b up to 96
    for b in range(3, 50):
        ds_params.append((2, b))
    # Very unbalanced large
    for b in [50, 60, 70, 80, 90]:
        ds_params.append((1, b))
    delta = run_family("Double stars D(a,b)", make_double_star, ds_params, results_all)
    if delta > overall_max:
        overall_max = delta
        overall_family = "Double stars D(a,b)"

    # --- Family c: Extended stars ES(s,k) ---
    es_params = []
    for s in range(2, 51):
        for k in [1, 2, 3, 5, 10]:
            if s + k <= 100:
                es_params.append((s, k))
    # Also some with large k
    for k in range(1, 51):
        for s in [2, 3, 5, 10]:
            if s + k <= 100:
                es_params.append((s, k))
    # Deduplicate
    es_params = sorted(set(es_params))
    delta = run_family("Extended stars ES(s,k)", make_extended_star, es_params, results_all)
    if delta > overall_max:
        overall_max = delta
        overall_family = "Extended stars ES(s,k)"

    # --- Family d: Subdivided stars SS(s,k) ---
    ss_params = []
    for s in range(2, 31):
        for k in range(1, 31):
            if s * k <= 150:
                ss_params.append((s, k))
    ss_params = sorted(set(ss_params))
    delta = run_family("Subdivided stars SS(s,k)", make_subdivided_star, ss_params, results_all)
    if delta > overall_max:
        overall_max = delta
        overall_family = "Subdivided stars SS(s,k)"

    # --- Family e: Tripod stars T(a,b,c) ---
    tp_params = []
    # Balanced: a=b=c
    for a in range(1, 27):
        tp_params.append((a, a, a))
    # Two equal, one different
    for a in range(1, 30):
        for c in [1, 2, 5, 10]:
            if a + a + c <= 80:
                tp_params.append((a, a, c))
    # Very unbalanced
    for a in range(1, 40):
        tp_params.append((a, 1, 1))
    for a in range(1, 30):
        for b in range(1, 20):
            if a + b + 1 <= 80 and (a, b, 1) not in tp_params:
                tp_params.append((a, b, 1))
    tp_params = sorted(set(tp_params))
    delta = run_family("Tripod stars T(a,b,c)", make_tripod_star, tp_params, results_all)
    if delta > overall_max:
        overall_max = delta
        overall_family = "Tripod stars T(a,b,c)"

    # --- Summary ---
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"OVERALL SUMMARY")
    print(f"{'='*70}")
    print(f"Overall max |delta_mu| = {float(overall_max):.8f}")
    print(f"Achieved in family: {overall_family}")
    exceeds = float(overall_max) > 0.53
    print(f"Exceeds 0.53? {'YES' if exceeds else 'NO'}")
    print()

    # Trend for each family: max |delta_mu| at selected n values
    print("Trend of max |delta_mu| by family (last 5 entries by n):")
    for fname, fdata in results_all.items():
        entries = fdata["entries"]
        if not entries:
            continue
        # Sort by n, take last 5
        by_n = sorted(entries, key=lambda e: e["n"])
        tail = by_n[-5:]
        trend_str = "  ".join(
            f"n={e['n']}:{e['max_delta_mu']:.6f}" for e in tail
        )
        print(f"  {fname}: {trend_str}")

    print(f"\nTotal time: {elapsed:.1f}s")

    # Save results
    summary = {
        "overall_max_delta_mu": float(overall_max),
        "overall_family": overall_family,
        "exceeds_053": exceeds,
        "elapsed_seconds": elapsed,
        "families": {},
    }
    for fname, fdata in results_all.items():
        summary["families"][fname] = {
            "family_max_delta_mu": fdata["family_max_delta_mu"],
            "family_max_info": fdata["family_max_info"],
            "num_entries": len(fdata["entries"]),
            "entries": fdata["entries"],
        }

    outpath = Path(__file__).parent / "results" / "ecms_families.json"
    with open(outpath, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nResults saved to {outpath}")


if __name__ == "__main__":
    main()
