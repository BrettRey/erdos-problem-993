#!/usr/bin/env python3
"""Attempt algebraic proof of the combined tail condition.

From A <= (1+x)I:
  A_{d+2} <= I_{d+2} + I_{d+1}
  A_{d+1} <= I_{d+1} + I_d

Combined tail at d+1: (I+A)_{d+2} <= (I+A)_{d+1}
<=> I_{d+2} + A_{d+2} <= I_{d+1} + A_{d+1}
<=> A_{d+2} - A_{d+1} <= I_{d+1} - I_{d+2}

From A <= (1+x)I:
  A_{d+2} - A_{d+1} <= (I_{d+2} + I_{d+1}) - A_{d+1}

So we need: (I_{d+2} + I_{d+1}) - A_{d+1} <= I_{d+1} - I_{d+2}
<=> 2*I_{d+2} <= A_{d+1}

Question: is A_{d+1} >= 2*I_{d+2} always?

Also check: is I_d + I_{d+1} >= 2*I_{d+2} always? (from LC of I)
By LC: I_{d+1}^2 >= I_d * I_{d+2}, so I_{d+2} <= I_{d+1}^2/I_d.
I_d + I_{d+1} >= 2*I_{d+1}^2/I_d <=> I_d^2 + I_d*I_{d+1} >= 2*I_{d+1}^2
<=> (I_d - I_{d+1})(I_d + 2*I_{d+1}) >= 0.
Since I_d > I_{d+1} (descending at d), this is TRUE.

So I_d + I_{d+1} >= 2*I_{d+2}. But A_{d+1} <= I_d + I_{d+1} is an UPPER bound.
We need A_{d+1} >= 2*I_{d+2} which is a LOWER bound.
The upper bound exceeds the lower bound, but that doesn't prove A_{d+1} >= 2*I_{d+2}.

Let's check computationally whether A_{d+1} >= 2*I_{d+2} always holds in gap cases.
"""
import subprocess
from collections import deque
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


def split_at_edge(n, adj, u, v):
    A = set()
    queue = deque([u])
    A.add(u)
    while queue:
        x = queue.popleft()
        for y in adj[x]:
            if y not in A and y != v:
                A.add(y)
                queue.append(y)
    return A, set(range(n)) - A


def rooted_is_poly(adj, vertices, root):
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
    return dp_in[root], dp_out[root]


def main():
    print("Testing A_{d+1} >= 2*I_{d+2} in gap cases (d_A >= d+2)", flush=True)
    print("=" * 70, flush=True)

    total = 0
    gap_total = 0
    bound_fail = 0  # A_{d+1} < 2*I_{d+2}
    tightest = None

    # Also test: is A_{d+1} >= I_{d+2} always? (weaker bound)
    weak_fail = 0

    # Even weaker: is (I+A)_{d+1} >= 2*I_{d+2} + A_{d+2}? (direct combined tail)
    # This is just the combined tail itself, verified 0 failures above.

    # New bound attempt: can we bound A_{d+1} from below using structure?
    # A_{d+1} = (P_u*P_v)_{d+1} + (R_u*R_v)_d
    # Lower bound: A_{d+1} >= (R_u*R_v)_d since both terms nonneg
    # Is (R_u*R_v)_d >= 2*I_{d+2}?
    rr_bound_fail = 0

    for n in range(3, 20):
        cmd = [GENG, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        n_gap = 0
        n_fail = 0

        for line in proc.stdout:
            g6 = line.decode().strip()
            nn, adj = parse_graph6(g6)
            I_T = independence_poly(nn, adj)
            d = first_descent(I_T)

            edges = [(u, v) for u in range(nn) for v in adj[u] if u < v]
            for u, v in edges:
                total += 1
                side_u, side_v = split_at_edge(nn, adj, u, v)
                P_u, R_u = rooted_is_poly(adj, side_u, u)
                P_v, R_v = rooted_is_poly(adj, side_v, v)

                A = _polyadd(_polymul(P_u, P_v), [0] + _polymul(R_u, R_v))
                d_A = first_descent(A)

                if d_A >= d + 2:
                    n_gap += 1
                    gap_total += 1

                    A_d1 = A[d + 1] if d + 1 < len(A) else 0
                    I_d2 = I_T[d + 2] if d + 2 < len(I_T) else 0

                    if A_d1 < 2 * I_d2:
                        n_fail += 1
                        bound_fail += 1
                        ratio = A_d1 / (2 * I_d2) if I_d2 > 0 else float('inf')
                        if tightest is None or ratio < tightest[0]:
                            tightest = (ratio, g6, nn, u, v, d, A_d1, I_d2)

                    if A_d1 < I_d2:
                        weak_fail += 1

                    # Check RR bound
                    RR = _polymul(R_u, R_v)
                    rr_d = RR[d] if d < len(RR) else 0
                    if rr_d < 2 * I_d2:
                        rr_bound_fail += 1

        proc.wait()
        print(f"n={n:2d}: gap_cases={n_gap:5d}  A[d+1]<2I[d+2]={n_fail}", flush=True)

    print()
    print(f"Total gap cases: {gap_total}")
    print(f"  A[d+1] < 2*I[d+2]: {bound_fail} ({100*bound_fail/max(gap_total,1):.1f}%)")
    print(f"  A[d+1] < I[d+2]:   {weak_fail}")
    print(f"  (RR)_d < 2*I[d+2]: {rr_bound_fail}")

    if tightest:
        ratio, g6, nn, u, v, d, A_d1, I_d2 = tightest
        print(f"\n  Tightest: A[d+1]/2I[d+2] = {ratio:.4f}")
        print(f"    n={nn} edge=({u},{v}) d={d} A[d+1]={A_d1} I[d+2]={I_d2}")
        print(f"    tree={g6}")


if __name__ == "__main__":
    main()
