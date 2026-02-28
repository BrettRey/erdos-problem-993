"""Diagnose SV1 failures and Delta profile details from Round 4 scan."""

from scan_p2_round4 import (parse_g6, dp_rooted, dp_rooted_with_AB,
                              get_mode, coeff, count_sign_changes, sign_pattern,
                              check_p2, check_p3, _polyadd)
from fractions import Fraction
import subprocess


def diagnose_tree(g6):
    """Print detailed diagnostics for a single tree."""
    n, adj = parse_g6(g6)
    E0, J0 = dp_rooted(n, adj, 0)[:2]
    poly = _polyadd(E0, [0] + J0)
    m = get_mode(poly)
    print(f"\n{'='*60}")
    print(f"Tree: {g6.strip()}  n={n}  I(T)={poly}  mode={m}")

    # Leaf counts
    leaf_count = [0] * n
    for v in range(n):
        for u in adj[v]:
            if len(adj[u]) == 1:
                leaf_count[v] += 1
    print(f"Leaf counts: {leaf_count}")

    for r in range(n):
        if leaf_count[r] == 0:
            continue

        E, J, A, B, children, ell = dp_rooted_with_AB(n, adj, r)
        d_seq = []
        for k in range(m):
            dk = coeff(A, k+1) * coeff(B, k) - coeff(A, k) * coeff(B, k+1)
            d_seq.append(dk)
        changes = count_sign_changes(d_seq)

        # Delta profile
        delta_seq = []
        for k in range(1, m):
            bkm1 = coeff(B, k-1)
            if bkm1 == 0:
                delta_seq.append(('skip', k))
                continue
            dk = d_seq[k] if k < len(d_seq) else 0
            dkm1 = d_seq[k-1] if k-1 < len(d_seq) else 0
            bk = coeff(B, k)
            bkp1 = coeff(B, k+1)
            akm1 = coeff(A, k-1)
            ck = bk*bk - bkm1*bkp1
            delta_val = Fraction(bkm1*dk + bk*dkm1 + akm1*ck, bkm1)
            delta_seq.append((float(delta_val), k, dk, dkm1, ck))

        # Only print interesting ones (with sign changes or negative d_k)
        if changes > 1 or any(dk < 0 for dk in d_seq):
            p2_ok, _ = check_p2(E, J, m)
            print(f"\n  root={r} ell={ell} P2={'OK' if p2_ok else 'FAIL'}")
            print(f"  A={A[:m+2]}")
            print(f"  B={B[:m+2]}")
            print(f"  d_k={d_seq}  signs={sign_pattern(d_seq)} changes={changes}")
            if delta_seq:
                print(f"  Delta profile:")
                for item in delta_seq:
                    if item[0] == 'skip':
                        print(f"    k={item[1]}: skip (b_{{k-1}}=0)")
                    else:
                        val, k, dk, dkm1, ck = item
                        print(f"    k={k}: Delta={val:.4f}  d_k={dk}  d_{{k-1}}={dkm1}  c_k={ck}")


# Diagnose the first SV1 failure: n=10, g6=I???C@obg
print("=== SV1 FAILURE EXAMPLES ===")
diagnose_tree("I???C@obg")

# Also check n=14 for 3-change case
geng = '/opt/homebrew/bin/geng'
proc = subprocess.Popen([geng, '-q', '14', '13:13', '-c'],
                        stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    g6 = line.strip()
    n, adj = parse_g6(g6)
    E0, J0 = dp_rooted(n, adj, 0)[:2]
    poly = _polyadd(E0, [0] + J0)
    m = get_mode(poly)

    leaf_count = [0] * n
    for v in range(n):
        for u in adj[v]:
            if len(adj[u]) == 1:
                leaf_count[v] += 1

    for r in range(n):
        if leaf_count[r] == 0:
            continue
        E, J, A, B, children, ell = dp_rooted_with_AB(n, adj, r)
        d_seq = []
        for k in range(m):
            dk = coeff(A, k+1)*coeff(B, k) - coeff(A, k)*coeff(B, k+1)
            d_seq.append(dk)
        changes = count_sign_changes(d_seq)
        if changes >= 3:
            diagnose_tree(g6)
            break
    else:
        continue
    break  # found one

proc.terminate()

# Non-degenerate Delta minimum: look for trees where Delta_k > 0 but small
print("\n\n=== NON-DEGENERATE DELTA MINIMUM (n <= 14) ===")
delta_min = float('inf')
delta_min_info = None

for nn in range(5, 15):
    proc = subprocess.Popen([geng, '-q', str(nn), f'{nn-1}:{nn-1}', '-c'],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        g6 = line.strip()
        n, adj = parse_g6(g6)
        E0, J0 = dp_rooted(n, adj, 0)[:2]
        poly = _polyadd(E0, [0] + J0)
        m = get_mode(poly)

        leaf_count = [0] * n
        for v in range(n):
            for u in adj[v]:
                if len(adj[u]) == 1:
                    leaf_count[v] += 1

        for r in range(n):
            if leaf_count[r] == 0:
                continue
            E, J, A, B, children, ell = dp_rooted_with_AB(n, adj, r)
            for k in range(1, m):
                bkm1 = coeff(B, k-1)
                if bkm1 == 0:
                    continue
                dk = coeff(A,k+1)*coeff(B,k) - coeff(A,k)*coeff(B,k+1)
                dkm1 = coeff(A,k)*coeff(B,k-1) - coeff(A,k-1)*coeff(B,k)
                bk = coeff(B, k)
                bkp1 = coeff(B, k+1)
                akm1 = coeff(A, k-1)
                ck = bk*bk - bkm1*bkp1
                delta_num = bkm1*dk + bk*dkm1 + akm1*ck
                delta_val = Fraction(delta_num, bkm1)
                if 0 < delta_val < delta_min:
                    delta_min = delta_val
                    delta_min_info = (g6, n, r, k, ell, dk, dkm1, ck, float(delta_val))
    proc.wait()

if delta_min_info:
    g6, n, r, k, ell, dk, dkm1, ck, val = delta_min_info
    print(f"  Min nonzero Delta: {val:.6f} at {g6} n={n} root={r} k={k} ell={ell}")
    print(f"  d_k={dk}  d_{{k-1}}={dkm1}  c_k={ck}")
    diagnose_tree(g6)
