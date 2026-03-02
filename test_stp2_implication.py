#!/usr/bin/env python3
"""
Test the implication:
  {LC(E) + LC(J) + J <= E coefficientwise + E ≽ J prefix} => STP2(E, J)?

All integer arithmetic.
"""

import random
import sys
import time

def is_lc(seq):
    """Check log-concavity: seq[k]^2 >= seq[k-1]*seq[k+1] for all valid k."""
    n = len(seq)
    for k in range(1, n - 1):
        if seq[k] * seq[k] < seq[k-1] * seq[k+1]:
            return False
    return True

def smallest_mode(seq):
    """Return index of the maximum value (smallest index if tie)."""
    mx = max(seq)
    for i, v in enumerate(seq):
        if v == mx:
            return i
    return 0

def check_prefix_dominance(E, J):
    """E ≽ J prefix: E[k+1]*J[k] >= E[k]*J[k+1] for k=0,...,mode(E)-1."""
    mode = smallest_mode(E)
    n = min(len(E), len(J))
    for k in range(min(mode, n - 1)):
        if k + 1 >= n:
            break
        if E[k+1] * J[k] < E[k] * J[k+1]:
            return False
    return True

def check_stp2(E, J):
    """STP2 diagonal: J[k+1]*E[k-1] <= J[k]*E[k] for k >= 1."""
    n = min(len(E), len(J))
    failures = []
    for k in range(1, n - 1):
        if J[k+1] * E[k-1] > J[k] * E[k]:
            failures.append(k)
    return failures

def gen_random_lc_seq(length, max_coeff):
    """Generate a random LC sequence with seq[0]=1."""
    seq = [1]
    seq.append(random.randint(2, max_coeff))
    for i in range(2, length):
        bound = (seq[i-1] * seq[i-1]) // seq[i-2]
        if bound <= 0:
            while len(seq) < length:
                seq.append(0)
            break
        val = random.randint(0, bound)
        seq.append(val)
        if val == 0:
            while len(seq) < length:
                seq.append(0)
            break
    return seq

def gen_j_basic(E):
    """Generate J <= E with J[0]=1, check LC."""
    n = len(E)
    J = [1]
    for k in range(1, n):
        J.append(random.randint(0, E[k]))
    return J if is_lc(J) else None

def gen_j_near_prefix(E):
    """Generate J that barely satisfies prefix dominance."""
    n = len(E)
    mode = smallest_mode(E)
    J = [1]
    for k in range(1, n):
        if k <= mode and E[k] > 0 and E[k-1] > 0:
            bound = min((E[k] * J[k-1]) // E[k-1], E[k])
            if bound <= 0:
                J.append(0)
            else:
                lo = max(0, bound - max(1, bound // 5))
                J.append(random.randint(lo, bound))
        elif E[k] > 0:
            J.append(random.randint(E[k] * 2 // 3, E[k]))
        else:
            J.append(0)
    return J if (len(J) == n and is_lc(J)) else None

def gen_j_slow_tail(E):
    """Generate J where J decays slower than E in tail."""
    n = len(E)
    mode = smallest_mode(E)
    J = [1]
    for k in range(1, n):
        if k <= mode and E[k] > 0 and E[k-1] > 0:
            bound = min((E[k] * J[k-1]) // E[k-1], E[k])
            if bound <= 0:
                J.append(0)
            else:
                J.append(random.randint(max(0, bound // 2), bound))
        elif E[k] > 0:
            J.append(random.randint(max(0, E[k] * 3 // 4), E[k]))
        else:
            J.append(0)
    return J if (len(J) == n and is_lc(J)) else None

def run_random_batch(num_trials, length_range, max_coeff_range, gen_j_func, label):
    """Run batch of random trials."""
    total = 0
    passed = 0
    stp2_fails = 0
    examples = []
    
    for _ in range(num_trials):
        length = random.randint(length_range[0], length_range[1])
        max_coeff = random.randint(max_coeff_range[0], max_coeff_range[1])
        E = gen_random_lc_seq(length, max_coeff)
        if len(E) < 3:
            continue
        J = gen_j_func(E)
        if J is None:
            continue
        total += 1
        if any(J[k] > E[k] for k in range(len(E))):
            continue
        if not check_prefix_dominance(E, J):
            continue
        passed += 1
        fails = check_stp2(E, J)
        if fails:
            stp2_fails += 1
            if len(examples) < 20:
                examples.append((list(E), list(J), fails))
    
    print(f"  [{label}] generated={total:,}  passed_all_4={passed:,}  STP2_fails={stp2_fails}", flush=True)
    return total, passed, stp2_fails, examples

def exhaustive_len4(max_a=30):
    """Exhaustive for length 4."""
    print(f"\n=== EXHAUSTIVE length=4, a in [1,{max_a}] ===", flush=True)
    total = 0
    passed = 0
    stp2_fails = 0
    examples = []
    
    for a in range(1, max_a + 1):
        for b in range(0, a * a + 1):
            if a * a < b:
                continue
            bound_c = (b * b) // a if a > 0 else 0
            for c in range(0, bound_c + 1):
                E = [1, a, b, c]
                mode_E = smallest_mode(E)
                for j1 in range(0, a + 1):
                    for j2 in range(0, b + 1):
                        if j1 * j1 < j2:
                            continue
                        for j3 in range(0, c + 1):
                            if j2 * j2 < j1 * j3:
                                continue
                            J = [1, j1, j2, j3]
                            total += 1
                            if not check_prefix_dominance(E, J):
                                continue
                            passed += 1
                            fails = check_stp2(E, J)
                            if fails:
                                stp2_fails += 1
                                if len(examples) < 10:
                                    examples.append((list(E), list(J), fails))
        if a % 5 == 0:
            print(f"    a={a}/{max_a}, pairs={total:,}, passed={passed:,}, fails={stp2_fails}", flush=True)
    
    print(f"  Total pairs: {total:,}  Passed: {passed:,}  STP2 failures: {stp2_fails}", flush=True)
    return total, passed, stp2_fails, examples

def exhaustive_len5(max_a=8):
    """Exhaustive for length 5, small coefficients."""
    print(f"\n=== EXHAUSTIVE length=5, a in [1,{max_a}] ===", flush=True)
    total = 0
    passed = 0
    stp2_fails = 0
    examples = []
    
    for a in range(1, max_a + 1):
        for b in range(0, a * a + 1):
            if a * a < b:
                continue
            bound_c = (b * b) // a if a > 0 else 0
            for c in range(0, bound_c + 1):
                if b > 0:
                    bound_d = (c * c) // b
                elif c == 0:
                    bound_d = 0
                else:
                    continue
                for d in range(0, bound_d + 1):
                    E = [1, a, b, c, d]
                    if not is_lc(E):
                        continue
                    mode_E = smallest_mode(E)
                    for j1 in range(0, a + 1):
                        for j2 in range(0, b + 1):
                            if j1 * j1 < j2:
                                continue
                            for j3 in range(0, c + 1):
                                if j2 > 0 and j2 * j2 < j1 * j3:
                                    continue
                                if j2 == 0 and j3 > 0:
                                    continue
                                for j4 in range(0, d + 1):
                                    if j3 > 0 and j3 * j3 < j2 * j4:
                                        continue
                                    if j3 == 0 and j4 > 0:
                                        continue
                                    J = [1, j1, j2, j3, j4]
                                    total += 1
                                    if not check_prefix_dominance(E, J):
                                        continue
                                    passed += 1
                                    fails = check_stp2(E, J)
                                    if fails:
                                        stp2_fails += 1
                                        if len(examples) < 10:
                                            examples.append((list(E), list(J), fails))
        print(f"    a={a}/{max_a}, pairs={total:,}, passed={passed:,}, fails={stp2_fails}", flush=True)
    
    print(f"  Total pairs: {total:,}  Passed: {passed:,}  STP2 failures: {stp2_fails}", flush=True)
    return total, passed, stp2_fails, examples

def main():
    random.seed(42)
    t0 = time.time()
    
    print("=" * 70)
    print("TESTING: LC(E) + LC(J) + J<=E + E≽J prefix  ==>  STP2(E,J)?")
    print("=" * 70, flush=True)
    
    all_examples = []
    grand = [0, 0, 0]
    
    # Phase 1: Random basic
    print("\n--- PHASE 1: Random basic J ---", flush=True)
    for num, lr, cr, label in [
        (3_000_000, (4, 8), (2, 50), "small-coeff"),
        (3_000_000, (4, 8), (50, 500), "med-coeff"),
        (2_000_000, (4, 10), (500, 10000), "large-coeff"),
        (2_000_000, (4, 6), (2, 100), "short-mixed"),
    ]:
        t, p, f, ex = run_random_batch(num, lr, cr, gen_j_basic, label)
        grand[0] += t; grand[1] += p; grand[2] += f
        all_examples.extend(ex)
    
    # Phase 2: Near-prefix
    print("\n--- PHASE 2: Near-prefix-boundary J ---", flush=True)
    for num, lr, cr, label in [
        (2_000_000, (4, 8), (2, 100), "near-pfx-small"),
        (2_000_000, (4, 8), (100, 1000), "near-pfx-med"),
        (1_000_000, (5, 10), (50, 5000), "near-pfx-long"),
    ]:
        t, p, f, ex = run_random_batch(num, lr, cr, gen_j_near_prefix, label)
        grand[0] += t; grand[1] += p; grand[2] += f
        all_examples.extend(ex)
    
    # Phase 3: Slow tail
    print("\n--- PHASE 3: Slow-tail J ---", flush=True)
    for num, lr, cr, label in [
        (2_000_000, (4, 8), (2, 100), "slow-tail-small"),
        (2_000_000, (5, 10), (50, 2000), "slow-tail-long"),
    ]:
        t, p, f, ex = run_random_batch(num, lr, cr, gen_j_slow_tail, label)
        grand[0] += t; grand[1] += p; grand[2] += f
        all_examples.extend(ex)
    
    # Phase 4: Exhaustive length 4
    t, p, f, ex = exhaustive_len4(max_a=30)
    grand[0] += t; grand[1] += p; grand[2] += f
    all_examples.extend(ex)
    
    # Phase 5: Exhaustive length 5
    t, p, f, ex = exhaustive_len5(max_a=8)
    grand[0] += t; grand[1] += p; grand[2] += f
    all_examples.extend(ex)
    
    elapsed = time.time() - t0
    
    print("\n" + "=" * 70)
    print("GRAND SUMMARY")
    print(f"  Total (E,J) pairs tested: {grand[0]:,}")
    print(f"  Passed all 4 conditions:  {grand[1]:,}")
    print(f"  STP2 failures:            {grand[2]}")
    print(f"  Elapsed:                  {elapsed:.1f}s")
    
    if all_examples:
        print(f"\n  COUNTEREXAMPLES ({len(all_examples)} found):")
        for E, J, fails in all_examples[:10]:
            print(f"    E = {E}")
            print(f"    J = {J}")
            mode_E = smallest_mode(E)
            print(f"    mode(E) = {mode_E}, STP2 fails at k = {fails}")
            for k in fails:
                lhs = J[k+1] * E[k-1]
                rhs = J[k] * E[k]
                print(f"      k={k}: J[{k+1}]*E[{k-1}] = {lhs} > J[{k}]*E[{k}] = {rhs}")
            # Also show prefix dominance check detail
            print(f"    Prefix dom check (k=0..{mode_E-1}):")
            for k in range(mode_E):
                if k + 1 < len(E):
                    l = E[k+1]*J[k]; r = E[k]*J[k+1]
                    print(f"      k={k}: E[{k+1}]*J[{k}]={l} >= E[{k}]*J[{k+1}]={r}? {l >= r}")
            print()
    else:
        print("\n  NO COUNTEREXAMPLES FOUND. Implication appears to hold.")
    print("=" * 70, flush=True)

if __name__ == "__main__":
    main()
