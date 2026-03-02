#!/usr/bin/env python3
"""
Test: LC(E) + LC(J) + J<=E + E≽J prefix  ==>  STP2(E,J)?
Result from v1: MASSIVELY FALSE. This version focuses on:
1. Quick exhaustive for tiny sequences (get smallest counterexamples)
2. Detailed analysis of counterexample structure
3. Summary statistics
"""

import random
import time

def is_lc(seq):
    for k in range(1, len(seq) - 1):
        if seq[k] * seq[k] < seq[k-1] * seq[k+1]:
            return False
    return True

def smallest_mode(seq):
    mx = max(seq)
    for i, v in enumerate(seq):
        if v == mx:
            return i
    return 0

def check_prefix_dom(E, J):
    mode = smallest_mode(E)
    n = min(len(E), len(J))
    for k in range(min(mode, n - 1)):
        if k + 1 >= n:
            break
        if E[k+1] * J[k] < E[k] * J[k+1]:
            return False
    return True

def check_stp2(E, J):
    n = min(len(E), len(J))
    failures = []
    for k in range(1, n - 1):
        if J[k+1] * E[k-1] > J[k] * E[k]:
            failures.append(k)
    return failures

def print_example(E, J, fails, idx):
    mode_E = smallest_mode(E)
    print(f"  #{idx}: E = {E}")
    print(f"        J = {J}")
    print(f"        mode(E) = {mode_E}, len = {len(E)}")
    print(f"        STP2 fails at k = {fails}")
    for k in fails:
        lhs = J[k+1] * E[k-1]
        rhs = J[k] * E[k]
        print(f"          k={k}: J[{k+1}]*E[{k-1}] = {J[k+1]}*{E[k-1]} = {lhs}"
              f"  vs  J[{k}]*E[{k}] = {J[k]}*{E[k]} = {rhs}"
              f"  (excess = {lhs - rhs})")
    # Show ratio J[k]/E[k]
    ratios = []
    for k in range(len(E)):
        if E[k] > 0:
            ratios.append(f"{J[k]}/{E[k]}")
        else:
            ratios.append("0/0")
    print(f"        J/E ratios: {ratios}")
    # Where does the failure sit relative to mode?
    for k in fails:
        rel = "TAIL" if k >= mode_E else ("PREFIX" if k < mode_E else "AT MODE")
        print(f"          k={k} is in the {rel} (mode={mode_E})")
    print()

def exhaustive_len4(max_a=15):
    """Exhaustive for length 4."""
    print(f"=== EXHAUSTIVE length=4, a in [1,{max_a}] ===", flush=True)
    total = 0
    passed = 0
    stp2_fails = 0
    examples = []
    
    for a in range(1, max_a + 1):
        for b in range(0, a * a + 1):
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
                            if not check_prefix_dom(E, J):
                                continue
                            passed += 1
                            fails = check_stp2(E, J)
                            if fails:
                                stp2_fails += 1
                                if len(examples) < 5:
                                    examples.append((list(E), list(J), fails))
        if a % 5 == 0:
            print(f"  a={a}/{max_a}, pairs={total:,}, passed={passed:,}, fails={stp2_fails}", flush=True)
    
    print(f"  TOTAL: pairs={total:,}, passed={passed:,}, STP2_failures={stp2_fails}", flush=True)
    return total, passed, stp2_fails, examples

def exhaustive_len5(max_a=6):
    """Exhaustive for length 5."""
    print(f"\n=== EXHAUSTIVE length=5, a in [1,{max_a}] ===", flush=True)
    total = 0
    passed = 0
    stp2_fails = 0
    examples = []
    
    for a in range(1, max_a + 1):
        for b in range(0, a * a + 1):
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
                                    if not check_prefix_dom(E, J):
                                        continue
                                    passed += 1
                                    fails = check_stp2(E, J)
                                    if fails:
                                        stp2_fails += 1
                                        if len(examples) < 5:
                                            examples.append((list(E), list(J), fails))
        print(f"  a={a}/{max_a}, pairs={total:,}, passed={passed:,}, fails={stp2_fails}", flush=True)
    
    print(f"  TOTAL: pairs={total:,}, passed={passed:,}, STP2_failures={stp2_fails}", flush=True)
    return total, passed, stp2_fails, examples

def random_phase(num_trials, label):
    """Run diverse random trials."""
    print(f"\n=== RANDOM: {label} ({num_trials:,} trials) ===", flush=True)
    total = 0
    passed = 0
    stp2_fails = 0
    examples = []
    
    for _ in range(num_trials):
        length = random.randint(4, 10)
        max_coeff = random.choice([20, 50, 200, 1000, 5000])
        
        # Generate E
        E = [1, random.randint(2, max_coeff)]
        for i in range(2, length):
            bound = (E[i-1] * E[i-1]) // E[i-2]
            if bound <= 0:
                while len(E) < length:
                    E.append(0)
                break
            E.append(random.randint(0, bound))
            if E[-1] == 0:
                while len(E) < length:
                    E.append(0)
                break
        if len(E) < 4:
            continue
        
        mode_E = smallest_mode(E)
        
        # Generate J with different strategies
        strategy = random.choice(["basic", "near_prefix", "slow_tail", "match_then_drop"])
        J = [1]
        for k in range(1, len(E)):
            if E[k] == 0:
                J.append(0)
                continue
            
            if strategy == "basic":
                J.append(random.randint(0, E[k]))
            elif strategy == "near_prefix":
                if k <= mode_E and E[k-1] > 0:
                    bnd = min((E[k] * J[k-1]) // E[k-1], E[k])
                    J.append(random.randint(max(0, bnd - max(1, bnd//5)), max(0, bnd)))
                else:
                    J.append(random.randint(E[k]*2//3, E[k]))
            elif strategy == "slow_tail":
                if k <= mode_E and E[k-1] > 0:
                    bnd = min((E[k] * J[k-1]) // E[k-1], E[k])
                    J.append(random.randint(max(0, bnd//2), max(0, bnd)))
                else:
                    J.append(random.randint(max(0, E[k]*3//4), E[k]))
            elif strategy == "match_then_drop":
                if k <= mode_E:
                    # nearly match E
                    J.append(random.randint(max(0, E[k]-2), E[k]))
                else:
                    # drop but not too fast
                    J.append(random.randint(E[k]//2, E[k]))
        
        if len(J) != len(E):
            continue
        if not is_lc(J):
            continue
        if any(J[k] > E[k] for k in range(len(E))):
            continue
        if not check_prefix_dom(E, J):
            continue
        
        total += 1
        passed += 1
        fails = check_stp2(E, J)
        if fails:
            stp2_fails += 1
            if len(examples) < 10:
                examples.append((list(E), list(J), fails))
    
    print(f"  passed_all_4={passed:,}, STP2_failures={stp2_fails}", flush=True)
    return total, passed, stp2_fails, examples

def main():
    random.seed(42)
    t0 = time.time()
    
    print("=" * 70)
    print("TESTING: LC(E) + LC(J) + J<=E + E≽J prefix  ==>  STP2(E,J)?")
    print("=" * 70, flush=True)
    
    all_examples = []
    grand = [0, 0, 0]
    
    # Exhaustive small
    for func, args in [(exhaustive_len4, (15,)), (exhaustive_len5, (6,))]:
        t, p, f, ex = func(*args)
        grand[0] += t; grand[1] += p; grand[2] += f
        all_examples.extend(ex)
    
    # Random phases
    for num, label in [
        (5_000_000, "diverse mixed strategies"),
        (5_000_000, "diverse round 2"),
    ]:
        t, p, f, ex = random_phase(num, label)
        grand[0] += t; grand[1] += p; grand[2] += f
        all_examples.extend(ex)
    
    elapsed = time.time() - t0
    
    print("\n" + "=" * 70)
    print("GRAND SUMMARY")
    print(f"  Total (E,J) pairs tested (all 4 conds): {grand[1]:,}")
    print(f"  STP2 failures:                          {grand[2]:,}")
    if grand[1] > 0:
        print(f"  Failure rate:                           {grand[2]/grand[1]*100:.2f}%")
    print(f"  Elapsed:                                {elapsed:.1f}s")
    
    if all_examples:
        # Sort by length then by sum of E
        all_examples.sort(key=lambda x: (len(x[0]), sum(x[0])))
        print(f"\n  SMALLEST COUNTEREXAMPLES ({min(len(all_examples), 15)} shown):\n")
        for idx, (E, J, fails) in enumerate(all_examples[:15]):
            print_example(E, J, fails, idx + 1)
    else:
        print("\n  NO COUNTEREXAMPLES FOUND.")
    
    # Additional analysis: where do failures concentrate?
    if all_examples:
        print("  FAILURE LOCATION ANALYSIS:")
        prefix_count = 0
        tail_count = 0
        at_mode_count = 0
        for E, J, fails in all_examples:
            mode = smallest_mode(E)
            for k in fails:
                if k < mode:
                    prefix_count += 1
                elif k == mode:
                    at_mode_count += 1
                else:
                    tail_count += 1
        print(f"    Failures in PREFIX (k < mode): {prefix_count}")
        print(f"    Failures AT MODE (k = mode):   {at_mode_count}")
        print(f"    Failures in TAIL (k > mode):   {tail_count}")
    
    print("=" * 70, flush=True)

if __name__ == "__main__":
    main()
