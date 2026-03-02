#!/usr/bin/env python3
"""
Test: LC(E) + LC(J) + J<=E + E≽J prefix  ==>  STP2(E,J)?

v4: REQUIRE all entries of E and J to be STRICTLY POSITIVE (≥1).
This matches the tree DP setting where E and J are IS polynomials with positive coefficients.

Also test with J[k] ≥ 1 for k ≤ support(J).
"""

import random
import time
import sys

def is_lc(seq):
    for k in range(1, len(seq) - 1):
        if seq[k] * seq[k] < seq[k-1] * seq[k+1]:
            return False
    return True

def all_positive(seq):
    return all(x >= 1 for x in seq)

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
        if k + 1 >= n: break
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
    for k in fails:
        lhs = J[k+1] * E[k-1]
        rhs = J[k] * E[k]
        rel = "TAIL" if k > mode_E else ("PREFIX" if k < mode_E else "AT_MODE")
        print(f"        STP2 fail k={k} [{rel}]: J[{k+1}]*E[{k-1}]={J[k+1]}*{E[k-1]}={lhs}"
              f" > J[{k}]*E[{k}]={J[k]}*{E[k]}={rhs} (excess={lhs-rhs})")
    ratios = [f"{J[k]}/{E[k]}" for k in range(len(E))]
    print(f"        J/E: {ratios}")
    # Also show J[k]/E[k] as floats for quick inspection
    frats = [f"{J[k]/E[k]:.4f}" for k in range(len(E))]
    print(f"        J/E floats: {frats}")
    sys.stdout.flush()

def gen_lc_positive(length, max_start):
    """Generate LC sequence with all entries >= 1, seq[0]=1."""
    E = [1, random.randint(2, max_start)]
    for i in range(2, length):
        bound = (E[i-1] * E[i-1]) // E[i-2]
        if bound < 1:
            return E  # can't extend positively
        E.append(random.randint(1, bound))
    return E

def exhaustive_pos_len4(max_a=20):
    """Exhaustive for length 4, all entries >= 1."""
    print(f"=== EXHAUSTIVE POSITIVE len=4, a in [1,{max_a}] ===", flush=True)
    total = passed = stp2_fails = 0
    examples = []
    for a in range(1, max_a + 1):
        for b in range(1, a * a + 1):
            bc = (b * b) // a
            for c in range(1, bc + 1):
                E = [1, a, b, c]
                for j1 in range(1, a + 1):
                    for j2 in range(1, b + 1):
                        if j1 * j1 < j2: continue
                        for j3 in range(1, c + 1):
                            if j2 * j2 < j1 * j3: continue
                            J = [1, j1, j2, j3]
                            total += 1
                            if not check_prefix_dom(E, J): continue
                            passed += 1
                            fails = check_stp2(E, J)
                            if fails:
                                stp2_fails += 1
                                if len(examples) < 20:
                                    examples.append((list(E), list(J), fails))
        if a % 5 == 0:
            print(f"  a={a}/{max_a}, pairs={total:,}, passed={passed:,}, fails={stp2_fails}", flush=True)
    print(f"  TOTAL: pairs={total:,}  passed={passed:,}  STP2_fails={stp2_fails}", flush=True)
    return passed, stp2_fails, examples

def exhaustive_pos_len5(max_a=5):
    """Exhaustive for length 5, all entries >= 1."""
    print(f"\n=== EXHAUSTIVE POSITIVE len=5, a in [1,{max_a}] ===", flush=True)
    total = passed = stp2_fails = 0
    examples = []
    for a in range(1, max_a + 1):
        for b in range(1, a * a + 1):
            bc = (b * b) // a
            for c in range(1, bc + 1):
                bd = (c * c) // b if b > 0 else 0
                for d in range(1, bd + 1):
                    E = [1, a, b, c, d]
                    if not is_lc(E): continue
                    for j1 in range(1, a + 1):
                        for j2 in range(1, b + 1):
                            if j1*j1 < j2: continue
                            for j3 in range(1, c + 1):
                                if j2*j2 < j1*j3: continue
                                for j4 in range(1, d + 1):
                                    if j3*j3 < j2*j4: continue
                                    J = [1, j1, j2, j3, j4]
                                    total += 1
                                    if not check_prefix_dom(E, J): continue
                                    passed += 1
                                    fails = check_stp2(E, J)
                                    if fails:
                                        stp2_fails += 1
                                        if len(examples) < 20:
                                            examples.append((list(E), list(J), fails))
        print(f"  a={a}/{max_a}, pairs={total:,}, passed={passed:,}, fails={stp2_fails}", flush=True)
    print(f"  TOTAL: pairs={total:,}  passed={passed:,}  STP2_fails={stp2_fails}", flush=True)
    return passed, stp2_fails, examples

def random_positive_phase(num_trials):
    """Random trials with all-positive sequences."""
    print(f"\n=== RANDOM POSITIVE ({num_trials:,} trials) ===", flush=True)
    passed = stp2_fails = 0
    examples = []
    
    for trial in range(num_trials):
        length = random.randint(4, 10)
        max_start = random.choice([5, 10, 30, 100, 500, 2000])
        
        E = gen_lc_positive(length, max_start)
        if len(E) < 4: continue
        if not all_positive(E): continue
        
        mode_E = smallest_mode(E)
        strat = random.choice(["basic", "near_prefix", "slow_tail", "match_drop"])
        
        J = [1]
        for k in range(1, len(E)):
            ek = E[k]
            if ek < 1: break  # shouldn't happen
            if strat == "basic":
                J.append(random.randint(1, ek))
            elif strat == "near_prefix":
                if k <= mode_E and E[k-1] > 0:
                    bnd = min((ek * J[k-1]) // E[k-1], ek)
                    bnd = max(bnd, 1)
                    lo = max(1, bnd - max(1, bnd//5))
                    J.append(random.randint(lo, bnd))
                else:
                    J.append(random.randint(max(1, ek*2//3), ek))
            elif strat == "slow_tail":
                if k <= mode_E and E[k-1] > 0:
                    bnd = min((ek * J[k-1]) // E[k-1], ek)
                    bnd = max(bnd, 1)
                    J.append(random.randint(max(1, bnd//2), bnd))
                else:
                    J.append(random.randint(max(1, ek*3//4), ek))
            elif strat == "match_drop":
                if k <= mode_E:
                    J.append(random.randint(max(1, ek-3), ek))
                else:
                    J.append(random.randint(max(1, ek//2), ek))
        
        if len(J) != len(E): continue
        if not all_positive(J): continue
        if not is_lc(J): continue
        if any(J[k] > E[k] for k in range(len(E))): continue
        if not check_prefix_dom(E, J): continue
        
        passed += 1
        fails = check_stp2(E, J)
        if fails:
            stp2_fails += 1
            if len(examples) < 30:
                examples.append((list(E), list(J), fails))
        
        if trial % 2_000_000 == 0 and trial > 0:
            print(f"  ...{trial:,} trials, {passed:,} passed, {stp2_fails} STP2 fails", flush=True)
    
    print(f"  passed={passed:,}  STP2_fails={stp2_fails}", flush=True)
    return passed, stp2_fails, examples

def targeted_sharp_drop(num_trials):
    """E has sharp drop after mode, J has slow decay -- all positive."""
    print(f"\n=== TARGETED SHARP-DROP E ({num_trials:,} trials) ===", flush=True)
    passed = stp2_fails = 0
    examples = []
    
    for trial in range(num_trials):
        length = random.randint(5, 10)
        peak = random.randint(50, 5000)
        
        # Build E: ramp up to peak, then drop sharply
        mode_pos = random.randint(1, length - 2)
        E = [1]
        # Ramp up
        for k in range(1, mode_pos + 1):
            if k == 1:
                E.append(random.randint(2, peak))
            else:
                # LC: E[k] <= E[k-1]^2 / E[k-2], also want E[k] > E[k-1]
                bound = (E[k-1] * E[k-1]) // E[k-2]
                lo = E[k-1] + 1
                if lo > bound:
                    break  # can't increase further
                E.append(random.randint(lo, bound))
        
        if len(E) < mode_pos + 1:
            continue
        
        # Drop: fast
        for k in range(mode_pos + 1, length):
            bound = (E[k-1] * E[k-1]) // E[k-2]
            if bound < 1:
                break
            # Pick small to make sharp drop
            hi = max(1, min(bound, E[k-1] // 2))
            E.append(random.randint(1, hi))
        
        if len(E) < 4: continue
        if not is_lc(E): continue
        if not all_positive(E): continue
        
        mode_E = smallest_mode(E)
        
        # J: match ratio in prefix, slow decay in tail
        J = [1]
        for k in range(1, len(E)):
            if k <= mode_E and E[k-1] > 0:
                bnd = min((E[k] * J[k-1]) // E[k-1], E[k])
                bnd = max(bnd, 1)
                J.append(random.randint(max(1, bnd * 3 // 4), bnd))
            else:
                # Slow decay: J[k] close to E[k]
                J.append(random.randint(max(1, E[k] * 4 // 5), E[k]))
        
        if len(J) != len(E): continue
        if not all_positive(J): continue
        if not is_lc(J): continue
        if any(J[k] > E[k] for k in range(len(E))): continue
        if not check_prefix_dom(E, J): continue
        
        passed += 1
        fails = check_stp2(E, J)
        if fails:
            stp2_fails += 1
            if len(examples) < 20:
                examples.append((list(E), list(J), fails))
    
    print(f"  passed={passed:,}  STP2_fails={stp2_fails}", flush=True)
    return passed, stp2_fails, examples

def main():
    random.seed(42)
    t0 = time.time()
    
    print("=" * 70)
    print("TESTING (POSITIVE): LC(E) + LC(J) + J<=E + E>=J prefix => STP2?")
    print("All entries of E and J are >= 1 (strictly positive).")
    print("=" * 70, flush=True)
    
    all_examples = []
    grand_passed = 0
    grand_fails = 0
    
    # Exhaustive
    p, f, ex = exhaustive_pos_len4(max_a=12)
    grand_passed += p; grand_fails += f; all_examples.extend(ex)
    
    p, f, ex = exhaustive_pos_len5(max_a=5)
    grand_passed += p; grand_fails += f; all_examples.extend(ex)
    
    # Random positive
    p, f, ex = random_positive_phase(10_000_000)
    grand_passed += p; grand_fails += f; all_examples.extend(ex)
    
    # Targeted sharp drop
    p, f, ex = targeted_sharp_drop(5_000_000)
    grand_passed += p; grand_fails += f; all_examples.extend(ex)
    
    elapsed = time.time() - t0
    
    print("\n" + "=" * 70)
    print("GRAND SUMMARY (ALL POSITIVE)")
    print(f"  Passed all 4 conditions: {grand_passed:,}")
    print(f"  STP2 failures:           {grand_fails:,}")
    if grand_passed > 0:
        print(f"  Failure rate:            {grand_fails/grand_passed*100:.4f}%")
    print(f"  Elapsed:                 {elapsed:.1f}s")
    
    if all_examples:
        all_examples.sort(key=lambda x: (len(x[0]), sum(x[0])))
        unique = []
        seen = set()
        for E, J, fails in all_examples:
            key = (tuple(E), tuple(J))
            if key not in seen:
                seen.add(key)
                unique.append((E, J, fails))
        
        print(f"\n  COUNTEREXAMPLES ({len(unique)} unique, showing first 20):\n")
        for idx, (E, J, fails) in enumerate(unique[:20], 1):
            print_example(E, J, fails, idx)
        
        # Location analysis
        prefix_c = tail_c = mode_c = 0
        for E, J, fails in unique:
            mode = smallest_mode(E)
            for k in fails:
                if k < mode: prefix_c += 1
                elif k == mode: mode_c += 1
                else: tail_c += 1
        print(f"\n  FAILURE LOCATIONS (across {len(unique)} unique examples):")
        print(f"    In PREFIX (k < mode): {prefix_c}")
        print(f"    AT MODE  (k = mode):  {mode_c}")
        print(f"    In TAIL  (k > mode):  {tail_c}")
    else:
        print("\n  NO COUNTEREXAMPLES FOUND. Implication appears to hold for positive sequences.")
    
    print("=" * 70, flush=True)

if __name__ == "__main__":
    main()
