#!/usr/bin/env python3
"""
Test: LC(E) + LC(J) + J<=E + E≽J prefix  ==>  STP2(E,J)?
v3: fast exhaustive (small range), big random phase.
"""

import random
import time
import sys

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
    for k in fails:
        lhs = J[k+1] * E[k-1]
        rhs = J[k] * E[k]
        rel = "TAIL" if k >= mode_E else ("PREFIX" if k < mode_E else "AT_MODE")
        print(f"        STP2 fail k={k} [{rel}]: J[{k+1}]*E[{k-1}]={J[k+1]}*{E[k-1]}={lhs}"
              f" > J[{k}]*E[{k}]={J[k]}*{E[k]}={rhs} (excess={lhs-rhs})")
    # J/E ratio
    ratios = [f"{J[k]}/{E[k]}" if E[k]>0 else "-" for k in range(len(E))]
    print(f"        J/E: {ratios}")
    sys.stdout.flush()

def exhaustive_len4(max_a=8):
    print(f"=== EXHAUSTIVE len=4, a in [1,{max_a}] ===", flush=True)
    total = passed = stp2_fails = 0
    examples = []
    for a in range(1, max_a + 1):
        for b in range(0, a * a + 1):
            bc = (b * b) // a if a > 0 else 0
            for c in range(0, bc + 1):
                E = [1, a, b, c]
                for j1 in range(0, a + 1):
                    for j2 in range(0, b + 1):
                        if j1 * j1 < j2: continue
                        for j3 in range(0, c + 1):
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
    print(f"  pairs={total:,}  passed={passed:,}  STP2_fails={stp2_fails}", flush=True)
    return passed, stp2_fails, examples

def exhaustive_len5(max_a=4):
    print(f"\n=== EXHAUSTIVE len=5, a in [1,{max_a}] ===", flush=True)
    total = passed = stp2_fails = 0
    examples = []
    for a in range(1, max_a + 1):
        for b in range(0, a * a + 1):
            bc = (b * b) // a if a > 0 else 0
            for c in range(0, bc + 1):
                bd = (c * c) // b if b > 0 else (0 if c == 0 else -1)
                if bd < 0: continue
                for d in range(0, bd + 1):
                    E = [1, a, b, c, d]
                    if not is_lc(E): continue
                    for j1 in range(a + 1):
                        for j2 in range(b + 1):
                            if j1*j1 < j2: continue
                            for j3 in range(c + 1):
                                if j2 > 0 and j2*j2 < j1*j3: continue
                                if j2 == 0 and j3 > 0: continue
                                for j4 in range(d + 1):
                                    if j3 > 0 and j3*j3 < j2*j4: continue
                                    if j3 == 0 and j4 > 0: continue
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

def random_phase(num_trials):
    print(f"\n=== RANDOM PHASE ({num_trials:,} trials) ===", flush=True)
    passed = stp2_fails = 0
    examples = []
    
    for trial in range(num_trials):
        length = random.randint(4, 10)
        max_coeff = random.choice([10, 30, 100, 500, 2000, 10000])
        
        E = [1, random.randint(2, max_coeff)]
        for i in range(2, length):
            bound = (E[i-1] * E[i-1]) // E[i-2]
            if bound <= 0:
                while len(E) < length: E.append(0)
                break
            E.append(random.randint(0, bound))
            if E[-1] == 0:
                while len(E) < length: E.append(0)
                break
        if len(E) < 4: continue
        
        mode_E = smallest_mode(E)
        strat = random.choice(["basic", "near_prefix", "slow_tail", "match_drop"])
        J = [1]
        ok = True
        for k in range(1, len(E)):
            if E[k] == 0:
                J.append(0); continue
            if strat == "basic":
                J.append(random.randint(0, E[k]))
            elif strat == "near_prefix":
                if k <= mode_E and E[k-1] > 0:
                    bnd = min((E[k] * J[k-1]) // E[k-1], E[k])
                    bnd = max(bnd, 0)
                    J.append(random.randint(max(0, bnd - max(1, bnd//5)), bnd))
                else:
                    J.append(random.randint(E[k]*2//3, E[k]))
            elif strat == "slow_tail":
                if k <= mode_E and E[k-1] > 0:
                    bnd = min((E[k] * J[k-1]) // E[k-1], E[k])
                    bnd = max(bnd, 0)
                    J.append(random.randint(max(0, bnd//2), bnd))
                else:
                    J.append(random.randint(max(0, E[k]*3//4), E[k]))
            elif strat == "match_drop":
                if k <= mode_E:
                    J.append(random.randint(max(0, E[k]-3), E[k]))
                else:
                    J.append(random.randint(E[k]//2, E[k]))
        
        if len(J) != len(E): continue
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

def main():
    random.seed(42)
    t0 = time.time()
    
    print("=" * 70)
    print("TESTING: LC(E) + LC(J) + J<=E + E≽J prefix  ==>  STP2(E,J)?")
    print("=" * 70, flush=True)
    
    all_examples = []
    grand_passed = 0
    grand_fails = 0
    
    # Exhaustive
    p, f, ex = exhaustive_len4(max_a=8)
    grand_passed += p; grand_fails += f; all_examples.extend(ex)
    
    p, f, ex = exhaustive_len5(max_a=4)
    grand_passed += p; grand_fails += f; all_examples.extend(ex)
    
    # Random
    p, f, ex = random_phase(10_000_000)
    grand_passed += p; grand_fails += f; all_examples.extend(ex)
    
    elapsed = time.time() - t0
    
    print("\n" + "=" * 70)
    print("GRAND SUMMARY")
    print(f"  Passed all 4 conditions: {grand_passed:,}")
    print(f"  STP2 failures:           {grand_fails:,}")
    if grand_passed > 0:
        print(f"  Failure rate:            {grand_fails/grand_passed*100:.3f}%")
    print(f"  Elapsed:                 {elapsed:.1f}s")
    
    if all_examples:
        # Sort by length, then sum of E
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
        print("\n  NO COUNTEREXAMPLES FOUND.")
    
    print("=" * 70, flush=True)

if __name__ == "__main__":
    main()
