#!/usr/bin/env python3
"""Deep investigation of the mean-mode relationship for IS distributions.

Questions:
  1. Is mode ≤ mean always (stronger than mode ≤ ceil(mean))?
  2. Does mode ≤ ceil(mean) fail for non-tree graphs?
  3. What are the skewness/variance patterns?
  4. How tight is the relationship? What's max(mode - mean)?
"""

import subprocess
import math
from indpoly import independence_poly
from graph6 import parse_graph6

def analyze_distribution(poly):
    """Compute mean, mode, variance, skewness of IS distribution."""
    total = sum(poly)
    if total == 0:
        return None
    
    # Normalized probabilities
    p = [c / total for c in poly]
    
    # Mean
    mean = sum(k * p[k] for k in range(len(p)))
    
    # Mode (index of max)
    mode = max(range(len(poly)), key=lambda k: poly[k])
    
    # Variance
    var = sum(k**2 * p[k] for k in range(len(p))) - mean**2
    
    # Skewness
    if var > 0:
        sd = math.sqrt(var)
        skew = sum((k - mean)**3 * p[k] for k in range(len(p))) / sd**3
    else:
        skew = 0.0
    
    return {
        'mean': mean,
        'mode': mode,
        'var': var,
        'sd': math.sqrt(var) if var > 0 else 0,
        'skew': skew,
        'mode_minus_mean': mode - mean,
        'alpha': len(poly) - 1,
        'n': None,  # set by caller
    }

def main():
    print("=" * 70)
    print("  DEEP INVESTIGATION: Mean vs Mode for Tree IS Distributions")
    print("=" * 70)
    print()
    
    # ── Part 1: Is mode ≤ mean (not just ceil)? ──
    print("PART 1: Is mode ≤ mean for all trees?")
    print()
    
    mode_gt_mean = 0
    max_excess = 0
    max_excess_info = None
    left_skewed = 0
    total = 0
    
    for n in range(3, 21):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_gt = 0
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            info = analyze_distribution(poly)
            info['n'] = tn
            total += 1
            
            if info['mode'] > info['mean']:
                mode_gt_mean += 1
                n_gt += 1
                excess = info['mode'] - info['mean']
                if excess > max_excess:
                    max_excess = excess
                    max_excess_info = info
                    max_excess_poly = poly
            
            if info['skew'] < -0.01:
                left_skewed += 1
        
        pct = n_gt / len(lines) * 100 if lines else 0
        print(f"  n={n:2d}: {n_gt:6d}/{len(lines):6d} trees have mode > mean ({pct:.1f}%)")
    
    print()
    print(f"  TOTAL: {mode_gt_mean}/{total} trees have mode > mean ({mode_gt_mean/total*100:.1f}%)")
    print(f"  Max(mode - mean) = {max_excess:.6f}")
    if max_excess_info:
        m = max_excess_info
        print(f"    At n={m['n']}: mode={m['mode']}, mean={m['mean']:.4f}, var={m['var']:.4f}, skew={m['skew']:.4f}")
        print(f"    Poly = {max_excess_poly}")
    print(f"  Left-skewed trees: {left_skewed}/{total}")
    print()
    
    # ── Part 2: Non-tree graphs ──
    print("PART 2: Does mode ≤ ceil(mean) fail for non-tree graphs?")
    print()
    print("  Testing all CONNECTED graphs (not just trees)...")
    
    nontree_violations = 0
    nontree_total = 0
    
    for n in range(3, 11):  # Keep small — lots of graphs
        cmd = f"/opt/homebrew/bin/geng {n} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        n_violations = 0
        for line in lines:
            tn, adj = parse_graph6(line)
            # Skip trees (n-1 edges)
            edge_count = sum(len(a) for a in adj) // 2
            if edge_count == tn - 1:
                continue
            
            poly = independence_poly(tn, adj)
            info = analyze_distribution(poly)
            nontree_total += 1
            
            if info['mode'] > math.ceil(info['mean']):
                n_violations += 1
                nontree_violations += 1
                if nontree_violations <= 5:
                    print(f"    VIOLATION: n={tn}, edges={edge_count}, mode={info['mode']}, "
                          f"mean={info['mean']:.4f}, ceil={math.ceil(info['mean'])}")
                    print(f"      Poly = {poly}")
        
        pct = n_violations / max(1, len(lines)) * 100
        print(f"  n={n}: checked {len(lines) - sum(1 for l in lines if sum(len(parse_graph6(l)[1][v]) for v in range(parse_graph6(l)[0])) // 2 == parse_graph6(l)[0] - 1)} non-tree graphs")
    
    print(f"\n  Total non-tree violations: {nontree_violations}/{nontree_total}")
    print()
    
    # ── Part 3: Variance and skewness patterns ──
    print("PART 3: Variance and skewness statistics")
    print()
    
    all_skews = []
    all_cv = []  # coefficient of variation = sd/mean
    
    for n in range(3, 19):
        cmd = f"/opt/homebrew/bin/geng {n} {n-1}:{n-1} -c -q"
        proc = subprocess.run(cmd, shell=True, capture_output=True)
        lines = [l for l in proc.stdout.split(b"\n") if l]
        
        skews_n = []
        cvs_n = []
        for line in lines:
            tn, adj = parse_graph6(line)
            poly = independence_poly(tn, adj)
            info = analyze_distribution(poly)
            skews_n.append(info['skew'])
            if info['mean'] > 0:
                cvs_n.append(info['sd'] / info['mean'])
            all_skews.append(info['skew'])
            all_cv.append(info['sd'] / info['mean'] if info['mean'] > 0 else 0)
        
        avg_skew = sum(skews_n) / len(skews_n)
        min_skew = min(skews_n)
        max_skew = max(skews_n)
        avg_cv = sum(cvs_n) / len(cvs_n) if cvs_n else 0
        
        print(f"  n={n:2d}: skew range [{min_skew:.4f}, {max_skew:.4f}], avg={avg_skew:.4f}; avg CV={avg_cv:.4f}")
    
    print()
    print(f"  Overall: {sum(1 for s in all_skews if s < 0)} left-skewed, "
          f"{sum(1 for s in all_skews if s >= 0)} right-skewed (out of {len(all_skews)})")
    
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
