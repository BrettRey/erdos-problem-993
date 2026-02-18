import json
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

def get_path_poly(n):
    """Computes independence polynomial of Path graph P_n."""
    if n == 0: return [1]
    if n == 1: return [1, 1]
    
    # Recurrence: P_n = P_{n-1} + x * P_{n-2}
    p_prev = [1, 1] # P_1
    p_prev2 = [1]   # P_0
    
    for i in range(2, n + 1):
        # x * P_{n-2} shifts coefficients by 1
        term2 = [0] + p_prev2
        
        # Add P_{n-1} + term2
        len1 = len(p_prev)
        len2 = len(term2)
        new_poly = [0] * max(len1, len2)
        
        for k in range(len1):
            new_poly[k] += p_prev[k]
        for k in range(len2):
            new_poly[k] += term2[k]
            
        p_prev2 = p_prev
        p_prev = new_poly
        
    return p_prev

def analyze_roots(poly_coeffs):
    """
    Analyzes the roots of a polynomial.
    """
    # Reverse coefficients for numpy (highest power first)
    roots = np.roots(poly_coeffs[::-1])
    
    real_parts = roots.real
    imag_parts = roots.imag
    
    max_real = np.max(real_parts)
    min_real = np.min(real_parts)
    max_imag = np.max(np.abs(imag_parts))
    
    all_real = np.allclose(imag_parts, 0)
    
    return {
        "roots": roots,
        "max_real": max_real,
        "min_real": min_real,
        "max_imag": max_imag,
        "all_real": all_real
    }

def main():
    json_path = "results/analysis_n26.json"
    if len(sys.argv) > 1:
        json_path = sys.argv[1]
        
    if not os.path.exists(json_path):
        print(f"Error: File {json_path} not found.")
        return

    print(f"Analyzing roots for trees in: {json_path}")
    
    with open(json_path, 'r') as f:
        data = json.load(f)
        
    n = data.get("n", 26) # Default to 26 if not found
    near_misses = data.get("top_near_misses", [])
    lc_failures = data.get("lc_failures", [])
    
    print(f"N={n}")
    print(f"Found {len(near_misses)} near misses and {len(lc_failures)} LC failures.")
    
    plt.figure(figsize=(10, 10))
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Plot reference Path graph
    path_poly = get_path_poly(n)
    path_stats = analyze_roots(path_poly)
    plt.scatter(path_stats["roots"].real, path_stats["roots"].imag, 
                color='green', marker='o', s=50, label=f'Path P_{n} (Real Roots)', alpha=0.6)
    
    # Plot Near Misses
    print("\n--- Near Misses (Top 5) ---")
    for i, entry in enumerate(near_misses[:10]):
        poly = entry.get("poly")
        nm_ratio = entry.get("nm_ratio", 0)
        stats = analyze_roots(poly)
        
        label = f'NM #{i+1} ({nm_ratio:.4f})' if i < 3 else None # Label only top 3
        color = 'blue'
        # Vary alpha to show density
        plt.scatter(stats["roots"].real, stats["roots"].imag, 
                    color=color, marker='x', s=30, alpha=0.5, label=label)
        
        if i < 5:
            print(f"Rank {i+1} (NM={nm_ratio:.4f}): MaxRe={stats['max_real']:.4f}, MaxIm={stats['max_imag']:.4f}")

    # Plot LC Failures
    print("\n--- LC Failures ---")
    for i, entry in enumerate(lc_failures):
        poly = entry.get("poly")
        lc_ratio = entry.get("lc_ratio", 0)
        stats = analyze_roots(poly)
        
        label = f'LC Fail ({lc_ratio:.4f})'
        plt.scatter(stats["roots"].real, stats["roots"].imag, 
                    color='red', marker='^', s=60, label=label)
        
        print(f"LC Fail {i+1} (Ratio={lc_ratio:.4f}): MaxRe={stats['max_real']:.4f}, MaxIm={stats['max_imag']:.4f}")

    plt.title(f"Roots of Independence Polynomials (n={n})")
    plt.xlabel("Real Part")
    plt.ylabel("Imaginary Part")
    plt.legend()
    
    output_file = f"roots_n{n}.png"
    plt.savefig(output_file, dpi=150)
    print(f"\nPlot saved to {output_file}")

if __name__ == "__main__":
    main()
