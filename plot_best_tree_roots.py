
import json
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def analyze_roots(poly_coeffs):
    roots = np.roots(poly_coeffs[::-1])
    return roots

def main():
    json_path = "best_roots_tree_n30.json"
    if len(sys.argv) > 1:
        json_path = sys.argv[1]
        
    with open(json_path, 'r') as f:
        data = json.load(f)
        
    poly = data.get("best_poly")
    if not poly:
        poly = data.get("poly")
        
    if not poly:
        print("No polynomial found in file.")
        return
        
    roots = analyze_roots(poly)
    
    plt.figure(figsize=(10, 10))
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(True, linestyle='--', alpha=0.3)
    
    plt.scatter(roots.real, roots.imag, color='blue', marker='x', s=50, label='Roots')
    
    nm = data.get("best_nm", data.get("nm", 0))
    plt.title(f"Roots of Best Near-Miss Tree (n=30, NM={nm:.4f})")
    plt.xlabel("Real Part")
    plt.ylabel("Imaginary Part")
    plt.legend()
    
    out_file = "roots_best_n30.png"
    plt.savefig(out_file, dpi=150)
    print(f"Plot saved to {out_file}")
    
    print("Root Stats:")
    print(f"Max Real: {np.max(roots.real):.4f}")
    print(f"Min Real: {np.min(roots.real):.4f}")
    print(f"Max Imag: {np.max(np.abs(roots.imag)):.4f}")

if __name__ == "__main__":
    main()
