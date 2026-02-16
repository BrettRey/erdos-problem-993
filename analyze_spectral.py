#!/usr/bin/env python3
"""Spectral analysis of independence polynomials.

computes roots of the independence polynomial and analyzes their distribution
in the complex plane to identify "dangerous" trees that are close to breaking
unimodality.
"""

import argparse
import glob
import json
import math
import os
import sys
import numpy as np

# Add current directory to path to import indpoly if needed
sys.path.append(os.getcwd())
try:
    import indpoly
except ImportError:
    pass

def compute_spectral_metrics(poly_coeffs):
    """Compute spectral metrics from polynomial coefficients.
    
    Args:
        poly_coeffs: List of coefficients [i_0, i_1, ..., i_k]
        
    Returns:
        dict of metrics
    """
    # np.roots expects coefficients from highest degree to lowest
    # poly_coeffs is i_0 + i_1*x + ... (lowest to highest)
    roots = np.roots(poly_coeffs[::-1])
    
    if len(roots) == 0:
        return None

    # Filter out roots that are practically real (numerical noise)
    # But strictly, all roots of these polynomials should be considered.
    
    # Angles in (-pi, pi]
    angles = np.angle(roots)
    
    # We are interested in roots closest to the positive real axis (angle = 0)
    # Since coefficients are positive, there are no positive real roots.
    # Roots come in conjugate pairs. We look at absolute angle.
    abs_angles = np.abs(angles)
    
    min_angle = np.min(abs_angles)
    max_real = np.max(roots.real)
    max_imag = np.max(np.abs(roots.imag))
    
    # "Danger Score": How close is the root to the positive real axis?
    # A simple score is 1 / min_angle.
    danger_score = 1.0 / (min_angle + 1e-9)
    
    return {
        "roots": roots,
        "min_angle": min_angle,
        "min_angle_deg": np.degrees(min_angle),
        "max_real": max_real,
        "max_imag": max_imag,
        "danger_score": danger_score
    }

def find_polys(data, filepath, depth=0):
    """Recursively find polynomials in a JSON object."""
    if depth > 2: # Don't go too deep
        return

    if isinstance(data, dict):
        # Check for direct polynomial keys
        for key in data:
            if key == 'poly' or key.endswith('_poly'):
                val = data[key]
                if isinstance(val, list) and len(val) > 0 and isinstance(val[0], (int, float)):
                    # Found a polynomial!
                    yield val, data
        
        # Check for adj/tree structure if no poly found in this dict
        if 'adj' in data and 'n' in data:
             try:
                 poly = indpoly.independence_poly(data['n'], data['adj'])
                 yield poly, data
             except:
                 pass

        # Recurse
        for key, val in data.items():
            if isinstance(val, (dict, list)):
                yield from find_polys(val, filepath, depth + 1)
                
    elif isinstance(data, list):
        for item in data:
            if isinstance(item, (dict, list)):
                yield from find_polys(item, filepath, depth + 1)

def analyze_file(filepath):
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
    except Exception as e:
        # print(f"Error reading {filepath}: {e}", file=sys.stderr)
        return []

    results = []
    # Use the recursive finder
    for poly, context in find_polys(data, filepath):
        metrics = compute_spectral_metrics(poly)
        if metrics:
            metrics['filename'] = os.path.basename(filepath)
            metrics['n'] = context.get('n', len(poly)-1) if isinstance(context, dict) else len(poly)-1
            metrics['lc_ratio'] = context.get('lc_ratio', 0.0) if isinstance(context, dict) else 0.0
            metrics['unimodal'] = context.get('unimodal', True) if isinstance(context, dict) else True
            results.append(metrics)
            
    return results

def main():
    parser = argparse.ArgumentParser(description="Spectral analysis of independence polynomials")
    parser.add_argument('files', nargs='*', help="JSON files to analyze (default: results/*.json)")
    parser.add_argument('--top', type=int, default=20, help="Number of top results to show")
    parser.add_argument('--sort', choices=['danger', 'real', 'angle'], default='danger', help="Sort metric")
    
    args = parser.parse_args()
    
    files = args.files
    if not files:
        files = glob.glob("results/*.json")
        
    results = []
    print(f"Analyzing {len(files)} files...")
    
    for f in files:
        file_results = analyze_file(f)
        results.extend(file_results)
            
    # Sort
    if args.sort == 'danger':
        results.sort(key=lambda x: x['danger_score'], reverse=True)
    elif args.sort == 'real':
        results.sort(key=lambda x: x['max_real'], reverse=True)
    elif args.sort == 'angle':
        results.sort(key=lambda x: x['min_angle']) # Smallest angle = most dangerous
        
    print(f"\nTop {args.top} 'Dangerous' Trees (Spectral Analysis):")
    print(f"{'Filename':<55} | {'N':<4} | {'Min Angle':<9} | {'Max Real':<9} | {'Danger':<8} | {'LC Ratio':<8}")
    print("-" * 115)
    
    for r in results[:args.top]:
        print(f"{r['filename']:<55} | {r['n']:<4} | {r['min_angle_deg']:<8.2f}° | {r['max_real']:<9.4f} | {r['danger_score']:<8.2f} | {r['lc_ratio']:.4f}")

    # Also check if any found roots have positive real parts (shouldn't happen for stable polynomials, but check)
    positive_real_roots = [r for r in results if r['max_real'] > 0]
    if positive_real_roots:
        print("\nWARNING: Found trees with roots having positive real parts (Stability Violation?!):")
        for r in positive_real_roots:
             print(f"{r['filename']} (Max Re: {r['max_real']})")

if __name__ == "__main__":
    main()
