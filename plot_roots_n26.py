#!/usr/bin/env python3
"""
Plot complex roots of the independence polynomials for the two LC-failing
trees at n=26.

These are the only two trees (out of 279,793,450) at n=26 whose independence
polynomials fail log-concavity, matching the count from Kadrawi & Levit (2023).

Output: paper/figures/roots_n26_lc_failures.pdf
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# House style
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / ".house-style"))
from plot_style import setup, COLORS

setup()
plt.rcParams.update({
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsmath}',
})

# Independence polynomial coefficients [i_0, i_1, ..., i_alpha]
TREES = [
    {
        "label": "Tree 1",
        "lc_ratio": 1.145,
        "lc_index": 13,
        "coeffs": [1, 26, 300, 2040, 9142, 28551, 63933, 103736,
                    121376, 100144, 55499, 18683, 2979, 51, 1],
    },
    {
        "label": "Tree 2",
        "lc_ratio": 1.030,
        "lc_index": 13,
        "coeffs": [1, 26, 300, 2037, 9089, 28147, 62183, 98968,
                    112870, 90178, 48086, 15498, 2372, 48, 1],
    },
]


def find_roots(coeffs):
    """Find all roots of I(T; x) = sum_k coeffs[k] x^k."""
    return np.roots(list(reversed(coeffs)))


def classify_roots(roots, tol=1e-8):
    """Separate roots into real and complex (non-real) sets."""
    real_mask = np.abs(roots.imag) < tol
    real_roots = np.sort(roots[real_mask].real)
    complex_roots = roots[~real_mask]
    return real_roots, complex_roots


def plot_tree_roots(ax, info, clip_radius=1.4):
    """Plot roots in the complex plane, zoomed to show structure."""
    all_roots = find_roots(info["coeffs"])
    real_roots, complex_roots = classify_roots(all_roots)

    # Unit circle
    theta = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(theta), np.sin(theta),
            color='#D0D0D0', linewidth=0.7, zorder=1)

    # Crosshairs through origin
    ax.axhline(0, color='#D8D8D8', linewidth=0.4, zorder=0)
    ax.axvline(0, color='#D8D8D8', linewidth=0.4, zorder=0)

    # Roots inside the clipping radius
    inside_complex = complex_roots[np.abs(complex_roots) <= clip_radius]
    inside_real = real_roots[np.abs(real_roots) <= clip_radius]
    n_omitted = len(all_roots) - len(inside_complex) - len(inside_real)

    # Complex roots
    if len(inside_complex) > 0:
        ax.scatter(inside_complex.real, inside_complex.imag,
                   c=COLORS['primary'], s=20, zorder=3,
                   edgecolors='white', linewidths=0.3,
                   alpha=0.9, label="Complex")

    # Real roots
    if len(inside_real) > 0:
        ax.scatter(inside_real, np.zeros_like(inside_real),
                   c=COLORS['secondary'], s=24, zorder=4,
                   marker='D', edgecolors='white', linewidths=0.3,
                   alpha=0.9, label="Real")

    # Collect omitted root info (returned for figure-level note)
    omitted_info = None
    if n_omitted > 0:
        outside = all_roots[np.abs(all_roots) > clip_radius]
        unique_moduli = sorted(set(round(abs(r), 0) for r in outside))
        omitted_info = (n_omitted, unique_moduli)

    # Formatting
    margin = clip_radius * 1.05
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)
    ax.set_aspect('equal')
    ax.set_xlabel(r"$\operatorname{Re}$")
    ax.set_ylabel(r"$\operatorname{Im}$")

    # Title
    title = (f"{info['label']}  "
             f"($i_{{k-1}} i_{{k+1}} / i_k^2 = {info['lc_ratio']:.3f}$, "
             f"$k = {info['lc_index']}$)")
    ax.set_title(title, pad=8)

    ax.legend(loc='upper right', markerscale=0.9)

    # Summary to terminal
    moduli_all = np.abs(all_roots)
    print(f"{info['label']}: {len(real_roots)} real, "
          f"{len(complex_roots)} complex roots  "
          f"(moduli: {moduli_all.min():.3f} to {moduli_all.max():.1f})",
          flush=True)
    return omitted_info


def main():
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.8))

    omitted = []
    for i, info in enumerate(TREES):
        om = plot_tree_roots(axes[i], info)
        if om:
            omitted.append(om)

    fig.tight_layout(w_pad=2.0, rect=[0, 0.06, 1, 1])

    # Single footnote for omitted outlier roots
    if omitted:
        per_tree = []
        for n_om, mods in omitted:
            per_tree.append(str(int(mods[0])))
        mod_str = " and ".join(per_tree)
        note = (f"One conjugate pair omitted from each panel "
                rf"($|z| \approx {mod_str}$, respectively).")
        fig.text(0.5, 0.01, note, ha='center', fontsize=9,
                 color=COLORS['dark'])

    out_dir = Path(__file__).parent / "paper" / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "roots_n26_lc_failures"
    fig.savefig(out_path.with_suffix('.pdf'))
    print(f"Saved: {out_path.with_suffix('.pdf')}", flush=True)
    plt.close(fig)


if __name__ == "__main__":
    main()
