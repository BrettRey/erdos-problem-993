#!/usr/bin/env python3
"""Plot independence polynomial root configurations for key trees.

Generates a 4-panel figure showing roots in the complex plane for:
  (a) n=26 LC-failure #1 (LC ratio 1.145)
  (b) n=26 LC-failure #2 (LC ratio 1.030)
  (c) Best near-miss tree at n=26 (nm=0.845)
  (d) Champion multi-arm star M(s; 5,5,4,2) at n=200
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from indpoly import independence_poly

# ── House style ──────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "figure.dpi": 200,
})

def poly_roots(coeffs):
    """Compute roots of polynomial given as coefficient list [a0, a1, ...]."""
    # numpy.roots expects highest degree first
    return np.roots(coeffs[::-1])

def build_multiarm_star(s, arms):
    """Build adjacency list for M(s; a1, a2, ..., ak)."""
    n = 1 + s + sum(arms)
    adj = [[] for _ in range(n)]
    hub = 0
    v = 1
    # pendant leaves
    for _ in range(s):
        adj[hub].append(v)
        adj[v].append(hub)
        v += 1
    # arms (paths from hub)
    for a in arms:
        prev = hub
        for _ in range(a):
            adj[prev].append(v)
            adj[v].append(prev)
            prev = v
            v += 1
    return n, adj

def plot_panel(ax, roots, title, highlight_real=True):
    """Plot roots in the complex plane on a given axis."""
    re = roots.real
    im = roots.imag
    ax.scatter(re, im, s=12, c="#2E86C1", alpha=0.7, edgecolors="k", linewidths=0.3, zorder=3)
    ax.axhline(0, color="#888", lw=0.5, zorder=1)
    ax.axvline(0, color="#888", lw=0.5, zorder=1)
    if highlight_real:
        real_roots = roots[np.abs(im) < 1e-10]
        if len(real_roots) > 0:
            ax.scatter(real_roots.real, [0]*len(real_roots), s=30, c="#E74C3C",
                       marker="D", edgecolors="k", linewidths=0.5, zorder=4)
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Re")
    ax.set_ylabel("Im")
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

def main():
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    fig.suptitle("Independence Polynomial Roots: Key Trees", fontsize=13, fontweight="bold")

    # ── Load n=26 analysis results ──
    with open("results/analysis_n26.json", "r") as f:
        data = json.load(f)

    # Panel (a): LC failure #1
    lc_trees = data.get("lc_failures", [])
    if len(lc_trees) >= 1:
        poly_a = lc_trees[0]["poly"]
        roots_a = poly_roots(poly_a)
        plot_panel(axes[0,0], roots_a, f"(a) LC failure #1, n=26\nLC ratio = {lc_trees[0].get('lc_ratio', 1.145):.3f}")
    else:
        axes[0,0].text(0.5, 0.5, "LC data not found", transform=axes[0,0].transAxes, ha="center")
        axes[0,0].set_title("(a) LC failure #1")

    # Panel (b): LC failure #2
    if len(lc_trees) >= 2:
        poly_b = lc_trees[1]["poly"]
        roots_b = poly_roots(poly_b)
        plot_panel(axes[0,1], roots_b, f"(b) LC failure #2, n=26\nLC ratio = {lc_trees[1].get('lc_ratio', 1.030):.3f}")
    else:
        axes[0,1].text(0.5, 0.5, "LC data not found", transform=axes[0,1].transAxes, ha="center")
        axes[0,1].set_title("(b) LC failure #2")

    # Panel (c): Best near-miss at n=26
    nm_trees = data.get("top_near_misses", [])
    if len(nm_trees) >= 1:
        poly_c = nm_trees[0]["poly"]
        nm_val = nm_trees[0].get("nm_ratio", 0.845)
        roots_c = poly_roots(poly_c)
        plot_panel(axes[1,0], roots_c, f"(c) Best near-miss, n=26\nnm = {nm_val:.4f}")
    else:
        axes[1,0].text(0.5, 0.5, "NM data not found", transform=axes[1,0].transAxes, ha="center")
        axes[1,0].set_title("(c) Best near-miss, n=26")

    # Panel (d): Champion multi-arm star M(s; 5,5,4,2) at n=200
    arms = [5, 5, 4, 2]
    s = 200 - 1 - sum(arms)  # s = 183
    n_d, adj_d = build_multiarm_star(s, arms)
    poly_d = independence_poly(n_d, adj_d)
    roots_d = poly_roots(poly_d)
    plot_panel(axes[1,1], roots_d, f"(d) $M({s};\\,5,5,4,2)$, n={n_d}")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = "paper/figures/root_configurations.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"Saved → {out}")

    out_png = "paper/figures/root_configurations.png"
    fig.savefig(out_png, bbox_inches="tight", dpi=200)
    print(f"Saved → {out_png}")

if __name__ == "__main__":
    main()
