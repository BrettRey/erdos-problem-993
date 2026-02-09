# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Role: Researcher / Coder

This is a computational mathematics project. The goal is to search for a counterexample to Erdős Problem #993: that the independent set sequence of every tree is unimodal.

## The Problem

Given a tree $T$ on $n$ vertices, let $i_k(T)$ = number of independent sets of size $k$ in $T$. The conjecture states that the sequence $(i_0, i_1, \ldots, i_\alpha)$ is unimodal (non-decreasing then non-increasing).

A **counterexample** is any tree where this sequence is not unimodal (i.e., it dips then rises again).

## Approach

1. Enumerate non-isomorphic trees (use nauty/geng if available, otherwise custom generation)
2. For each tree, compute the independent set polynomial (count independent sets by size)
3. Check whether the coefficient sequence is unimodal
4. Report any non-unimodal sequence immediately

## Source Grounding

- Primary source: Alavi, Malde, Schwenk, and Erdős (1987)
- Problem page: erdosproblems.com/993
- Do NOT fabricate results or claim verification without running actual computations

## Setup

Python 3.14 is available. No dependencies are pre-installed. Before running anything:

```bash
pip install networkx   # graph operations
# Optional:
pip install numpy      # faster array operations
pip install sympy      # symbolic polynomial manipulation
```

nauty (`geng`) is **not installed**. To use it for tree enumeration:

```bash
brew install nauty
# Then: geng -c <n> | filterg -e<n-1>  # connected graphs with exactly n-1 edges = trees
# Or:   geng <n> <n-1>:<n-1>  -c       # equivalent: connected, exactly n-1 edges
```

Without nauty, generate trees in Python (e.g., Prüfer sequences or recursive construction).

## Key Algorithms

### Independent set polynomial via dynamic programming

For a tree, root it at any vertex, then DP bottom-up:
- `dp[v][0]` = polynomial counting independent sets in subtree(v) where v is excluded
- `dp[v][1]` = polynomial counting independent sets in subtree(v) where v is included

```
dp[v][0] = product over children c of (dp[c][0] + dp[c][1])
dp[v][1] = x * product over children c of dp[c][0]
```

The total polynomial is `dp[root][0] + dp[root][1]`. Coefficients give $i_k$.

This runs in $O(n^2)$ time per tree (polynomial multiplication at each node). For large $n$, the polynomial multiplication dominates.

### Unimodality check

A sequence $(a_0, a_1, \ldots, a_m)$ is unimodal if there exists a peak index $p$ such that $a_0 \leq a_1 \leq \cdots \leq a_p \geq a_{p+1} \geq \cdots \geq a_m$. Equivalently: no index $i$ with $a_{i-1} > a_i < a_{i+1}$ (no valley after a descent).

### Tree counts (OEIS A000055)

| $n$ | Trees |
|-----|-------|
| 1-10 | 1, 1, 1, 2, 3, 6, 11, 23, 47, 106 |
| 15 | 7741 |
| 20 | 823,065 |
| 25 | 105,157,672 |

Exhaustive enumeration becomes expensive around $n = 20$-25. Beyond that, targeted search (e.g., caterpillars, spiders, double stars) or heuristic exploration may be needed.

## Running

```bash
# Run the main search (once implemented)
python3 search.py

# Run tests
python3 -m pytest tests/

# Run a single test
python3 -m pytest tests/test_file.py::test_name
```

## Performance Notes

- Pure Python polynomial multiplication is the bottleneck. Use `numpy.polymul` or represent polynomials as lists and multiply with convolution.
- For $n > 18$, consider C extensions or rewriting the core DP in C/Rust.
- Parallelism: tree enumeration and checking are embarrassingly parallel. Use `multiprocessing.Pool`.
- Rust and Cargo are **not installed** on this system. If needed: `brew install rust`.
