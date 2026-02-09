# Erdős Problem #993 — Independent Set Sequence Unimodality for Trees

**Problem:** Is the independent set sequence of every tree unimodal?

**Source:** Alavi, Erdős, Malde, and Schwenk (1987)
**Reference:** [erdosproblems.com #993](https://erdosproblems.com/993)

## Statement

For a tree $T$ on $n$ vertices, let $i_k(T)$ = number of independent sets of size $k$. The sequence $(i_0, i_1, \ldots, i_\alpha)$ (where $\alpha$ is the independence number) is conjectured to be **unimodal**: it goes up, reaches a peak, then goes down. Never up-down-up.

Known: the conjecture is **false for general graphs** — only the tree case is open.

## Strategy

Computational search for a counterexample. Enumerate all non-isomorphic trees up to increasing vertex counts, compute independent set sequences, check unimodality.

## Status

- **Phase:** Setup
- **Created:** 2026-02-09
