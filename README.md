# Erdős Problem #993

Computational search for a counterexample to the conjecture that the independent set sequence of every tree is unimodal.

## The problem

Given a tree *T* on *n* vertices, let *i_k(T)* = the number of independent sets of size *k* in *T*. The conjecture (Alavi, Malde, Schwenk, and Erdős 1987) states that the sequence *(i_0, i_1, ..., i_α)* is always unimodal: it increases to a peak and then decreases.

This is known to be false for general graphs. The tree case remains open.

- [erdosproblems.com #993](https://erdosproblems.com/993)

## Approach

1. Enumerate non-isomorphic trees on *n* vertices (using [nauty](https://pallini.di.uniroma1.it/)'s `geng` when available, otherwise custom generation)
2. Compute the independent set polynomial for each tree via dynamic programming on the rooted tree
3. Check whether the coefficient sequence is unimodal
4. Report any counterexample

## Setup

```bash
pip install networkx
```

Optional: install [nauty](https://pallini.di.uniroma1.it/) for fast tree enumeration (`brew install nauty` on macOS).

## License

MIT
