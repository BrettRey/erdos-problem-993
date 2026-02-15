# Proof of Private Neighbor Property for n < 3k Gap

## Problem Statement
We wish to prove the **Private Neighbor Property (PNP)**:
For any tree $T$ and any maximal independent set $S$ of size $k < \text{mode}(I(T))$, there exists $u \in S$ with $\text{priv}(u) \ge 2$.

The case $n \ge 3k$ is proved via the Private Neighbor Bound ($P \ge n - 2k + 1$) and pigeonhole principle.
The gap is $n < 3k$.

## Analytic Proof Step 1: PNP Failure Implies P ≤ k

If PNP fails for maximal IS $S$ of size $k$, then $\text{priv}(u) \le 1$ for all $u \in S$.
Therefore $P = \sum_{u \in S} \text{priv}(u) \le k$.

**Important correction (2026-02-14):** An earlier version cited Bollobás & Cockayne (1979) to claim $\text{priv}(u) \ge 1$ for all $u$ in any maximal IS (external private neighbors). This is **FALSE**. In K_{1,4}, the maximal IS $S = \{$all 4 leaves$\}$ has $\text{priv}(u) = 0$ for every leaf, because the center (the only non-$S$ vertex) has 4 $S$-neighbors.

The Bollobás-Cockayne result is about *closed* private neighbors (including the vertex itself) for *minimal dominating sets*. For independent dominating sets, $u$ is always its own closed private neighbor ($N[u] \cap S = \{u\}$), which is trivially true but doesn't give external private neighbors.

Empirically (`verify_external_priv.py`, n ≤ 18): many maximal IS below mode have vertices with $\text{priv} = 0$ (e.g., 18,262 such IS at n=17). But PNP still holds in all cases because some *other* vertex has $\text{priv} \ge 2$.

**The proof does not need $\text{priv}(u) \ge 1$ for all $u$.** It only needs $P \le k$, which follows directly from $\text{priv}(u) \le 1$ for all $u$.

## Analytic Proof Step 2: The Gap Connection

If PNP fails, then $P \le k$. Using the Private Neighbor Bound $P \ge n - 2k + 1$:
$$k \ge P \ge n - 2k + 1 \implies 3k \ge n + 1 \implies k \ge \frac{n+1}{3}$$

Since $k$ is an integer: $k \ge \lceil (n+1)/3 \rceil = \lfloor n/3 \rfloor + 1$.

This confirms that **PNP failure implies we are in the n < 3k gap**, and gives a lower bound on $k$.

## Analytic Proof Step 3: The Mode Conjecture

**1-Private Mode Conjecture (generalized):**
> For any tree $T$, if $S$ is a maximal IS with $\text{priv}(u) \le 1$ for all $u \in S$, then $|S| \ge \text{mode}(I(T))$.

This is equivalent to: PNP failure implies $k \ge \text{mode}$.

**Proof (two cases):**

**Case 1 (Low mode):** If $\text{mode}(I(T)) \le \lfloor n/3 \rfloor + 1$:
From Step 2, $k \ge \lfloor n/3 \rfloor + 1 \ge \text{mode}$. Done analytically.

**Case 2 (High mode):** If $\text{mode}(I(T)) > \lfloor n/3 \rfloor + 1$:
We need to show that no maximal IS with all $\text{priv} \le 1$ exists below mode.

**Empirical verification** (`verify_1private_gap.py` + `verify_external_priv.py`):
- Checked all trees through $n = 18$ (204,909 trees).
- No tree with mode $> \lfloor n/3 \rfloor + 1$ admits any maximal IS (of any size) with all $\text{priv} \le 1$.
- Specifically: 0 PNP failures and 0 non-1-Private PNP-like failures across 99,147 maximal IS below mode.

**Note:** The 1-Private IS check (priv = 1 for all vertices) is sufficient because:
- Through n=18, every IS with all priv ≤ 1 below mode is in fact 1-Private (pnp_fail_not_1priv = 0).
- Equivalently, no maximal IS below mode has some vertex with priv = 0 AND all other vertices with priv ≤ 1.
- The 11,265 1-Private IS found all have $k \ge \text{mode}$.

## Structural Insight

1-Private Maximal Independent Sets represent "spread out" coverings:
- Each $u_i \in S$ pairs with a unique private neighbor $v_i$.
- $S \cup P$ induces a matching (plus edges between $v$'s).
- Remaining nodes $R = V \setminus (S \cup P)$ are "shared" (connected to $\ge 2$ nodes in $S$).
This structure forces $S$ to be relatively large (close to $n/2$ unless $|R|$ is very large).
Large independent sets typically fall on the right side of the unimodal curve (or at the peak), satisfying $k \ge \text{mode}$.
Small independent sets (on the left slope, $k < \text{mode}$) correspond to "efficient" domination (stars), where centers cover many private leaves ($\text{priv} \ge 2$).

## Verified Structures in the Gap

The gap cases ($n < 3k, k < \text{mode}$) for $n \le 18$ are dominated by:
- **Caterpillars (Path-Spine):** Spine nodes $u_i$ cover multiple leaves. High priv.
- **Double Stars / D2 / D3:** Centers cover multiple leaves. High priv.
None of these admit a maximal IS with all priv ≤ 1 that is below the mode.
