# Analytic Proof of the Low-Priv Mode Conjecture

## Theorem (generalized from Gemini's original)
For any tree $T$, if $S$ is a maximal IS with $\text{priv}(u) \le 1$ for all $u \in S$, then $|S| \ge \text{mode}(I(T))$.

**Note:** The original "1-Private Mode Conjecture" assumed priv(u) = 1 for all u (citing Bollobás-Cockayne for priv ≥ 1). In fact, priv(u) = 0 can occur in maximal IS (counterexample: K_{1,4}, S = all leaves). The correct generalization uses priv ≤ 1. See `pnp_n_lt_3k_proof.md` for details.

## Proof Strategy
Divide trees into two classes based on mode:
1. **Low Mode Trees:** $\text{mode}(I(T)) \le \lfloor n/3 \rfloor + 1$.
2. **High Mode Trees:** $\text{mode}(I(T)) > \lfloor n/3 \rfloor + 1$.

---

## Part 1: The Lower Bound

If S is a maximal IS with all priv ≤ 1, then $P = \sum \text{priv}(u) \le k$.
Using the Private Neighbor Bound $P \ge n - 2k + 1$:
$$k \ge P \ge n - 2k + 1 \implies 3k \ge n + 1 \implies k \ge \lceil(n+1)/3\rceil = \lfloor n/3 \rfloor + 1$$

**This bound is tight:** low-mode trees (e.g., long paths with few branches) achieve equality: maximal IS with all priv = 0 and $k = \lceil(n+1)/3\rceil$ exist for n ≥ 8 (`check_low_priv_size.py`).

**Refinement with m_out:** Let $m_{out}$ = edges within $V \setminus S$. Then:
$$n - 1 = E(S, V \setminus S) + m_{out}, \quad E(S, V \setminus S) \ge 2(n-k) - P$$
$$\implies k \ge \frac{n + 1 + m_{out}}{3}$$

Special cases by residue (at the minimum $k = \lceil(n+1)/3\rceil$):
- $n \equiv 2 \pmod{3}$: $m_{out} = 0$ forced, so $V \setminus S$ is independent and $S$ is a color class of $T$, giving $k \ge \lfloor n/2 \rfloor$.
- $n \equiv 1 \pmod{3}$: $m_{out} \le 1$.
- $n \equiv 0 \pmod{3}$: $m_{out} \le 2$.

The $m_{out} = 0$ case is the most powerful: it forces $k \ge \lfloor n/2 \rfloor$, which exceeds $\lfloor n/3 \rfloor + 1$ for $n \ge 8$. However, for $n \not\equiv 2 \pmod{3}$, $m_{out} > 0$ is possible and the color class bound doesn't apply.

---

## Part 2: The Low Mode Case

If $\text{mode}(I(T)) \le \lfloor n/3 \rfloor + 1$: from Part 1, $k \ge \lfloor n/3 \rfloor + 1 \ge \text{mode}$. Done analytically.

---

## Part 3: The High Mode Case

If $\text{mode}(I(T)) > \lfloor n/3 \rfloor + 1$: the Part 1 bound is insufficient (mode exceeds the bound).

### Empirical verification

No tree with mode $> \lfloor n/3 \rfloor + 1$ admits a maximal IS **below mode** with all priv $\le 1$ (verified exhaustively through $n = 18$; 99,147 maximal IS below mode, zero PNP failures of any kind).

**Clarification:** high-mode trees DO have maximal IS with all priv ≤ 1, but these are always well ABOVE mode. Example: $K_{1,10}$ has mode = 5 and a low-priv IS of size 10 (all leaves, all priv = 0). Through n = 16, the minimum $k$ among low-priv IS in high-mode trees exceeds mode by at least 4 (`check_lowpriv_mink_highmode.py`).

### Mode bound

$\text{mode}(I(T)) \le \lfloor (n-1)/2 \rfloor$ for all trees (verified through $n = 20$, 823,065 trees; tight at $K_{1,n-1}$ for odd $n$). The star $K_{1,m}$ achieves $\text{mode} = \lfloor m/2 \rfloor = \lfloor (n-1)/2 \rfloor$.

### Structural analysis (`analyze_high_mode_structure.py`, `analyze_high_mode_leaves.py`)

High-mode trees are star-like:
- 100% have max degree $\ge 4$ (across all n tested)
- Average leaf count is 10-15 for n=14-18
- All have diameter $\le n/2$
- The max-degree vertex ("hub") always has $\ge 2$ leaf-neighbors

Every maximal IS below mode in a high-mode tree **contains the hub** (verified through n=16: without_hub = 0 in all cases). When the hub is in S, its leaf-neighbors are all private (each leaf's only neighbor is the hub, which is in S), giving priv(hub) $\ge$ leaf_deg(hub) $\ge 2$.

### Why the hub is always in the IS below mode

- If hub $\notin S$: all pendant leaf-neighbors of hub must be in $S$ (maximality). Each non-leaf neighbor $w$ of the hub leads to a subtree $T_w$. In each subtree, $S \cap T_w$ is a maximal IS of $T_w$ with $|S \cap T_w| \ge \gamma(T_w)$. So $|S| \ge \ell + \sum_w \gamma(T_w)$ where $\ell$ = leaf-degree of hub.
- If $\ell \ge \text{mode}$: excluding hub forces $|S| \ge \ell \ge \text{mode}$. (True for 57% of high-mode trees at n=17.)
- Even when $\ell < \text{mode}$: the subtree contributions push $|S|$ above mode (empirically, $|S| \ge \text{mode} + 4$).
- Minimum hub_priv = 5 (at n=14) among IS below mode in high-mode trees.

### Why the m_out approach falls short

The $m_{out}$ refinement does NOT close the high-mode gap:
- At $k = \lceil(n+1)/3\rceil$ with all priv $\le 1$: $m_{out} \le 3k - n - 1$ (upper bound from edge counting).
- The constraints are: $k \ge (n + 1 + m_{out})/3$ and $m_{out} \ge 0$. These are circular (substituting back gives $k \ge k$).
- The color class argument ($m_{out} = 0 \implies k \ge \lfloor n/2 \rfloor$) only applies when $m_{out} = 0$ is forced, which requires $n \equiv 2 \pmod{3}$.

### Dead end: "all priv ≤ 1 ⟹ k ≥ ⌊n/2⌋"

This claim is **FALSE** (`check_low_priv_size.py`). Thousands of maximal IS with all priv ≤ 1 have $k < \lfloor n/2 \rfloor$, starting at n=8. However, all such IS have $k \ge \text{mode}$ (consistent with PNP). The low-priv IS with small $k$ only occur in low-mode trees, where the analytic bound suffices.

### Hub exclusion bound FAILS at n=17

The natural proof strategy (hub inclusion lemma) does NOT work:

**Claim (FALSE):** Excluding the hub forces $|S| \ge \text{mode}$ in high-mode trees.

**Counterexample (`analyze_hub_exclusion_bound.py`):** At n=17, a high-mode tree with mode=7 has hub with ell=4 pendant leaves and non-leaf subtrees with $\sum \gamma(T_w) = 2$. The bound $\ell + \sum \gamma(T_w) = 6 < 7 = \text{mode}$. So IS of size 6 (below mode) can exclude the hub.

| n | Bound ok? | min gap |
|---|-----------|---------|
| 11-16 | YES | +1 to +6 |
| 17 | **NO** | -1 |
| 18 | YES | +1 |
| 19 | **NO** | -1 |
| 20 | **NO** | -2 |

Yet PNP still holds at these sizes: the IS that exclude the hub have max_priv $\ge 2$ from **subtree IS vertices**, not from the hub. When the hub is excluded, each non-leaf subtree $T_w$ is a star (since $\gamma(T_w) = 1$), and the center of that star has priv $\ge |T_w| - 1 \ge 2$ (assuming $|T_w| \ge 3$).

### Dead end: Color class bound

**Claim (FALSE):** mode(I(T)) ≤ min(|A|, |B|) where (A, B) is the bipartition of T.

**Reality** (`check_color_class_mode.py`): Fails badly. At n=17: 2,730 trees with mode > min_class (min gap = -7). The star $K_{1,n-1}$ has min_class = 1 but mode = $\lfloor(n-1)/2\rfloor$.

Even though $m_{out} = 0$ forces S to be a color class, knowing S is a color class doesn't help because mode can far exceed the smaller color class size.

### Dead end: Leaf Support Hypothesis

**Claim (FALSE):** Every maximal IS below mode has some $u \in S$ with ≥ 2 leaf-neighbors.

**Counterexample** (Gemini, n=18): Spider with center degree 5, four arms of length 4, one pendant leaf. IS S = {center, four distant feet}. Center has only 1 leaf-neighbor (the pendant), but priv(center) = 5 (all 5 neighbors are private, 4 via internal nodes). See `notes/injection_exploration.md` for details.

Held through n=16 (13,042 cases, 0 failures) but is NOT a universal property.

### Degree Bound Lemma (proved)

**Lemma.** In a tree, for any $u$ in maximal IS $S$ of size $k$: $\text{priv}(u) \ge \deg(u) - (k-1)$.

**Proof.** A neighbor $v$ of $u$ is "shared" if $|N(v) \cap S| \ge 2$, meaning some other $w \in S \setminus \{u\}$ is also adjacent to $v$. By acyclicity, each shared neighbor $v$ maps to a distinct $w \in S \setminus \{u\}$ at distance 2 from $u$ (if two shared neighbors mapped to the same $w$, we'd have a cycle $u$-$v_1$-$w$-$v_2$-$u$). So shared$(u) \le |S \setminus \{u\}| = k-1$, giving priv$(u) = \deg(u) - \text{shared}(u) \ge \deg(u) - (k-1)$. ∎

**Corollary.** If max_deg_in_S $\ge k+1$, then priv $\ge 2$ for the max-degree vertex.

### Disjunction: pigeonhole OR degree bound

The disjunction "n ≥ 3k OR max_deg_in_S ≥ k+1" covers **100% of IS below mode through n=16** (`check_degree_bound_pnp.py`).

**But FAILS at n=17** (`check_degree_bound_n17.py`): 43 uncovered cases where k=6, mode=7, max_deg_in_S = 6 (= k), n = 17 < 3k = 18. In all 43 cases, max_priv ∈ {4, 5, 6} (PNP holds strongly). The high priv comes from a hub vertex with many neighbors that happen to be private (mostly leaf private neighbors, but also internal private neighbors via Gemini's spider mechanism).

### What's actually needed for a proof

Three approaches have been ruled out:
1. **Hub inclusion** (fails n=17)
2. **Color class bound** (mode > min_class for many trees)
3. **Leaf support** (fails n=18)

The Degree Bound Lemma + pigeonhole covers 100% through n=16 but fails at n=17. What remains is to close the gap where deg_in_S = k and n < 3k.

Possible approaches:

1. **Characterize 1-Private sets at the boundary:** If we can show that all trees admitting a maximal IS with all priv ≤ 1 and k = ⌈(n+1)/3⌉ have mode ≤ k, the proof is complete. These tight cases appear to be exclusively double-stars.

2. **Refined shared-neighbor bound:** In the 43 gap cases, the actual shared count is much less than k-1. The degree bound is conservative because it assumes each shared neighbor maps to a distinct S-vertex, but in practice the distance-2 S-vertices are concentrated. A tighter analysis of tree structure may improve the bound.

3. **Global edge structure:** When k is in the gap [⌈(n+1)/3⌉, mode-1], the tree structure prevents all privs from being ≤ 1 by forcing enough edge concentration on the highest-degree S-vertex.

---

## Final Conclusion

1. $\text{mode} \le \lfloor n/3 \rfloor + 1$: $k \ge \lfloor n/3 \rfloor + 1 \ge \text{mode}$. (Analytic)
2. $\text{mode} > \lfloor n/3 \rfloor + 1$: no IS with all priv $\le 1$ exists below mode. (Empirical, n $\le 18$)

PNP is fully verified through n=18 (99,147 maximal IS, 0 failures).
mode $\le \lfloor(n-1)/2\rfloor$ verified through n=20 (823,065 trees, 0 exceptions).

## Scripts
- `check_mode_vs_half.py`: mode ≤ ⌊n/2⌋ verification (n ≤ 20)
- `check_low_priv_size.py`: "all priv ≤ 1 ⟹ k ≥ ⌊n/2⌋" test (FAILS)
- `check_highmode_lowpriv.py`: low-priv IS in high-mode trees (exist but above mode)
- `check_lowpriv_mink_highmode.py`: min k gap in high-mode trees (≥ mode + 4)
- `analyze_high_mode_structure.py`: structural characterization of high-mode trees
- `analyze_high_mode_leaves.py`: hub leaf-degree and inclusion analysis
- `analyze_hub_exclusion_bound.py`: hub exclusion bound (FAILS at n=17)
- `analyze_recursive_pnp.py`: source of priv ≥ 2 (hub provides, degree bound)
- `check_degree_bound_pnp.py`: pigeonhole + degree bound disjunction (100% through n=16)
- `check_degree_bound_n17.py`: degree bound at n=17-18 (43 uncovered at n=17)
- `check_color_class_mode.py`: color class bound (FAILS: mode > min_class)
- `analyze_gap_cases.py`: structural analysis of 43 gap cases at n=17
- `verify_external_priv.py`: priv=0 exists below mode but PNP holds (n ≤ 18)
