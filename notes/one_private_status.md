# Status of the 1-Private Mode Conjecture

## Current Status: Conditional Two-Conjecture PNP Framework (Corrected 2026-07-10)

**Conjecture (PNP):** For any tree $T$ and any maximal IS $S$ with $|S| < \text{mode}(I(T))$, some $u \in S$ has $\text{priv}(u) \ge 2$.

Equivalently: if every $u \in S$ has $\text{priv}(u) \le 1$, then $|S| \ge \text{mode}(I(T))$.

**Current dependency correction.** Hub Exclusion and Transfer are proved, but
PNP does not reduce to Conjecture A alone. The current manuscript correctly
requires both:

1. Conjecture A for the $d_{\text{leaf}}\le1$ lane; and
2. a separate Case-B hub bound comparing $\operatorname{mode} I(T)$ with the
   explicit hub-reduction size bound.

Transfer preserves the 1-Private property on the residual forest; it does not
by itself compare the mode of $I(T)$ with modes or size bounds of residual
components. Even a proof of PNP would constrain below-mode maximal independent
sets rather than directly prove unimodality. The dependency map in
`paper/main_v2.tex` is authoritative.

---

## Proof Framework (Revised 2026-02-14)

### Part 1: Size Bound (PROVED)

**Lemma.** Any maximal IS with all $\text{priv}(u) \le 1$ has $k \ge \lceil (n+1)/3 \rceil$.

*Proof:* $P = \sum \text{priv}(u) \le k$. Private Neighbor Bound: $P \ge n - 2k + 1$.
So $k \ge n - 2k + 1$, giving $k \ge \lceil (n+1)/3 \rceil$. $\square$

### Part 2: Low-Mode Case (PROVED)

**Lemma.** If $\text{mode}(I(T)) \le \lfloor n/3 \rfloor + 1$, then $k \ge \text{mode}$.

*Proof:* From Part 1, $k \ge \lceil (n+1)/3 \rceil = \lfloor n/3 \rfloor + 1 \ge \text{mode}$. $\square$

### Leaf Injection Lemma (PROVED)

**Lemma.** For any tree $T$ and any maximal IS $S$ with all $\text{priv}(u) \le 1$,
$|S| \ge \ell(T)$ (number of leaves).

*Proof:* Define $f: \text{Leaves}(T) \to S$ by:
- If leaf $v \in S$: $f(v) = v$.
- If leaf $v \notin S$: $f(v) = u$ where $u = N(v)$ (the support vertex; $u \in S$ since $S$ dominates and $v$ has no other neighbor).

$f$ is injective:
1. Two leaves both in $S$ map to distinct vertices (trivially).
2. A leaf in $S$ (degree 1) maps to itself; a leaf not in $S$ maps to its support (degree $\ge 2$ for $n \ge 3$). These are distinct.
3. Two leaves $v_1, v_2 \notin S$ with $f(v_1) = f(v_2) = u$ would give $d_{\text{leaf}}(u) \ge 2$. But $u \in S$ and all leaf-children of $u$ are outside $S$ (adjacent to $u$), so $\text{priv}(u) \ge d_{\text{leaf}}(u) \ge 2$, contradicting $\text{priv} \le 1$.

So $|S| \ge |\text{Leaves}(T)|$. $\square$

### Hub Exclusion Lemma (PROVED)

**Lemma.** For any 1-Private maximal IS $S$, if $v$ has $d_{\text{leaf}}(v) \ge 2$ leaf-children, then $v \notin S$. Consequently, all $d_{\text{leaf}}(v)$ leaf-children of $v$ are in $S$.

*Proof:* If $v \in S$, all leaf-children are adjacent to $v$ hence outside $S$. Each leaf-child $w$ has $N(w) \cap S = \{v\}$, so $w$ is private to $v$. Thus $\text{priv}(v) \ge d_{\text{leaf}}(v) \ge 2$, contradicting 1-Private.

Since $v \notin S$, each leaf-child $w$ has $N(w) = \{v\}$ and $v \notin S$. For domination, $w$ must be in $S$. $\square$

### Part 3: High-Mode Case (Two Sub-Cases)

For trees with $\text{mode} > \lfloor n/3 \rfloor + 1$:

#### Case A: All $d_{\text{leaf}} \le 1$ (CONJECTURE)

**Conjecture.** If every vertex of $T$ has $d_{\text{leaf}}(v) \le 1$, then $\text{mode}(I(T)) \le \lfloor n/3 \rfloor + 1$.

*Status:* Verified exhaustively through $n = 22$ (227,678 such trees, **zero** with mode above threshold). Also verified for large constructed examples (caterpillars, lobsters, etc.) up to $n = 1000$ with mode/threshold ratios always $< 0.95$.

**Mean-bound status:** The current manuscript proves
$\mu(I(T))<n/3$ for every $d_{\text{leaf}}\le1$ tree by Steiner peeling.
The remaining implication is the independent mode--mean conjecture
$\operatorname{mode} I(T)\le\lceil\mu(T)\rceil$, or the narrower
tie-fugacity condition $\mu(\lambda_m)\ge m-1$. Neither follows from positive
log-concavity alone. Moreover, $d_{\text{leaf}}\le1$ log-concavity is false:
the July 2026 Ramos--Sun stress corpus contains many such non-log-concave
trees.

**Extremal family:** Spiders $S(2^k)$ (hub + $k$ arms of length 2). Algebraically:
$$I(S(2^k); x) = (1+2x)^k + x(1+x)^k, \quad n = 2k+1$$
$$\mu = \frac{2k \cdot 3^{k-1} + 2^k + k \cdot 2^{k-1}}{3^k + 2^k}, \quad \frac{n}{3} - \mu = \frac{3^k + 2^{k-1}(k-4)}{3(3^k + 2^k)} \to \frac{1}{3}$$
Positive for $k \ge 5$ ($n \ge 11$); direct check for $k \le 4$. Second-tightest: $S(2^k, 1^1)$, gap $\to 1/6$.

**Per-vertex bound $P(v \in S) \le 1/3$ FAILS** (leaves reach $P \approx 0.49$). The bound $\mu < n/3$ is a global property, not a per-vertex one. Inductive arguments via vertex deletion are problematic because deleting a vertex can break the $d_{\text{leaf}} \le 1$ property.

If the mode--mean or tie-fugacity bridge holds, the proved mean bound and
Part 2 handle this lane.

#### Case B: Some $d_{\text{leaf}} \ge 2$ (SEPARATE HUB BOUND OPEN)

**Setup:** Let $v$ have $d_{\text{leaf}}(v) = d \ge 2$. By Hub Exclusion, $v \notin S$ and all $d$ leaf-children of $v$ are in $S$.

**Transfer Lemma (PROVED).** Let $T' = T - \{v, w_1, \ldots, w_d\}$ (residual forest). Then $S' = S \cap V(T')$ is a 1-Private maximal IS in $T'$.

*Proof:*
1. **Maximal:** For any $w \in V(T') \setminus S$, adding $w$ to $S$ creates an adjacency in $T$ with some $u \in S$. If $u$ is a leaf-child of $v$, then $w$ is adjacent to that leaf, so $w = v$ (leaves have degree 1). But $v \notin V(T')$. So $u \in S' \subseteq V(T')$, and $w$ conflicts within $T'$. $\square$
2. **1-Private transfers:** For $u \in S'$, its private neighbors in $T$ are $\text{priv}_T(u) = \{w \in N_T(u) : N_T(w) \cap S = \{u\}\}$. Since leaf-children of $v$ have $N(w_i) = \{v\}$ and $v \notin S$, leaf-children of $v$ are NOT private neighbors of any $u \in S'$. Since $v \notin S$, $v$ is not private to any $u$ either (even if $u$ is adjacent to $v$, because $v$ has multiple $S$-neighbors: the $d \ge 2$ leaf-children). So $\text{priv}_{T'}(u) = \text{priv}_T(u) \cap V(T')$, and $|\text{priv}_{T'}(u)| \le |\text{priv}_T(u)| \le 1$. $\square$

**Bound.** By Part 1 applied to each component of $T'$ (each component $C_j$ has $|S' \cap C_j| \ge \lceil(|C_j|+1)/3\rceil$):
$$k = d + |S'| \ge d + \sum_j \lceil(|C_j|+1)/3\rceil \ge d + \lceil(n - d - h + c)/3\rceil$$
where $h$ = number of hubs, $c$ = components of $T'$, and $c \ge h \ge 1$.

**Verification:** $\text{mode}(I(T)) \le F + \lceil(n - F - h + c)/3\rceil$ holds for all 8,710,881 trees with $d_{\text{leaf}} \ge 2$ through $n = 22$ (zero violations, minimum surplus $= 1$).

**Open Case-B comparison.** The displayed computational inequality is the
separate Case-B hub-bound conjecture. Recursively transferring $S$ to smaller
components proves lower bounds on $|S|$, but it does not supply the missing
upper bound on $\operatorname{mode} I(T)$. Conjecture A therefore does not
close this case by itself.

---

## RETRACTED: Leaf-Mode Inequality

~~**Conjecture.** For any tree $T$ with $\text{mode} > \lfloor n/3 \rfloor + 1$, $\ell(T) \ge \text{mode}$.~~

**FALSE.** Counterexample found (Gemini, confirmed independently):

- $S(41, 1^{20})$: $n = 62$, $\ell = 21$, $\text{mode} = 22$, threshold $= 21$.
- Mode exceeds both $\ell$ and threshold.
- 812 violations found among spiders through $n = 150$, worst gap $= -7$.
- Pattern: brooms/spiders with one long arm + many unit arms.

**But PNP still holds for all counterexample trees:** the Hub Exclusion Lemma forces $k \ge d + \gamma(T') \gg \ell$ for exactly these trees.

---

## Combining the Proved Parts

| Case | Condition | Why $k \ge \text{mode}$ | Status |
|------|-----------|------------------------|--------|
| Low-mode | $\text{mode} \le \lfloor n/3 \rfloor + 1$ | Part 1: $k \ge \lfloor n/3 \rfloor + 1$ | PROVED |
| High-mode, all $d_{\text{leaf}} \le 1$ | mode $> \lfloor n/3 \rfloor + 1$ | Conjecture A: can't happen | CONJECTURE |
| High-mode, some $d_{\text{leaf}} \ge 2$ | mode $> \lfloor n/3 \rfloor + 1$ | Transfer + separate Case-B hub comparison | CONJECTURE |

**PNP is conditional on both Conjecture A and the Case-B hub bound.** The
proved Transfer Lemma is a necessary combinatorial reduction, not the missing
mode comparison.

---

## Computational Evidence Summary

| Check | Range | Result |
|-------|-------|--------|
| PNP (priv $\ge 2$ below mode) | $n \le 18$, 99,147 IS | 0 failures |
| 1-Private below mode | $n \le 18$, 2,366,860 IS | 0 exist |
| $k \ge \ell$ for 1-Private | $n \le 18$, 2,366,860 IS | 0 violations |
| $d_{\text{leaf}} \le 1 \Rightarrow$ low-mode | $n \le 22$, 227,678 trees | 0 violations |
| $d_{\text{leaf}} \le 1 \Rightarrow$ low-mode (large) | families to $n = 1000$ | mode/thr $\le 0.944$ |
| $\mu < n/3$ for $d_{\text{leaf}} \le 1$ | $n \le 20$, 43,029 trees | max $\mu/(n/3) = 0.973$ |
| $\mu < n/3$ for $S(2^k)$ | algebraic, all $k$ | gap $= [3^k + 2^{k-1}(k{-}4)]/[3(3^k{+}2^k)] > 0$ |
| $F + \lceil(n{-}F{-}h{+}c)/3\rceil \ge \text{mode}$ | $n \le 22$, 8,710,881 trees | 0 violations, min surplus $= 1$ |
| Spider extremality ($\mu$) | $n \le 20$, 43,029 trees | best $d_{\text{leaf}} \le 1$ always spider |
| Core contraction $\Rightarrow$ $\mu/n$ increase | $n \le 16$ | **FAILS** (18% decrease rate) |
| Leaf-Mode (spiders, large $n$) | $n \le 150$ | **812 violations** (RETRACTED) |

---

## What Remains

**Two open steps for PNP: Conjecture A and the Case-B hub bound.**

**Conjecture A.** The mean inequality $\mu(I(T))<n/3$ is proved. What remains
is to prove mode--mean localization, or the narrower tie-fugacity condition,
for $d_{\text{leaf}}\le1$ trees. Real-rootedness cannot be invoked: general
tree independence polynomials are not real-rooted, and the lane is not even
log-concave in general.

**Case B.** Prove
$\operatorname{mode} I(T)\le F+\lceil(n-F-h+c)/3\rceil$ for trees with
multi-leaf hubs, or find a replacement comparison strong enough for the PNP
corollary.

Finally, PNP is not itself a proof of independence-sequence unimodality; it is
one structural ingredient in the manuscript's broader program.

**Failed approaches:**
- Per-vertex bound $P(v \in S) \le 1/3$: FAILS (leaves reach $P \approx 0.49$)
- Induction via vertex deletion: FAILS (can break $d_{\text{leaf}} \le 1$)
- Matching-based edge decomposition: gives $\mu \le 2n/3$ (too weak by factor 2)
- Core edge contraction increasing $\mu/n$: FAILS (18% decrease rate from $n = 9$)
- Spider extremality via contraction: dead end

**Current directions:**
- (a) Prove the tie-fugacity condition
  $\mu(\lambda_m)\ge m-1$ for $d_{\text{leaf}}\le1$ trees. This is the
  narrowest live bridge from the proved mean bound to Conjecture A.
- (b) Prove the broader mode--mean inequality
  $\operatorname{mode} I(T)\le\lceil\mu(T)\rceil$, without importing
  real-rootedness or log-concavity.
- (c) Prove the separate Case-B hub comparison. The signed
  Poisson-binomial reserve program is a bounded analytic subroute for
  hub-bouquet families. Its quarter-scale effective-drop and raw-reserve
  theorem is now proved for every finite Poisson-binomial law, closing the
  signed bridge and the product term `A=(1+x)^sQ`. The hub-included
  perturbation `xR` remains open for growing arms, and the result has not been
  lifted to arbitrary multi-hub trees or the Case-B mode bound.

---

## Previous Approaches (Superseded)

- **Leaf-Mode Inequality**: FALSE at $n = 62$. Brooms with long handles are counterexamples.
- **Core edge contraction**: Contracting an edge between two non-leaf vertices does NOT always increase $\mu/n$. Failures at $n \ge 9$ (18% rate), always involving degree-2 vertex adjacent to degree-${\ge}3$. Raw $\mu$ also decreases (not just ratio). Dead end for spider extremality.
- **Expansion Property**: Proved analytically but addresses a scenario (pendant IS below mode) that never arises.
- **Hub inclusion / degree bound / color class**: Dead ends at various $n$.
- **Gemini polynomial decomposition**: Only valid for pendant-pair trees (special case).
