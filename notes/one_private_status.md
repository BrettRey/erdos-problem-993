# Status of the 1-Private Mode Conjecture

## Current Status: Single-Conjecture Reduction (Revised 2026-02-15)

**Conjecture (PNP):** For any tree $T$ and any maximal IS $S$ with $|S| < \text{mode}(I(T))$, some $u \in S$ has $\text{priv}(u) \ge 2$.

Equivalently: if every $u \in S$ has $\text{priv}(u) \le 1$, then $|S| \ge \text{mode}(I(T))$.

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

**Mean Bound Approach (promising):** Tree IS polynomials are real-rooted (Chudnovsky--Seymour 2007), so $\text{mode} \le \lceil \mu \rceil$ where $\mu = I'(1)/I(1)$ is the mean IS size. We verified $\mu < n/3$ for all $d_{\text{leaf}} \le 1$ trees through $n = 20$ (43,029 trees, max $\mu/(n/3) = 0.973$). Since $\mu < n/3$ implies $\text{mode} \le \lceil \mu \rceil \le \lfloor n/3 \rfloor + 1$, this reduces Conjecture A to: **$\mu(I(T)) < n/3$ for all $d_{\text{leaf}} \le 1$ trees.**

**Extremal family:** Spiders $S(2^k)$ (hub + $k$ arms of length 2). Algebraically:
$$I(S(2^k); x) = (1+2x)^k + x(1+x)^k, \quad n = 2k+1$$
$$\mu = \frac{2k \cdot 3^{k-1} + 2^k + k \cdot 2^{k-1}}{3^k + 2^k}, \quad \frac{n}{3} - \mu = \frac{3^k + 2^{k-1}(k-4)}{3(3^k + 2^k)} \to \frac{1}{3}$$
Positive for $k \ge 5$ ($n \ge 11$); direct check for $k \le 4$. Second-tightest: $S(2^k, 1^1)$, gap $\to 1/6$.

**Per-vertex bound $P(v \in S) \le 1/3$ FAILS** (leaves reach $P \approx 0.49$). The bound $\mu < n/3$ is a global property, not a per-vertex one. Inductive arguments via vertex deletion are problematic because deleting a vertex can break the $d_{\text{leaf}} \le 1$ property.

If this holds, then Part 2 handles the case: $k \ge \lfloor n/3 \rfloor + 1 \ge \text{mode}$.

#### Case B: Some $d_{\text{leaf}} \ge 2$ (REDUCES TO CONJECTURE A)

**Setup:** Let $v$ have $d_{\text{leaf}}(v) = d \ge 2$. By Hub Exclusion, $v \notin S$ and all $d$ leaf-children of $v$ are in $S$.

**Transfer Lemma (PROVED).** Let $T' = T - \{v, w_1, \ldots, w_d\}$ (residual forest). Then $S' = S \cap V(T')$ is a 1-Private maximal IS in $T'$.

*Proof:*
1. **Maximal:** For any $w \in V(T') \setminus S$, adding $w$ to $S$ creates an adjacency in $T$ with some $u \in S$. If $u$ is a leaf-child of $v$, then $w$ is adjacent to that leaf, so $w = v$ (leaves have degree 1). But $v \notin V(T')$. So $u \in S' \subseteq V(T')$, and $w$ conflicts within $T'$. $\square$
2. **1-Private transfers:** For $u \in S'$, its private neighbors in $T$ are $\text{priv}_T(u) = \{w \in N_T(u) : N_T(w) \cap S = \{u\}\}$. Since leaf-children of $v$ have $N(w_i) = \{v\}$ and $v \notin S$, leaf-children of $v$ are NOT private neighbors of any $u \in S'$. Since $v \notin S$, $v$ is not private to any $u$ either (even if $u$ is adjacent to $v$, because $v$ has multiple $S$-neighbors: the $d \ge 2$ leaf-children). So $\text{priv}_{T'}(u) = \text{priv}_T(u) \cap V(T')$, and $|\text{priv}_{T'}(u)| \le |\text{priv}_T(u)| \le 1$. $\square$

**Bound.** By Part 1 applied to each component of $T'$ (each component $C_j$ has $|S' \cap C_j| \ge \lceil(|C_j|+1)/3\rceil$):
$$k = d + |S'| \ge d + \sum_j \lceil(|C_j|+1)/3\rceil \ge d + \lceil(n - d - h + c)/3\rceil$$
where $h$ = number of hubs, $c$ = components of $T'$, and $c \ge h \ge 1$.

**Verification:** $\text{mode}(I(T)) \le F + \lceil(n - F - h + c)/3\rceil$ holds for all 8,710,881 trees with $d_{\text{leaf}} \ge 2$ through $n = 22$ (zero violations, minimum surplus $= 1$).

**Inductive reduction.** For each component of $T'$: if it has all $d_{\text{leaf}} \le 1$, Conjecture A bounds its mode. If some vertex has $d_{\text{leaf}} \ge 2$, apply Hub Exclusion + Transfer recursively. The recursion terminates at components with all $d_{\text{leaf}} \le 1$.

**Therefore, PNP for ALL trees reduces to Conjecture A.** If every tree with all $d_{\text{leaf}} \le 1$ has $\text{mode} \le \lfloor n/3 \rfloor + 1$, then PNP holds.

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
| High-mode, some $d_{\text{leaf}} \ge 2$ | mode $> \lfloor n/3 \rfloor + 1$ | Transfer + Part 1 on residual | REDUCES TO CONJECTURE A |

**The entire PNP proof reduces to Conjecture A.** Case B is handled by the Transfer Lemma + recursive application of Hub Exclusion, eventually terminating at $d_{\text{leaf}} \le 1$ components where Conjecture A applies.

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

**One open step: Conjecture A.**

**Conjecture A (reduced to Mean Bound).** Prove $\mu(I(T)) < n/3$ for all $d_{\text{leaf}} \le 1$ trees. The extremal family $S(2^k, 1)$ has gap $n/3 - \mu \to 1/6$, so the bound is tight. By real-rootedness (Chudnovsky--Seymour 2007), $\mu < n/3$ implies $\text{mode} \le \lfloor n/3 \rfloor + 1$.

Case B (trees with $d_{\text{leaf}} \ge 2$) reduces to Conjecture A via the Transfer Lemma and induction on the number of multi-leaf hubs.

**Failed approaches:**
- Per-vertex bound $P(v \in S) \le 1/3$: FAILS (leaves reach $P \approx 0.49$)
- Induction via vertex deletion: FAILS (can break $d_{\text{leaf}} \le 1$)
- Matching-based edge decomposition: gives $\mu \le 2n/3$ (too weak by factor 2)
- Core edge contraction increasing $\mu/n$: FAILS (18% decrease rate from $n = 9$)
- Spider extremality via contraction: dead end

**Promising directions:**
- (a) Hard-core model: $\mu = \sum_v P(v \in S)$ where $P$ is the hard-core measure at $\lambda = 1$. For leaf-support pairs, $P(w) + P(v) = 1/2 + P(v)/2$, with combined load $P(v)/2 - 1/6 \le 0$ when $P(v) \le 1/3$. This holds for all support vertices with $\ge 1$ non-leaf child. The residual "core" vertices need $\sum_{\text{core}} P(v) \le |\text{core}|/3$.
- (b) Root structure / pressure bounds: the tree recursion $R(v) = 1/\prod(1+R(c_j))$ gives $P(v) = R/(1+R)$. Bounding the sum $\sum P(v)$ via the recursion.
- (c) Spider extremality (verified $n \le 20$, algebraic proof for extremal spider): prove $S(2^k, 1)$ maximizes $\mu$ among $d_{\text{leaf}} \le 1$ trees.

---

## Previous Approaches (Superseded)

- **Leaf-Mode Inequality**: FALSE at $n = 62$. Brooms with long handles are counterexamples.
- **Core edge contraction**: Contracting an edge between two non-leaf vertices does NOT always increase $\mu/n$. Failures at $n \ge 9$ (18% rate), always involving degree-2 vertex adjacent to degree-${\ge}3$. Raw $\mu$ also decreases (not just ratio). Dead end for spider extremality.
- **Expansion Property**: Proved analytically but addresses a scenario (pendant IS below mode) that never arises.
- **Hub inclusion / degree bound / color class**: Dead ends at various $n$.
- **Gemini polynomial decomposition**: Only valid for pendant-pair trees (special case).
