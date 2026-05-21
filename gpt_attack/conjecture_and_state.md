# Conjecture, definitions, and current state — Erdős #993

## Definitions

- **Tree:** connected acyclic simple graph.
- **Independent set** in a graph G: a subset of vertices no two of which are adjacent.
- **Independence count:** $i_k(G) = |\{S \subseteq V(G) : |S| = k, S \text{ independent}\}|$, with $i_0 = 1$.
- **Independence number:** $\alpha(G) = \max\{k : i_k(G) > 0\}$.
- **Independence polynomial:** $I(G; x) = \sum_{k=0}^{\alpha(G)} i_k(G) x^k$.
- **Mean:** $\mu(G) = \sum_k k \cdot i_k(G) / \sum_k i_k(G) = I'(G; 1) / I(G; 1)$.
- **Mode:** any $k^*$ achieving $\max_k i_k$. May not be unique.
- **Unimodal sequence:** $a_0, a_1, \ldots, a_m$ such that there exists $p$ with $a_0 \leq \cdots \leq a_p \geq \cdots \geq a_m$.
- **Log-concave sequence:** $a_k^2 \geq a_{k-1} a_{k+1}$ for all $1 \leq k \leq m-1$. Log-concavity plus positivity implies unimodality.
- **Leaf-neighbour count** $d_{\mathrm{leaf}}(v)$: the number of leaves adjacent to $v$.
- **Leaf-neighbour bound:** a tree $T$ has $d_{\mathrm{leaf}} \leq 1$ if every vertex has at most one leaf-neighbour.
- **Hub:** a vertex $v$ with $d_{\mathrm{leaf}}(v) \geq 2$.
- **Private neighbour / PNP:** a neighbour of $v$ that is not a neighbour of any other fixed vertex in the context. "PNP framework" refers to the structural lemmas using 1-private maximal independent sets.
- **Edge subdivision** $T_e$: replace edge $e = uv$ with a path $u - w - v$ (one new vertex).
- **Edge contraction** $T/e$: identify $u$ and $v$ into a single vertex with the union of their neighbourhoods.
- **ECMS (Edge Contraction Mode Stability):** a conjectural property stating that mode positions are preserved in a specific sense under contraction of a designated edge class. See manuscript §4 for the precise statement.

## The conjecture

**Conjecture (Alavi, Malde, Schwenk, Erdős 1987).** For every tree T, the independence sequence $i_0(T), i_1(T), \ldots, i_{\alpha(T)}(T)$ is unimodal.

## Related weaker statements

**Conjecture A (Reynolds 2026, manuscript).** For every tree T on n vertices with
$d_{\mathrm{leaf}}(v) \leq 1$ for every vertex v,
$$\mathrm{mode}(I(T)) \leq \lfloor n/3 \rfloor + 1.$$
This is a reduced-regime statement, not a statement about all trees. The global bound
$\mathrm{mode}(I(T)) \leq \lfloor n/3 \rfloor + 1$ is false: for the star $K_{1,d}$,
$I(K_{1,d};x)=(1+x)^d+x$, so the mode is near $d/2$, exceeding
$\lfloor(d+1)/3\rfloor+1$ for large d. Conjecture A is relevant because the PNP
framework and Case B hub bound try to reduce the hard mode-control cases to the
$d_{\mathrm{leaf}}\leq 1$ regime.

**Log-concavity conjecture (historically).** Every tree has a log-concave independence sequence. **Status: FALSE.** Kadrawi–Levit–Yosef–Mizrachi (2023) found exactly two log-concavity failures among the 279 million trees on 26 vertices. Ramos–Sun (2025) found tens of thousands more via ML search.

## What is PROVED in Reynolds (2026)

Use these as black-box lemmas in your attack.

### P1. Subdivision-contraction identity (Theorem 4.3)
For every tree T and edge e of T:
$$I(T_e; x) = I(T; x) + x \cdot I(T/e; x).$$

### P2. Mean bound for d_leaf ≤ 1 trees (Corollary 3.14)
For every tree T on $n \geq 3$ vertices with $d_{\mathrm{leaf}}(v) \leq 1$ for all v (with a stated exception for K_2):
$$\mu(T) < n/3.$$

Proof ingredients: a decimation identity reducing the mean deficit to weighted local deficits, plus Steiner peeling on the heavy-vertex compensation function $F_{\mathrm{gap}}$.

**Note:** the manuscript review has flagged that the theorem statement needs the $n \geq 3$ hypothesis; $K_1$ gives $\mu = 1/2 > 1/3$ and $K_2$ gives $\mu = 2/3 = n/3$ (not strictly less). Handle these small cases separately.

### P3. Hub Exclusion Lemma (Lemma 3.3)
Specific structural property stated in manuscript §3. Used in the reduction to $d_{\mathrm{leaf}} \leq 1$ trees.

### P4. Transfer Lemma (Lemma 3.4)
Specific structural property stated in manuscript §3. Used in the reduction to $d_{\mathrm{leaf}} \leq 1$ trees.

### P5. Spider and caterpillar results (cited from Li et al. 2025, Alavi–Malde–Schwenk–Erdős 1987)
- Every spider has a log-concave (hence unimodal) independence polynomial.
- Every path has a unimodal (even log-concave) independence sequence.
- Every regular caterpillar has a unimodal independence sequence.

### P6. Levit (2006) descending tail
For every tree T, $i_k(T)$ is strictly decreasing for $k \geq \lceil (2\alpha(T) - 1)/3 \rceil$.

### P7. Basit–Galvin (2020) and Heilman (2025) shape results
For a uniformly random labelled tree on $n$ vertices, almost surely the initial ~46.8%–49.5% of the sequence is increasing and the terminal ~38.8% is decreasing.

## What is CONJECTURED (do not assume)

### C1. Conjecture A for the $d_{\mathrm{leaf}}\leq 1$ regime
Stated above. Verified computationally for all $d_{\mathrm{leaf}}\leq 1$ trees through $n = 23$ in the manuscript/certificates. Do not restate this as a global mode bound for all trees; stars are immediate counterexamples to the global statement.

### C2. Conjecture 3.6 (mode-mean)
$\mathrm{mode}(I(T)) \leq \lceil \mu(T) \rceil$ for every tree T. Verified for all trees through $n = 23$ using exact multiprecision arithmetic. If true for $d_{\mathrm{leaf}} \leq 1$ trees, combines with P2 to give Conjecture A without log-concavity.

### C3. ECMS (Edge Contraction Mode Stability)
Mode positions behave stably under designated edge contractions. Verified for 24.7 million edges. If true, subdivision preserves unimodality, so every minimal counterexample is homeomorphically irreducible.

### C4. Combined tail condition
A specific numerical condition on the tail of $I(T; x)$ coefficients (see manuscript §4). Verified in finite ranges. Pairs with ECMS to give the conditional closure of the subdivision route.

### C5. Reduction hypothesis
The informal statement that controlling 1-private maximal independent sets controls the mode position. The manuscript treats this as a reduction; the review board has flagged that a formal theorem is needed. **Producing a formal proof of this reduction would itself be a substantial contribution.**

## What is COMPUTATIONALLY VERIFIED

- Unimodality holds for **all trees with $n \leq 29$ vertices** (8.69 billion trees, exhaustive, on Modal cloud compute with exact integer arithmetic).
- Log-concavity fails for exactly **2 trees at $n = 26$** (Kadrawi–Levit witnesses).
- Log-concavity fails for exactly **19 trees at $n = 28$**, all at index $k = 14$, worst ratio $i_{13} i_{15} / i_{14}^2 = 1.5028$.
- Log-concavity is **clean at $n = 27$** (no failures among ~751 million trees).
- Conjecture 3.6 (mode-mean) verified for all trees through $n = 23$.
- ECMS verified for 24.7 million edges.
- 145,362 structured trees (subdivided stars, caterpillars, spiders, brooms, random perturbations) up to $n = 500$: zero unimodality failures, 378 log-concavity failures (all but 2 in subdivided stars).
- Multi-arm stars identified empirically as the extremal family; near-miss ratio $\mathrm{nm}(s) = 1 - C/s + O(1/s^2)$ with $C \in [4, 8)$ (this constant range is flagged by the review as needing a parity-dependent restatement).

## Existing Lean formalization (`Formal/` directory)

| File | Lines | Sorries | What it holds |
|------|-------|---------|---------------|
| `Formal/Basic.lean` | 610 | 0 | Algebraic infrastructure: polynomial operations, independence-polynomial definition, coefficient accessors. |
| `Formal/P3.lean` | 201 | 1 | Leaf-swap injection. The one sorry is `tree_has_pendant`, which Aristotle previously closed via `SimpleGraph.IsTree.exists_vert_degree_one_of_nontrivial`. |
| `Formal/JleE.lean` | 58 | 0 | Subgraph monotonicity: $J \leq E$. |
| `Formal/Algebra.lean` | 88 | 0 | Star + star degree-2 sum + binomial log-concavity. |
| `Formal/STP2Closure.lean` | 245 | 2 | The open STP2 closure theorem (the main sorry) and a multi-child corollary. Aristotle has closed the base cases (leaf × leaf, degree 0) but not the general proof. |

Mathlib version pinned at v4.28.0. New proofs should target this version.

**Known Lean pitfall:** the `isLC` definition using ℕ subtraction at $k = 0$ makes `lc_conv` (log-concavity preserved under convolution) FALSE because $k = 0$ produces the constraint $f(0) f(1) \leq f(0)^2$ which is strictly stronger than the usual $k \geq 1$ version. Guard with $k \geq 1$ or work over ℤ.

## Known failed approaches (do not repeat)

- Bridge-mode-sum scans (`attack1_bridge_mode_sum_scan_*`)
- Local tuple obstruction (`attack2_local_tuple_obstruction_*`)
- Mean-lift scans (`attack3_mean_lift_scan_*`)
- Approach-A through Approach-D subtree/product/bridge-shift scans (`attack4_*`)
- Weighted compensation probes (`attack_conjA_weighted_compensation.py`)
- ECMS deletion-bridge attacks (`attack_ecms_deletion_bridge.py`)
- Four separate Aristotle runs on STP2 closure (March 2026)
- Brute symbolic manipulation of the subdivision-contraction identity beyond what's in the manuscript

If you find yourself reinventing any of these, stop and pivot.
