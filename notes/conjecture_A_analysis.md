# Conjecture A Analysis: d_leaf ≤ 1 Trees Are Low-Mode

## Statement

**Conjecture A.** If every vertex of tree $T$ has $d_{\text{leaf}}(v) \le 1$, then $\text{mode}(I(T)) \le \lfloor n/3 \rfloor + 1$.

## Current Reduction (Corrected 2026-07-10)

### Step 1: Mean bound (PROVED)

The current manuscript proves

$$
\mu(I(T))<n/3
$$

for every $d_{\text{leaf}}\le1$ tree by a decimation identity and Steiner
peeling. Spider extremality is no longer needed to establish this inequality.

### Step 2: Mode--mean or tie-fugacity bridge (OPEN)

Conjecture A would follow from either

$$
\operatorname{mode} I(T)\le\lceil\mu(T)\rceil
$$

or the narrower tie-fugacity condition

$$
\mu(\lambda_m)\ge m-1.
$$

Neither statement follows from positive log-concavity alone. For example,
$(1,2,3,4,5)$ is strictly log-concave but has mode $4$ and mean $8/3$.
Moreover, the July 2026 Ramos--Sun stress corpus contains many
$d_{\text{leaf}}\le1$ trees whose independence sequences are not
log-concave. Thus no log-concavity shortcut remains.

### Historical Step: Spider Extremality (CONJECTURE, verified n ≤ 20)

**Conjecture A1.** Among all $d_{\text{leaf}} \le 1$ trees on $n$ vertices, spiders maximize $\mu$.

Verified exhaustively through $n = 20$ (43,029 trees at $n = 20$). The best $d_{\text{leaf}} \le 1$ tree is always a spider, at every $n$ from 5 to 20.

### Historical Step: Optimal Spider Configuration (verified n ≤ 20)

Among $d_{\text{leaf}} \le 1$ spiders on $n$ vertices, $\mu$ is maximized by:
- $S(2^k, 1)$ for even $n = 2k + 2$
- $S(3, 2^{k-1}, 1)$ for odd $n = 2k + 3$

Both have one arm of length 1 (giving $d_{\text{leaf}}(\text{hub}) = 1$) and remaining arms of length 2 (plus one of length 3 for odd $n$).

The arm-1 "trick": replacing two arms of length 2 with one of length 3 and one of length 1 increases $\mu$, because a pendant edge has disproportionately high per-vertex IS contribution. The $d_{\text{leaf}} \le 1$ constraint limits this to at most one arm of length 1.

### Calibration: Mean Bound for Optimal Spider (PROVED)

**Lemma.** $\mu(S(2^k, 1)) < n/3$ for all $k \ge 1$.

*Proof.* Let $n = 2k + 2$. The IS polynomial is:
$$I(x) = (1 + x)(1 + 2x)^k + x(1 + x)^k$$

Computing:
$$I(1) = 2 \cdot 3^k + 2^k, \quad I'(1) = 3^k + 4k \cdot 3^{k-1} + 2^k + k \cdot 2^{k-1}$$

The gap:
$$\frac{n}{3} - \mu = \frac{3^k + 2^{k-1}(k - 2)}{3(2 \cdot 3^k + 2^k)}$$

For $k \ge 3$: numerator $= 3^k + 2^{k-1}(k - 2) > 0$. Direct check for $k = 1, 2$. $\square$

Asymptotically: gap $\to 1/6$ as $k \to \infty$.

**For odd $n$:** $S(3, 2^{k-1}, 1)$ has even higher gap (the arm-3 contribution adds to $\mu$ less efficiently than arm-2). The analysis is similar.

### Summary of Proof Status

| Step | Statement | Status |
|------|-----------|--------|
| 1 | $\mu(T)<n/3$ for $d_{\text{leaf}}\le1$ | PROVED (Steiner peeling) |
| 2 | $\text{mode} \le \lceil \mu \rceil$ or tie-fugacity bridge | OPEN |
| 3 | Spiders maximize $\mu$ among $d_{\text{leaf}} \le 1$ | HISTORICAL CONJECTURE (not needed for Step 1) |
| 4 | $\mu(S(2^k, 1)) < n/3$ | PROVED calibration |

**Open:** Step 2. The proved mean bound plus either mode bridge would prove
Conjecture A. Conjecture A is only the $d_{\text{leaf}}\le1$ part of the
conditional PNP framework; a separate Case-B hub bound remains open.

---

## Historical Failed Approaches for Direct Proof of μ < n/3

### Per-vertex bound P(v ∈ S) ≤ 1/3

**FAILS.** Leaves can have $P(v \in S) \approx 0.49$ in $d_{\text{leaf}} \le 1$ trees. The bound $\mu < n/3$ is a global property (the sum is bounded) even though individual terms exceed $1/3$.

### Induction via vertex deletion

**FAILS.** Deleting a vertex from a $d_{\text{leaf}} \le 1$ tree can create $d_{\text{leaf}} \ge 2$ vertices. Example: in path $P_4 = a{-}b{-}c{-}d$, deleting leaf $d$ makes $c$ a leaf, giving $d_{\text{leaf}}(b) = 2$ (both $a$ and $c$ are now leaf-children of $b$).

### Leaf-support pair bound ≤ 2/3

For adjacent pair $(w, v)$ where $w$ is a leaf: $P(w) + P(v) \le 2/3$ would give $\mu \le n/3$ by pairing. But this bound FAILS: for paths, $P(w) + P(v) \approx 0.62$, which is close to $2/3$ but not always below it in non-path trees.

### Matching-based edge decomposition

Gives only $\mu \le 2n/3$ (each edge contributes $\le 1$). Too weak by factor 2.

---

## Why Spiders Are Extremal (Intuition)

A spider concentrates all branching at a single hub. This creates a product structure $I = \prod f_i + x \prod g_i$ that maximizes large IS (pushing the mean up). A non-spider tree distributes branching across multiple vertices, which constrains IS more and reduces $\mu$.

Formally: "collapsing" the core of a $d_{\text{leaf}} \le 1$ tree (contracting internal edges until one hub remains) should increase $\mu$. This would follow from showing that edge contraction in the core increases $\mu$, but a clean proof is elusive.

---

## Promising Proof Directions for Spider Extremality

1. **Compression argument:** Show that contracting any edge $\{u,v\}$ in the core (where both $u,v$ have $d_{\text{leaf}} = 0$ and degree $\ge 2$) increases $\mu$. This would give: any $d_{\text{leaf}} \le 1$ tree can be "compressed" to a spider while monotonically increasing $\mu$.

2. **Galvin-type extremality:** Galvin (2011) showed that among $d$-regular graphs, the $d$-regular tree maximizes $|I(G;1)|$. An analogous result for $\mu$ might establish spider extremality.

3. **Convexity of the arm contribution:** Show that the arm's contribution to $\mu$ is concave in the arm length, so concentrating all arm length at $a = 2$ (or splitting into $2 + 1$) is optimal.

4. **Hard-core model entropy:** In the hard-core model on trees, there are known results about occupation probabilities. The mean $\mu$ is the expected IS size at fugacity $\lambda = 1$. Bounds on the pressure/free energy might give $\mu < n/3$.
