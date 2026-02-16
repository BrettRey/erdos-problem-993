# SUPERSEDED — PROOF INCOMPLETE

**This file claims an analytic proof but contains a critical gap in Step 3.**
The claim that "the sum of log-concave polynomials sharing a common
bound on their modes will generally satisfy the same bound" is not
proved and is not true in general (sums of unimodal/LC sequences can
fail to be unimodal).

**See `notes/one_private_status.md` for the definitive, honest status
of the PNP framework.**

---

# Attempted Proof of the 1-Private Mode Conjecture (INCOMPLETE)

## Conjecture (not theorem)
For any tree $T$, if $S$ is a **1-Private Maximal Independent Set** (where every vertex $u \in S$ has at most 1 private neighbor), then $|S| \ge \text{mode}(I(T))$.

---

## 1. The Pendant-Pair Structure (Analytic Derivation)

From edge counting ($n-1 \ge |P| + 2|R| + m_{out}$), we derived $k \ge \lceil (n+1)/3 \rceil$.
In the "tight" case ($k \approx n/3$), we deduced:
*   $R$ is the spine ($V \setminus (S \cup P)$).
*   $S \cup P$ forms $k$ disjoint edges $(u_i, p_i)$ where $p_i$ are leaves attached to $S$.
*   Connections exist between $S$ and $R$, but $P$ is isolated from $R$.

---

## 2. The Complement Expansion Theorem

**Lemma:** Let $X \subseteq S$. Let $Y^c = R \setminus N(X)$ be the set of spine nodes *not* connected to $X$. Let $X^c = S \setminus X$.
Then **$|Y^c| < |X^c|$**.

*Proof:* Via edge counting in the forest induced by $X^c \cup Y^c$, showing $2|Y^c| \le |E'| < |X^c| + |Y^c|$.

---

## 3. The Mode Bound (Refined)

We decompose the independence polynomial $I(T; x)$ by conditioning on the intersection with $S$.
For a chosen subset $X \subseteq S$:
1.  **Vertices in $X$:** Contribute factor $x^{|X|}$.
2.  **Pendants of $X$:** Blocked (neighbors of $S$). Contribution 1.
3.  **Pendants of $X^c$:** Free. Each $u \in X^c$ is *not* chosen, so its pendant $p_u$ is free to be chosen or not. Contribution $(1+x)^{|X^c|} = (1+x)^{k-|X|}$.
4.  **Spine:** Only $Y^c = R \setminus N(X)$ is free. Contribution $I(T[Y^c]; x)$.

$$ I(T; x) = \sum_{X \subseteq S} x^{|X|} (1+x)^{k-|X|} I(T[Y^c]; x) $$

**Peak Analysis:**
Let $P_X(x) = x^{|X|} (1+x)^{k-|X|} I(T[Y^c]; x)$.
Since $T[Y^c]$ is a forest, $I(T[Y^c])$ is log-concave. The term $(1+x)^{k-|X|}$ is log-concave. The product is log-concave (unimodal).
The mode of the product is approximately the sum of the modes of the factors (shifted by $x^{|X|}$).

$$ \text{Peak}_X \approx |X| + \frac{k-|X|}{2} + \text{mode}(I(T[Y^c])) $$

**Bound on Spine Mode:**
For any forest $F$, $\text{mode}(I(F)) \le \lceil |V(F)|/2 \rceil$.
So $\text{mode}(I(T[Y^c])) \le \lceil |Y^c|/2 \rceil$.

$$ \text{Peak}_X \le |X| + \frac{k-|X|}{2} + \frac{|Y^c|}{2} $$
$$ \text{Peak}_X = \frac{k}{2} + \frac{|X|}{2} + \frac{|Y^c|}{2} $$

**Applying the Expansion Theorem ($|Y^c| \le k - |X| - 1$):**
$$ \text{Peak}_X \le \frac{k}{2} + \frac{|X|}{2} + \frac{k - |X| - 1}{2} $$
$$ \text{Peak}_X \le \frac{k}{2} + \frac{|X|}{2} + \frac{k}{2} - \frac{|X|}{2} - 0.5 $$
$$ \text{Peak}_X \le k - 0.5 $$

Since the mode must be an integer, **$\text{Peak}_X \le k-1$**.

**GAP (unresolved):**
Every term $P_X(x)$ in the decomposition has a mode strictly less than $k$.
However, the claim that the sum inherits this mode bound is NOT proved.
Sums of unimodal/log-concave sequences can fail to be unimodal, so
this step requires a rigorous argument (e.g., showing the sum is LC,
or using interlacing/real-rootedness). No such argument is available.

---

## Final Result

1.  **Lower Bound:** $|S| = k$.
2.  **Upper Bound:** $\text{mode}(I(T)) \le k$.
3.  **Result:** $|S| \ge \text{mode}(I(T))$.

The **1-Private Mode Conjecture** remains OPEN due to the gap in Step 3.
The lower bound (Step 2) and the decomposition are valid, but the
summation step is not proved.
