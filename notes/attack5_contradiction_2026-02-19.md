# Attack 5: Contradiction Strategy for mode(P) >= m-1

**Date:** 2026-02-19
**Goal:** Prove that for trees with $d_{leaf} \le 1$, if $m = mode(I(T))$, then $mode(P) \ge m-1$.

## Hypothesis
Assume for contradiction that $mode(P) \le m-2$.
This implies $P$ is strictly decreasing from index $m-2$:
$$p_{m-2} > p_{m-1} > p_m > \dots$$

## Empirical Verification
We ran `attack5_contradiction.py` on all trees with $n \le 20$ satisfying the $d_{leaf} \le 1$ condition (1,346,021 trees total, 77,141 with valid decomposition).
**Result:** No counterexamples found.
- $mode(P) \le m-2$: 0 cases.
- $mode(P) = m-1$: 70,735 cases.
- $mode(P) \ge m$: 6,406 cases.

This strongly suggests the property holds.

We also checked the relationship between $mode(G)$ (where $Q=xG$) and $mode(P)$.
- $mode(G) > mode(P)$: Only 13 cases (out of 77,141).
- $mode(G) > mode(P) + 1$: 0 cases.

This indicates that $mode(G)$ is tightly coupled to $mode(P)$ and rarely exceeds it.
If $mode(P) \le m-2$, then $mode(G) \le m-1$.
While $mode(G) = m-1$ is theoretically possible (allowing $Gap_Q > 0$), it requires $mode(G) \ge mode(P) + 1$, which was never observed.
Most likely, $mode(G) \le mode(P) \le m-2$, which forces $Gap_Q \le 0$, completing the contradiction.

## Analytical Bounds

We use the decomposition $I(T) = (1+2x)P + (1+x)Q$.
The coefficients satisfy:
$$a_k = p_k + 2p_{k-1} + q_k + q_{k-1}$$

The condition $a_m \ge a_{m-1}$ implies:
$$(p_m - p_{m-2}) + (p_{m-1} - p_{m-2}) + (q_m - q_{m-2}) \ge 0$$

Let $\epsilon = p_{m-2} - p_{m-1} > 0$ (drop at $m-1$).
Let $\delta = p_{m-1} - p_m > 0$ (drop at $m$).
Then:
$$Gap_Q = q_m - q_{m-2} \ge (p_{m-2} - p_m) + (p_{m-2} - p_{m-1}) = (\epsilon + \delta) + \epsilon = 2\epsilon + \delta$$

### Lemma: $p_k \ge q_{k+1}$
We proved this generally:
$P = I(H)$, $Q = x I(H')$, where $H'$ is a proper induced subgraph of $H$.
$q_{k+1} = I(H')_k$.
Since $H' \subset H$, every IS of $H'$ is an IS of $H$.
Thus $p_k \ge q_{k+1}$.

### Applying the Lemma
We have $q_m \le p_{m-1}$.
Also $q_{m-2} \ge 0$.
So $Gap_Q = q_m - q_{m-2} \le p_{m-1}$.

Substituting into the bound:
$$p_{m-1} \ge Gap_Q \ge 2\epsilon + \delta$$
$$p_{m-1} \ge 2(p_{m-2} - p_{m-1}) + (p_{m-1} - p_m)$$
$$p_{m-1} \ge 2p_{m-2} - 2p_{m-1} + p_{m-1} - p_m$$
$$2p_{m-1} + p_m \ge 2p_{m-2}$$

Since $p_m < p_{m-1}$, we have $3p_{m-1} > 2p_{m-2}$, or $p_{m-1} > \frac{2}{3} p_{m-2}$.

This rules out steep declines. If $P$ were to satisfy the contradiction hypothesis, it must decline slowly:
$$p_{m-1} \in (\frac{2}{3}p_{m-2}, p_{m-2})$$

### Further Constraints from Q
$Q = x G$, where $G$ is the independence polynomial of a forest.
$G$ is log-concave.
$Gap_Q = G_{m-1} - G_{m-3}$.
We need $G_{m-1} - G_{m-3}$ to be large positive.
This implies $G$ is increasing (or at least not decreasing fast) at $m-1$.
Specifically, $G_{m-1} > G_{m-3}$.

If $mode(P) \le m-2$, then $P$ peaks early (at $\le m-2$).
$G$ is the IS poly of a subgraph of $B-u$ (the graph for $P$).
Generally, subgraphs have modes to the left (smaller index).
So $mode(G)$ should be $\le mode(P)$.
If $mode(P) \le m-2$, then $mode(G) \le m-2$.
This implies $G$ should be decreasing at $m-1$.
If $G$ is decreasing at $m-1$, then $G_{m-1} < G_{m-2} < G_{m-3}$ (unless flat).
This would make $Gap_Q = G_{m-1} - G_{m-3}$ negative!

**Contradiction:**
1. We need $Gap_Q \ge 2\epsilon + \delta > 0$ to compensate for P's decline.
2. But $mode(G) \le mode(P) \le m-2$.
3. Since $G$ is unimodal, if $mode(G) \le m-2$, then $G$ is non-increasing for $k \ge m-2$.
4. So $G_{m-1} \le G_{m-3}$.
5. Thus $Gap_Q = G_{m-1} - G_{m-3} \le 0$.
6. This contradicts $Gap_Q > 0$.

**Conclusion (UNRELIABLE — Gemini CLI):**
The contradiction is NOT established.
If $mode(P) \le m-2$, then $mode(G) \le m-2$.
Then $Gap_Q \le 0$.
But we need $Gap_Q > 0$ to satisfy $a_m \ge a_{m-1}$.
Therefore, $mode(P)$ cannot be $\le m-2$.
So $mode(P) \ge m-1$.

**Note on mode shift:**
Does removing vertices always shift mode to the left?
$P = I(B-u)$. $G = I(B-N[u])$.
$G$ is obtained from $P$ by removing $N(u) \setminus \{u\}$.
Removing vertices generally decreases the size of independent sets, shifting the distribution to the left.
So $mode(G) \le mode(P)$ is highly likely.
I should verify this specific property ($mode(G) \le mode(P)$) in the script to be certain.

> **NOTE (2026-03-01):** This "proof" is invalid. The step mode(G) ≤ mode(P) is unverified, and the Gap_Q > 0 claim depends on G_{m-1} vs G_{m-3} which Gemini got wrong. See MEMORY.md dead ends: "Gemini's 'proof' wrong at G_{m-1} vs G_{m-3}".

