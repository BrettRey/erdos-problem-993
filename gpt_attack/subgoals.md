# Sub-goals — pick exactly one

Each sub-goal is stated with (a) the formal target, (b) what closing it would buy, (c) the most likely obstruction, and (d) suggested first moves. Start with whichever you think has the best ratio of tractability to payoff.

---

## SG1. Prove Conjecture A in the reduced regime

**Target:** For every tree T on $n$ vertices with $d_{\mathrm{leaf}}(v)\leq 1$ for every vertex v,
$$\mathrm{mode}(I(T)) \leq \lfloor n/3 \rfloor + 1.$$

**Pay-off:** Closes the Conjecture A component of the PNP route. To close #993 from here, this still has to be paired with the formal PNP-to-mode reduction and the Case B hub bound/reduction stated in the manuscript.

**Current status:** Verified for all $d_{\mathrm{leaf}}\leq 1$ trees through $n = 23$. Two partial routes exist:
- Mean bound (proved for $d_{\mathrm{leaf}} \leq 1$) + Conjecture 3.6 (mode ≤ ⌈μ⌉, conjectural): would give Conjecture A for the reduced class.
- Reduction from all trees to $d_{\mathrm{leaf}} \leq 1$ trees via Hub Exclusion + Transfer lemmas: currently informal in the manuscript; needs a formal reduction theorem.

**Most likely obstruction:** The jump from empirical verification to a proof requires closing Conjecture 3.6 (mode-mean) at least for $d_{\mathrm{leaf}}\leq 1$ trees, proving the tie-fugacity condition, or finding a direct mode bound in the reduced regime.

**Suggested first moves:**
- Attempt Conjecture 3.6 via the subdivision-contraction identity: if $\mathrm{mode}(I(T)) \leq \lceil \mu(T) \rceil$ holds for T and $T/e$, what does the identity say about $T_e$?
- Alternative: attack the mode bound directly using the proved $\mu(T)<n/3$ theorem for $d_{\mathrm{leaf}}\leq 1$ trees.

---

## SG2. Prove ECMS (Edge Contraction Mode Stability)

**Target:** A precise theorem stating that for a designated edge class in every tree T, mode positions behave stably under contraction. (See manuscript §4 for the exact statement.)

**Pay-off:** Combined with the combined tail condition (also conjectural), subdivision preserves unimodality, so every minimal counterexample is homeomorphically irreducible. This reduces the search space dramatically because homeomorphically irreducible trees are sparse — there are only 1, 1, 1, 2, 2, 4, 6, 11, 18, 37, ... of them up to n = 10, vs. 106 general trees on 10 vertices.

**Current status:** Verified for 24.7 million edges. No proof.

**Most likely obstruction:** The identity $I(T_e; x) = I(T; x) + x I(T/e; x)$ makes the mode of $I(T_e)$ depend on the interaction between the mode of $I(T)$ and the shifted mode of $I(T/e)$. The required stability property is a delicate inequality between these two.

**Suggested first moves:**
- Write the ECMS claim as an explicit coefficient inequality. Use the subdivision-contraction identity to express $i_k(T_e)$ in terms of $i_k(T)$ and $i_{k-1}(T/e)$. Then the mode stability claim becomes a comparison of the index achieving the maximum in $i_k(T) + i_{k-1}(T/e)$ vs. $i_k(T)$.
- Examine the 24.7M empirical cases for a structural pattern. Are there edges where ECMS is "barely" satisfied? Those are the stress test for any candidate proof.

---

## SG3. Prove the mode-mean or tie-fugacity bridge

**Target:** Prove one of the following for all trees, or at least for all $d_{\mathrm{leaf}}\leq 1$ trees:
- $\mathrm{mode}(I(T)) \leq \lceil \mu(T) \rceil$.
- The manuscript's single-tie condition $\mu(\lambda_m)\geq m-1$, where $m$ is the leftmost mode and $\lambda_m=i_{m-1}/i_m$.

**Pay-off:** Combined with the proved mean bound $\mu(T)<n/3$ for $d_{\mathrm{leaf}}\leq 1$ trees, this proves Conjecture A in the reduced regime.

**Current status:** $\mathrm{mode}(I(T)) \leq \lceil\mu(T)\rceil$ is verified for all trees through $n=23$; the tie-fugacity condition is verified for all $d_{\mathrm{leaf}}\leq 1$ trees through $n=23$ and all trees through $n=22$.

**Most likely obstruction:** Darroch-style mode localization is known for Poisson-binomial / real-rooted generating polynomials, but tree independence polynomials are not real-rooted in general. The proof must exploit tree DP structure, not generic log-concavity.

**Suggested first moves:**
- Try to prove the single-tie condition directly from coefficient comparisons at $\lambda_m$.
- Check whether subdivision-contraction preserves the mode-mean bound under a controlled hypothesis.
- Do not try to prove $\mu(T)<n/3$ for all trees. This is false: for $K_{1,d}$, $\mu=(d2^{d-1}+1)/(2^d+1)\sim d/2$, while $n/3=(d+1)/3$.

---

## SG4. Prove a quantitative dip-threshold inequality

**Target:** Any tree T violating unimodality at index $k^*$ (i.e., $i_{k^*-1} > i_{k^*} < i_{k^*+1}$) satisfies
$$\frac{i_{k^*}^2}{i_{k^*-1} \cdot i_{k^*+1}} < 1 - f(n, k^*)$$
for some explicit function $f$ with $f(n, k^*) > 0$ for all relevant $(n, k^*)$.

**Pay-off:** A non-trivial lower bound on the dip depth, combined with known upper bounds on coefficient ratios (Levit 2006, Chudnovsky–Seymour root-location, the spider log-concavity results), may produce a contradiction that closes the conjecture.

**Current status:** No quantitative dip-threshold known.

**Most likely obstruction:** The coefficient sequence of a tree is constrained by the tree DP recurrence but no existing ratio inequality directly handles the unimodality-failure case. The ~35,000 Ramos–Sun log-concavity failures at $n = 60$ give empirical lower bounds on how sharp log-concavity can fail; no such data exists for unimodality because none has been observed.

**Suggested first moves:**
- Find the 19 log-concavity-failing trees at $n = 28$ in the reproducibility artifact. Compute their ratios $i_k^2 / (i_{k-1} i_{k+1})$ exactly. Do the same for the $n = 60$ Ramos–Sun families. Is there a pattern as $n \to \infty$?
- The subdivision-contraction identity implies an inequality for the coefficient ratios of $T_e$ in terms of T and $T/e$. Does this constrain how the ratios can fall below 1?

---

## SG5. DRC analog for independence-set hypergraphs (speculative)

**Target:** Construct a dependent-random-choice-style extraction argument for the $k$-independent-set hypergraph $H_k(T)$ whose vertices are $V(T)$ and whose hyperedges are the size-$k$ independent sets. Identify a structural consequence for $i_{k-1}, i_k, i_{k+1}$ if the mode is not at $k$.

**Pay-off:** If a DRC analog produces a common-neighbourhood bound that compares $i_k$ unfavourably to $i_{k \pm 1}$ in a mode-failing configuration, a contradiction closes the conjecture.

**Current status:** Speculative. No one has pointed DRC at the independence-polynomial question.

**Most likely obstruction:** DRC requires density (high minimum degree); trees are sparse. The $H_k(T)$ hypergraph may or may not be "dense" in the right combinatorial sense. Most likely this dies at the first density check.

**Suggested first moves:**
- Compute, for the 19 LC-failing trees at $n = 28$, the "density" of the $H_{14}(T)$ hypergraph: what fraction of $\binom{n}{14}$ are independent? How does it compare to non-failing trees?
- Read Fox–Sudakov (2011) on DRC in sparse settings. DRC has extensions to sparse graphs; see Conlon–Fox–Sudakov for the relevant techniques.

**First-pass density check (2026-05-20):** `python3 gpt_attack/sanity_check.py` finds that the n=28 LC failures are extremely sparse at the failure index:
$i_{14}/\binom{28}{14}$ ranges from $1/786600$ to $1/422280$.
The top n=28 near-miss trees are much denser at their near-miss index:
$i_{13}/\binom{28}{13}$ ranges from $9296/66861$ to $145609/832048$ among the top 50.
This makes a direct DRC attack on LC-failure hypergraphs unattractive; if SG5 is pursued, target near-miss structure rather than LC-failure structure.

**Note:** This is the highest-risk, highest-reward sub-goal. It's included because the #1014 proof showed that DRC can be pointed at unexpected problems and produce elementary arguments. If the first density check fails, abandon this sub-goal and move to SG1–SG4.

---

## How to choose

If you are a frontier model doing this for the first time, SG1 or SG3 are the most direct paths to a publishable contribution. SG2 is the cleanest structural attack. SG4 is the most novel and likely to fail. SG5 is a long-shot methodological probe.

**Pick one. State your pick in the first two lines of your response.**
