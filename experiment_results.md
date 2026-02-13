# Evolutionary Search Experiment Results

## Objective
To find a counterexample to the unimodality conjecture by evolving trees that maximize the log-concavity violation ratio:
$$ R = \max_k \frac{i_{k-1} i_{k+1}}{i_k^2} $$
The hypothesis is that maximizing this ratio locally might eventually break the unimodality constraint ($i_{k-1} > i_k < i_{k+1}$).

## Setup
- **Algorithm:** Genetic Algorithm with elitism and tournament selection.
- **Representation:** NetworkX Graphs (converted to adjacency lists for polynomial computation).
- **Population:** 100-200 individuals.
- **Seeds:** The 2 known log-concavity failures at $N=26$ and top "near misses".
- **Mutations:** Leaf addition, subtree grafting, leaf rewiring, broom feature injection.
- **Constraint:** We experimented with disabling leaf removal to force tree growth ($N > 26$).

## Observations
1. **Initial Seed Dominance:** The seed tree at $N=26$ with ratio $\approx 1.145$ (at the tail) is extremely dominant. Random mutations tend to destroy this specific structural feature (a specific sparse tail), dropping the fitness immediately.
2. **Fitness Landscape:** The landscape appears to be a "needle on a plateau" or very rugged. Small changes to the high-fitness tree result in drastic fitness drops.
3. **Performance:** Python-based polynomial computation limits the search speed to ~20 generations per minute for population size 100. This is too slow for deep search without parallelization or C++ optimization.
4. **Growth:** When forced to grow, the LC ratio tends to decrease, suggesting the $N=26$ failure might be a "peak" anomaly for small $N$, or that the failure mode requires precise "tuning" that is harder to maintain as $N$ grows.

## Next Steps
- **Optimize:** Implement the polynomial computation in C or optimize the Python version (profiling needed).
- **Diversity:** Use a fitness sharing or novelty search approach to maintain population diversity and stop the $N=26$ seed from dominating.
- **Target:** Instead of global max ratio, target specific indices (e.g., "maximize $i_{10}$ while minimizing $i_{11}$").
