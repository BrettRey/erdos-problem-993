# Landslide Strategy: A Physics-Inspired Search for Unimodality Counterexamples

## The Analogy

We treat the sequence of independent set counts, $I(T) = (i_0, i_1, \dots, i_{\alpha})$, as a physical terrain or a pile of granular material.

### 1. The Terrain (The Sequence)
- **Elevation:** The magnitude of $i_k$.
- **Slope:** The ratio $r_k = i_k / i_{k-1}$.
- **Unimodality:** A single peak. The slope $r_k$ should start $>1$, cross $1$ exactly once (the peak), and then stay $\le 1$.
- **Landslide (Counterexample):** A disruption in this profile. Specifically, a localized "collapse" or "bulge" that creates a second peak or a valley (a sequence like ... High, Low, High ...).

### 2. The Physics of Failure
Real landslides occur when:
1.  **Accumulation:** Material builds up, steepening the slope.
2.  **Critical Angle (Angle of Repose):** The steepest angle a pile can maintain before collapsing. In our case, this is the "Log-Concavity" or "Unimodality" limit.
3.  **Trigger:** A small perturbation (rain, earthquake, erosion) that pushes the stress past the strength.
4.  **Propagation:** The failure spreads.

## The Mathematical Mapping

| Landslide Concept | Graph Theory Equivalent |
| :--- | :--- |
| **Material** | Independent sets of a specific size $k$. |
| **Adding Mass** | Grafting subtrees that are "heavy" in specific coefficients. |
| **Steep Slope** | A rapid decrease in coefficients ($i_k \gg i_{k+1}$). |
| **Weak Layer** | A range of $k$ where the sequence is nearly flat ($i_k \approx i_{k+1}$), susceptible to becoming a valley. |
| **Trigger** | Local tree operations: Edge subdivision, leaf movement, grafting a "heavy" branch at a specific node. |

## The Search Algorithm: "Slope Destabilization"

Instead of random evolution, we will build a "Landslide Generator":

### Phase 1: Pile Construction (Accumulation)
Build a "base slope" that is **marginally stable**.
- Use "Broom" graphs or similar structures known to have steep descents and "fat" tails.
- We want a sequence that descends, but *slowly* or *irregularly* in a specific region, creating a "terrace" on the mountainside.

### Phase 2: Debris Loading (The Trigger)
Identify a target index $k_{target}$ on the "terrace".
- We want to *inflate* $i_{k_{target}}$ relative to its neighbors $i_{k_{target}-1}$ and $i_{k_{target}+1}$.
- **Technique:** Graft a specific subtree $S$ to a vertex $v$ in the base tree $T$.
- The new polynomial is roughly the convolution of the base polynomial and the graft polynomial (modulo vertex constraints).
- We select $S$ such that it has a "mass concentration" that aligns with $k_{target}$.

### Phase 3: The Slide (Verification)
Check if the "terrace" has buckled into a "valley" (local minimum).

## Concrete Plan

1.  **Analyze Brooms:** Confirm they act as "steep slopes".
2.  **Catalog "Debris":** Generate small trees ($n \approx 10-15$) and profile their "spectral mass" (where are their coefficients largest?).
3.  **Implement `landslide_search.py`:**
    - **Step A:** Generate a "Base" (e.g., Broom-like or large random tree).
    - **Step B:** Find a "Weak Layer" (indices where $i_k \approx i_{k+1}$).
    - **Step C:** "Bombard" the weak layer. Try grafting "Debris" trees at various points to see if we can push *up* the coeff at $k+1$ while pushing *down* (or keeping steady) the coeff at $k$.

## Why this is different
Previous approaches (Genetic Algorithms, MCMC) search blindly or follow smooth gradients. This approach is **structural engineering**: it deliberately constructs a specific flaw (a flat shelf) and then applies a precise load to break it.
