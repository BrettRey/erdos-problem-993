# What is already known — read before attacking

Sourced from the companion paper (doi:10.5281/zenodo.19100781) and its
verified computations. Build past this; don't re-derive it.

## Proved (unimodality holds) — don't re-prove these classes
- **Paths and centipedes** (Alavi–Malde–Schwenk–Erdős 1987).
- **Regular caterpillars** (Galvin–Hilyard 2018).
- **Fibonacci trees** (Bencs 2018).
- **Spiders** — their independence polynomials are log-concave (Li et al. 2025),
  which covers **brooms** as a special case.
- Two further families containing the 2025 log-concavity-failing trees were
  proved unimodal (2026).
- **Tail is always safe:** `i_k` is strictly decreasing for
  `k >= ceil((2*alpha - 1)/3)` (Levit–Mandrescu 2006). So any valley must sit
  in roughly the first two-thirds of the sequence.
- **Random trees:** for a uniformly random labelled tree, a.a.s. the first
  ~49.5% of the sequence is increasing and the last ~38.8% decreasing
  (Basit–Galvin 2020).

## Proved in the companion paper (positive, but not the full conjecture)
- **Mean bound:** `mean(T) < n/3` for every tree with at most one private leaf
  per support vertex (`dleaf <= 1`), via a decimation identity and Steiner
  peeling. THEOREM. It would give a mode bound for that class *if* the
  mode–mean localization below held.
- **Leaf-attachment asymptotics:** attaching `s` leaves to a fixed vertex of a
  fixed base gives near-miss ratio `nm(s) = 1 - C/s + O(1/s^2)` with a
  parity-dependent constant `4 <= C <= 8` (for pure stars, `C=6` on even `s`,
  `8` on odd). THEOREM. Consequence: every such family is unimodal for all large
  `s`, and its near-miss ratio approaches 1 but the gap shrinks like `C/s` and
  never closes. This is *why* the champions creep toward the wall without
  crossing.

## Conditional hinges (proving/refuting either is high-value = T3)
- **Mode–mean localization (OPEN for trees):** is `mode <= ceil(mean)` (or the
  stronger `mode in {floor(mean), ceil(mean)}`) for tree independence sequences?
  Darroch (1964) proved the analogue for Poisson–binomial distributions.
  Verified for all 931,596 trees with `dleaf <= 1` on `n <= 23`, and for the
  known log-concavity breakers individually. Not proved in general.
- **Edge Contraction Mode Stability (ECMS) + a combined tail condition
  (verified only through n <= 19):** together they imply a subdivision–
  contraction identity forcing any *minimal* counterexample to be
  **homeomorphically irreducible** (no degree-2 vertices). Both are conjectural.
  Discharging ECMS is a clean sub-target.

## The search frontier (where a counterexample must live) — this is EVIDENCE
- **Exhaustive, zero counterexamples:** all **8,691,747,673** trees on
  `n <= 29` vertices checked in exact arithmetic. So any counterexample has
  **n >= 30**.
- **Log-concavity already fails** (but unimodality survives): exactly **2**
  trees at n=26 (edge lists + exact sequences embedded as DATA in erdos993_handoff.py), **0** at n=27, **19** at n=28, plus infinite
  known non-log-concave families. Unimodality has never failed. The gap between
  "log-concave" and "unimodal" is where the action is.
- **Best near-miss among all trees through n=28:** `nm = 0.8571` (at n=27).
  Beyond exhaustive range, tuned **multi-arm stars** `M(s; a1..ak)` are the
  empirically extremal family: `nm = 0.9437` (n=75), `0.9575` (n=100),
  `0.9792` (n=200), `0.9959` (n=1000). All obey the `1 - C/s` law above.

## Reading of the landscape (orientation, not proof)
The small regime is exhausted and the extremal families are provably bounded
away from the wall. So a counterexample, if one exists, is **large** and lies
**outside** every family listed here — or the conjecture is true. Either a new
structural class (T2), a discharge of a conditional hinge (T3), or an extremal
family that provably beats `1 - C/n` (T4) is a real contribution. A blind
random search over trees `n < 30` is wasted effort; that space is closed.
