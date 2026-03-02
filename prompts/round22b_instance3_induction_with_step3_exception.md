# Round 22b (Instance 3): Induction Bridge with Step-3 Exception Regime

Use only your prior outputs plus this prompt.

## Locked context update

- BL-adj bridge is false (already known).
- At `n=19`, negatives are no longer only step 2:
  - `xneg_step_not2 = 12`
  - all observed non-step2 negatives are at `step=3`, `k=1`, pair `(a,b)=(2,2)` in the current witness set.
- fixed-`lambda0` global reserve has `33` failures at `n=19`.

## Task

Produce a non-circular induction theorem with explicit exception handling:

1. Step-2 branch: reserve closure (possibly upgraded constant/channel).
2. Step>=3 branch: state a bridge invariant that is true except for a sharply defined micro-regime.
3. Add a third micro-regime branch for the observed step-3 exception pattern and show how to close it.

Avoid adjacent-minor monotonicity assumptions.

## Output format

1. `Three-branch induction theorem`
2. `Bridge invariant for generic step>=3`
3. `Micro-regime exception lemma`
4. `Dependency graph (no circularity)`
5. `Minimal finite verification plan`
