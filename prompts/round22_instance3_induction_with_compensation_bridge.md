# Round 22 (Instance 3): Build Induction Bridge from Compensation, Not Adjacent Minors

Use only your prior outputs plus this prompt.

## Locked falsifications and updates

- BL-adj bridge (`E(r+1)J(r) <= E(r)J(r+1)`) is false in the corpus.
- old local odd bridge forms are false.
- new candidate invariant from exhaustive diagnostics:
  local/slot-level neighbor compensation `d_t <= e_t + e_{t+1}` has 0 failures through `n<=18`.

Definitions:

- `d_t := max(0, err_{2t+1} - lambda0*Lambda_t*(10+01+11))`
- `e_t := max(0, Lambda_t*P(k-t)Q(k-t) - err_{2t})`

## Task

Give a two-branch induction package where the bridge is a compensation invariant rather than an adjacent-minor monotonicity.

1. Step-2 branch: reserve closure (as before).
2. t>=3 branch: state bridge invariant in cumulative/neighbor-compensation form.
3. Show non-circular dependency graph: which lemmas depend on which IH pieces.
4. State one minimal missing lemma in this new framework.

## Output format

1. `Two-branch induction theorem`
2. `Compensation bridge invariant`
3. `Dependency graph (no circularity)`
4. `Single missing lemma`
5. `Finite validation checklist`

Do not use adjacent-minor bridge assumptions.
