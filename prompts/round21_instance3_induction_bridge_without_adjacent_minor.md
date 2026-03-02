# Round 21 (Instance 3): Rebuild the Bridge (BL-adj is Empirically False)

You proposed bridge lemma:
`E(r+1)J(r) <= E(r)J(r+1)` on step-2 accumulators.

Use only your prior output + this prompt.

## Locked falsification (new exhaustive scan)

Across all step-2 support-root accumulators through `n<=18`:

- instances checked: `410,832`
- violations of BL-adj: `1,221,519`
- min slack of `E(r)J(r+1)-E(r+1)J(r)`: `-1550`

So BL-adj cannot be the bridge.

## Still-locked facts

- `X<0` only at `t=2` through `n<=18`
- step-2 global reserve closure holds with `lambda0` and shifted channels.

## Task

Produce a new two-branch induction package that does **not** use BL-adj.

1. Keep `t=2` reserve branch.
2. Replace `t>=3` branch with a weaker, plausibly true bridge condition that is consistent with BL-adj failures.

The new bridge should be one of:

- aggregated odd-budget monotonicity across steps,
- a bound on increment of odd debt from `t` to `t+1`,
- or a cumulative reserve invariant.

## Output format

1. `Revised two-branch induction theorem`
2. `New bridge invariant (exact inequality)`
3. `Why it survives BL-adj counterexamples`
4. `How to test bridge invariant on finite data`

Do not propose any adjacent-minor monotonicity as required assumption.
