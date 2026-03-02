# Round 20 (Instance 3): Build the Induction Bridge with a Step-2 Exception

You previously gave the clean induction interface and identified the adjacent-ratio type obstruction.

Use only this prompt + your prior output. Do not claim fresh computation.

## Locked context

- Goal remains chain STP2 at all rootings.
- Empirically, all `X_k<0` occur at step `t=2`; none at `t>=3` through n<=18.
- Global closure uses shifted reserve:
  - `sum_err <= D + lambda*(C10+C01+C11)`, `lambda*=0.05201381704686925`.
- Hard step-2 classes: `(2,14),(3,13),(2,13)`.

Definitions fixed as in prior prompt:

`D`, `C10`, `C01`, `C11`, `err_s`, `sum_err`.

## Task

Produce an induction package with exactly two branches:

1. `t>=3` branch: state the strongest claim you think is actually derivable from your prior obstruction analysis (either strict `X_k>=0` or weaker reserve-free bound).
2. `t=2` branch: state explicit reserve theorem using shifted channels and pair-class constants.

Then isolate the **single bridge lemma** needed so these two branches imply STP2 closure without circularity.

Your bridge lemma must be one of:

- adjacent-ratio/adjacent-minor control on accumulators, or
- equivalent diagonal-centre monotonicity statement.

## Output format

1. `Induction theorem (two-branch statement)`
2. `Bridge lemma (exact inequality)`
3. `Why bridge lemma is not implied by ladder STP2 alone`
4. `Minimal finite verification plan for bridge lemma`

Avoid broad discussion; give a proof-ready skeleton.
