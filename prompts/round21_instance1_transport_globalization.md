# Round 21 (Instance 1): Globalize the Transport Bound (Local Version Falsified)

You previously proposed:

- global theorem candidate: `sum_err <= D + lambda0*(C10+C01+C11)` with `lambda0=0.05201381704686925`,
- missing local lemma: per odd diagonal `err_{2t+1} <= lambda0*(C10(t)+C01(t)+C11(t))`.

Use only your prior output + this prompt.

## Locked context updates (new scans; treat as facts)

Conventions unchanged: support-root, step-prefix (`k<mode(I_new)`), boundary-correct indexing.

1. Global bound still holds through `n<=18`:
   - `sum_err <= D + lambda0*(C10+C01+C11)` with same `lambda0`.

2. Your local odd lemma is false in this corpus:
   - odd-diagonal checks on `X<0` cases: `587,674`
   - failures of
     `err_{2t+1} <= lambda0*Lambda_old(t)*(10+01+11)`: `38,642` (`6.575%`)
   - worst ratio `err/rhs`: `15.2138`.

3. A stronger local variant is still not enough:
   - with `(Lambda_old(t)+Lambda_old(t+1))*(10+01+11)` at same `lambda0`, failures drop to `1,600` but do not vanish.

Therefore the next step must be **global odd-budget control**, not per-diagonal control.

## Task

Replace the failed per-diagonal lemma with one aggregate theorem that can still imply the observed global inequality.

1. Define an explicit odd aggregate functional `OddDebt(k)` from `{err_{2t+1}}`.
2. State a new missing lemma in global form, e.g.
   - `OddDebt(k) <= lambda0*(C10+C01+C11)`
   - or equivalent weighted cumulative form.
3. Give a proof skeleton that uses diagonal-Abel transport but allows compensation between odd diagonals.
4. Explain why this avoids the local counterexamples above.

## Output format

1. `Revised theorem statement`
2. `Old local lemma failure mode`
3. `New global lemma (exact inequality)`
4. `Proof skeleton with compensation mechanism`
5. `Minimal falsification checklist`

No claims that per-odd-diagonal bound holds.
