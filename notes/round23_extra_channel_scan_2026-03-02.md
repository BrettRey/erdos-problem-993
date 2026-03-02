# Round 23 Local Extra-Channel Scan (2026-03-02)

Purpose: identify a minimal explicit `Extra(k)` channel family for the corrected inequality

`sum_err <= D + lambda0*R_shift + c*Extra`,

using locked conventions (support-root, boundary-correct indexing, step-prefix `k < mode(I_new)`, exhaustive trees up to `n<=19`).

## Constants

- `lambda0 = 0.05201381704686925`
- Baseline required on old reserve only:
  - `sum_err <= D + lambda*R_shift` needs
  - `lambda = 0.08144365672607116`
  - witness: `n=19`, class `(a,b)=(3,14)`, step `2`, `k=5`.

## New result (full X<0 corpus, n<=19)

Scanned all `X<0` step-prefix instances (`428,434` total under this scan definition) and measured

`deficit := (sum_err - D - lambda0*R_shift)_+`.

Tested candidate extras of the form `Extra = sum C_ab` where

`C_ab(k) = sum_i Lambda_old(i) P(k-i-a)Q(k-i-b)`.

### Ranked candidates (smaller `needed_c` is better)

- `Extra = C20 + C02 + C21 + C12`
  - `needed_c = 0.06853082706728619`
- `Extra = C20 + C02 + C21`
  - `needed_c = 0.08264069901816704`
- `Extra = C20 + C02 + C12`
  - `needed_c = 0.08310634062872929`
- `Extra = C20 + C02 + C30`
  - `needed_c = 0.09364289753090758`
- `Extra = C20 + C02`
  - `needed_c = 0.1048066620468871`

Single channels are much weaker:

- `C20`: `0.20453408622242417`
- `C02`: `0.2149512536697084`

### Common worst witness

For all top candidates above, worst case is the same:

- `n=19`, `g6=R?????????????????C??w?A^_?_~?`, `root=0`, `step=2`, `k=5`, `(a,b)=(3,14)`.

## Practical theorem candidates

Two concrete options now have exhaustive support through `n<=19`:

1. Keep old reserve only (single-constant upgrade):

`sum_err <= D + 0.08144365672607116 * (C10+C01+C11)`.

2. Keep `lambda0` and add one extra family:

`sum_err <= D + lambda0*(C10+C01+C11) + 0.06853082706728619*(C20+C02+C21+C12)`.

The second option preserves the original `lambda0` narrative but adds an explicit correction channel.

## Artifacts

- `results/round23_extra_channel_fit_n19_full.json`
- `results/round23_extra_channel_miniset_n19_full.json`

