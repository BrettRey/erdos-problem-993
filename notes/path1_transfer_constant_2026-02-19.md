# Path 1 Transfer Constant (2026-02-19)

Goal: test whether Route-1 (`mu_P(lambda_m(T)) >= m-2`) can be deduced directly from
Route-2-type bounds on `mu_B(lambda_m(T))`, where:

- `B = T-{l,s}` for canonical degree-2 bridge leaf (`l`: min-support-degree rule),
- `P = dp_B[u][0]`,
- `m = mode(I(T))`,
- `lambda = lambda_m(T) = i_{m-1}(T)/i_m(T)`.

## Quantities scanned

Define:

- `D := mu_B(lambda) - mu_P(lambda)`
- `a := lambda/(1+lambda)`
- `half_excess := D - 1/2`
- `exact_excess := D - (1-a)`

If `exact_excess <= 0` always, then the exact Route-2 threshold

`mu_B >= m - 1 - a`

would imply Route-1 immediately:

`mu_P = mu_B - D >= m - 1 - a - (1-a) = m-2`.

## New scanner

- script: `conjecture_a_route1_transfer_scan.py`

Artifacts:

- `results/whnc_route1_transfer_scan_n20.json`
- `results/whnc_route1_transfer_scan_n22.json`
- `results/whnc_route1_transfer_scan_n23_mod8_r0.json` ... `r7.json`
- merged full frontier:
  `results/whnc_route1_transfer_scan_n23_merged.json`

`n=23` was run in nauty split mode (`r/8`) and merged, to avoid partial-run loss.

## Full-frontier result (d_leaf<=1, n<=23, canonical leaf)

From `results/whnc_route1_transfer_scan_n23_merged.json`:

- checked trees: `931,596`
- `max_D = 0.5059567659780138`
- `max_half_excess = 0.005956765978013756`
- `max_exact_excess = 0.005107340588772491`
- `half_excess > 0` count: `10`
- `exact_excess > 0` count: `4`

So the direct implication `mu_B >= m-1-a => mu_P >= m-2` is **not exact**,
but it fails only by a very small additive constant.

Equivalent calibrated statement from data:

`mu_B >= m - 1 - a + kappa*  =>  mu_P >= m - 2`

with empirical `kappa* = 0.005107340588772491` on this frontier.

## Extremal witness for `exact_excess`

- `g6`: `V?????????????_?G?@??C??G??G??E???o?oB?@|S??`
- `n=23`, `m=8`, `lambda=0.9972305185989683`
- `mu_P=6.4302934756640955`
- `mu_B=6.9360941466859956`
- `D=0.5058006710219001`
- `a=0.4993066695668724`
- `1-a=0.5006933304331276`
- `exact_excess = 0.005107340588772491`

## Interpretation for Path 1

This does not close Path 1 algebraically yet, but it narrows it:

1. A Route-2-style proof with a tiny explicit additive slack (`~0.006`) would force Path 1.
2. The transfer gap is quantitatively tiny and now explicitly measured on full `n<=23`.

