# Attack 1: Mode Superadditivity for PLC Products

Date: 2026-02-19

## Goal

Test whether, for positive log-concave (PLC) coefficient sequences,

`mode(f*g) >= mode(f) + mode(g)`

(where `*` is polynomial multiplication/convolution, and mode is the **leftmost** max index).

If true, this would imply:

`mode(P) >= sum_i mode(I(T_{c_i}))`

for `P = prod_i I(T_{c_i})`, and potentially reduce the bridge target to
`sum_i mode(I(T_{c_i})) >= m-1`.

---

## New Scripts

- `attack1_mode_superadditivity_plc_2026_02_19.py`
  - PLC family checks (simple/binomial/adversarial)
  - tree IS polynomial pair checks
- `attack1_bridge_mode_sum_scan_2026_02_19.py`
  - full `d_leaf<=1`, canonical degree-2 bridge scan through `n<=23`
  - checks:
    1. `mode(P) >= m-1`
    2. `sum_i mode(I(T_{c_i})) >= m-1`
    3. `mode(P) >= sum_i mode(I(T_{c_i}))`

Outputs:

- `results/attack1_mode_superadditivity_plc_2026_02_19.json`
- `results/attack1_bridge_mode_sum_scan_n23_2026_02_19.json`

---

## Commands Run

```bash
python3 attack1_mode_superadditivity_plc_2026_02_19.py
python3 attack1_bridge_mode_sum_scan_2026_02_19.py --min-n 4 --max-n 23 --no-resume
```

---

## Results

## 1) Global PLC superadditivity is false

Explicit PLC counterexample:

- `f = g = [1,2,3,4]` (positive LC)
- `mode(f)=3`, `mode(g)=3`
- `f*g = [1,4,10,20,25,24,16]`, `mode(f*g)=4`
- so `4 < 3+3=6` (deficit `2`).

From adversarial PLC scan:

- checked pairs: `261,425`
- left-mode failures: `85,727`
- max left-mode deficit observed: `9`

Worst witness (left mode):

- `f=g=[22,49,61,70,71,72,73,74,75,76,77]`
- `mode(f)=10`, `mode(g)=10`, `mode(f*g)=11`
- deficit `=20-11=9`.

Sanity check family:

- Binomials `(1+x)^n` (left mode) showed **0 failures** in tested range (`n<=40`).

---

## 2) Superadditivity also fails for tree IS polynomials

Exhaustive unique tree IS polynomials (`n<=12`):

- checked pairs: `382,375`
- left-mode failures: `56,562`
- max left-mode deficit: `1`

First failure:

- edge tree polynomial `f=g=[1,2]` (graph6 `A_`)
- `f*g=[1,4,4]`
- `mode(f*g)=1 < 2=mode(f)+mode(g)`.

Random tree-pair stress (`250,000` pairs from trees up to `n<=16`) also shows many failures (`22,338` left-mode failures).

---

## 3) Bridge decomposition consequences (`d_leaf<=1`, full `n<=23`)

Total checked: `931,596` trees.

1. `mode(P) >= m-1`:
   - failures: `0` (still fully holds computationally).

2. `sum_i mode(I(T_{c_i})) >= m-1`:
   - failures: `32,181`.
   - minimum margin: `-1`.
   - witness: graph6 `G?`DE_` (`n=8`), `m=3`, child-mode sum `1 < 2=m-1`.

3. `mode(P) >= sum_i mode(I(T_{c_i}))`:
   - failures: `31,124`.
   - minimum margin: `-3`.
   - witness: graph6 `R??????_A?C?C?A??_?C??O??_F~??` (`n=19`), `k=8`, child modes all `1`, sum `8`, but `mode(P)=5`.

Split by number of children at `u`:

- `k=1`: no failures for (2) or (3)
- `k>=2`: almost all failures occur here (dominant at `k=2,3`).

---

## Conclusion

The proposed mode-superadditivity route is not viable:

- `mode(f*g) >= mode(f)+mode(g)` is false for PLC in general.
- It is also false for tree IS polynomials.
- In the bridge setting, both derived inequalities
  `mode(P) >= sum_i mode(I(T_{c_i}))` and
  `sum_i mode(I(T_{c_i})) >= m-1`
  fail frequently.

What remains true (computationally through full `n<=23`) is the direct target:

`mode(P) >= m-1`.

So Attack 1 should be treated as a negative route-closure result, not a proof path.
