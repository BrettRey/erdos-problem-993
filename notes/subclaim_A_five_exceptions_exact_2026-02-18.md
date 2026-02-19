# Sub-claim A finite exceptions: exact closure record (2026-02-18)

When using the lane-2 helper surrogate

- `lb_u2(k,j) := mu_{k,j+2}(u2_{k,j}) - mu_{k,j}(lambda_{k,j}) >= 1`,

the only observed failures are the five pairs:

- `(k,j) = (6,0), (6,2), (7,1), (8,0), (8,2)`.

Sub-claim A itself needs

- `Delta_total(k,j) := mu_{k,j+2}(lambda_{k,j+2}) - mu_{k,j}(lambda_{k,j}) >= 1`

(equivalently margin-step `margin(k,j+2)-margin(k,j) >= 0` when mode-step is `+1`).

## Exact-rational verification on the five pairs

Script:

- `prove_subclaim_A_five_exceptions_exact.py`

Run:

- `python3 prove_subclaim_A_five_exceptions_exact.py --out results/subclaimA_five_exceptions_exact_2026-02-18.json`

Output:

- mode-step failures: `0`,
- `Delta_total > 1` failures: `0`.

Exact values:

- `(6,0)`: `Delta_total = 1.0064935539379898...`,
- `(6,2)`: `Delta_total = 1.0048412183126203...`,
- `(7,1)`: `Delta_total = 1.0051330265619831...`,
- `(8,0)`: `Delta_total = 1.0037541618244004...`,
- `(8,2)`: `Delta_total = 1.0053181506584281...`.

Minimum margin above `1` on this finite set:

- `min(Delta_total-1) = 0.0037541618244004...` at `(k,j)=(8,0)`.

So these five helper exceptions do not obstruct Sub-claim A: each satisfies the
actual target inequality strictly.
