# Round 24 (Instance 2): Minimal Extra Channel Design Under n<=19 Data

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Base inequality shape:

`sum_err <= D + lambda0*R_shift + c*Extra`,

with `lambda0 = 0.05201381704686925`, `R_shift = C10+C01+C11`.

From exhaustive `n<=19` X<0 data (`428,434` cases):
- Baseline (no Extra): needs `lambda = 0.08144365672607116` on `R_shift`.

Candidate `Extra` families (all with `C_ab(k)=sum_i Lambda_old(i)P(k-i-a)Q(k-i-b)`):
- `Extra = C20+C02+C21+C12` needs `c = 0.06853082706728619`.
- `Extra = C20+C02+C21` needs `c = 0.08264069901816704`.
- `Extra = C20+C02` needs `c = 0.1048066620468871`.

Additional falsification facts at fixed `c=lambda0`:
- `Extra_axis = C20+C02` fails (`12` cases; max ratio `2.014977...`).
- `Extra_ring = C20+C02+C21+C12+C22` fails (`4` cases; max ratio `1.235550...`).

Common worst witness for top extras:
- `n=19`, `g6=R?????????????????C??w?A^_?_~?`, `root=0`, `step=2`, `k=5`, `(a,b)=(3,14)`.

## Task

Give a principled minimality argument for Extra-channel choice.

You must return:
1. A **primary** Extra formula you think is best (minimal complexity subject to robustness).
2. A **fallback** Extra formula that is strictly stronger but still interpretable.
3. For each, provide a constant target (`c_primary`, `c_fallback`) and explain why these are plausible near-optimal.
4. One explicit "cannot-do-better" obstruction argument from the common worst witness.

Do not just list channels; argue why each shift is needed or redundant.

## Output format

1. `Primary Extra formula and constant`
2. `Fallback Extra formula and constant`
3. `Why each shift is necessary (or removable)`
4. `Witness-based lower bound argument`
5. `Finite verification recipe`
