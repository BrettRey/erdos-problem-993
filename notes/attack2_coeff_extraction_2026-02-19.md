# Attack 2: Direct Coefficient Extraction from `I(T)=(1+2x)P+(1+x)Q`

Date: 2026-02-19

## 1. Goal and setup

For canonical degree-2 bridge decomposition (`d_leaf <= 1` trees):

- choose leaf `l` with degree-2 support `s`,
- `u` is the other neighbor of `s`,
- `B = T - {l,s}`,
- `P = dp_B[u][0]`, `Q = dp_B[u][1]`,
- `I(T) = (1+2x)P + (1+x)Q`.

Write `p_k=[x^k]P`, `q_k=[x^k]Q`, `a_k=[x^k]I(T)`, and `m=mode(I(T))` (leftmost mode).

Coefficient identity:

`a_k = p_k + 2p_{k-1} + q_k + q_{k-1}`.

Target:

`mode(P) >= m-1`, equivalently `g := p_{m-1}-p_{m-2} >= 0` (since `P` is PLC/unimodal).

## 2. Direct mode inequalities at `m`

From `a_m >= a_{m-1}`:

`I1 := a_m-a_{m-1} = (p_m-p_{m-2}) + (p_{m-1}-p_{m-2}) + (q_m-q_{m-2}) >= 0`.

So

`g >= -[(p_m-p_{m-2}) + (q_m-q_{m-2})]`.

From `a_m >= a_{m+1}`:

`I2 := a_m-a_{m+1} = (2p_{m-1}-p_m-p_{m+1}) + (q_{m-1}-q_{m+1}) >= 0`.

### Contradiction template (`g<0`)

Assume `g<0`. Since `P` is LC with nonnegative coefficients, ratios
`rho_k = p_k/p_{k-1}` are nonincreasing where defined.
Then `rho_{m-1}<1` implies `rho_m <= rho_{m-1}<1`, so `p_m < p_{m-1} < p_{m-2}` and:

`p_m - p_{m-2} <= g(1+rho_{m-1}) < g < 0`.

Hence the P-part of `I1` is strictly negative:

`(p_m-p_{m-2}) + g < 2g < 0`.

So `I1>=0` forces a positive two-step rise in `Q`:

`q_m - q_{m-2} >= -[(p_m-p_{m-2})+g]`.

This is exactly the compensation mechanism to test computationally.

## 3. Full extraction scan (`n<=23`)

Script:

- `attack2_coeff_scan_2026_02_19.py`

Command:

```bash
python3 attack2_coeff_scan_2026_02_19.py --min-n 4 --max-n 23 --out results/attack2_coeff_scan_n23.json
```

### Headline results

- total trees seen (`geng`): `23,942,356`
- `d_leaf<=1` checked: `931,596`
- `mode(P) >= m-1` failures (`g<0`): `0`
- ties (`g=0`): `0`
- minimum `g`: `1`
- `I1` failures: `0`
- `I2` failures: `0`
- `P` not LC: `0`
- `Q` not LC: `0`
- `P` not unimodal: `0`
- `q`-shift identity failures (`q_k=p'_{k-1}`): `0`
- coefficientwise dominance `P >= P'` failures: `0`

### Term-sign / compensation profile

Let

- `p_term1 = (p_m-p_{m-2}) + g`,
- `q_term1 = q_m-q_{m-2}`,
- `p_term2 = 2p_{m-1}-p_m-p_{m+1}`,
- `q_term2 = q_{m-1}-q_{m+1}`.

Observed:

- `p_term1 < 0`: `115,218` cases (`12.37%`)
- `q_term1 < 0`: `359` cases (`0.0385%`)
- `p_term2 < 0`: `57` cases (`0.0061%`)
- `q_term2 < 0`: `208,367` cases (`22.37%`)

Exact compensation counts:

- `p_term1<0` compensated by `q_term1`: `115,218/115,218`
- `q_term1<0` compensated by `p_term1`: `359/359`
- `p_term2<0` compensated by `q_term2`: `57/57`
- `q_term2<0` compensated by `p_term2`: `208,367/208,367`

Tightness:

- minimum compensation ratio on `I1` (`q_term1/(-p_term1)` when `p_term1<0`): `1.001057...`
- minimum compensation ratio on `I2` (`q_term2/(-p_term2)` when `p_term2<0`): `1.0` (exact tie cases)

Largest required `Q` compensation:

- from `I1`: `2766` (witness in `results/attack2_coeff_scan_n23.json`)
- from `I2`: `92`

## 4. Candidate sufficient condition check

Requested candidate: `2p_{m-1} >= a_{m-1}`.

Observed:

- true in only `12 / 931,596` cases (`0.0013%`)
- false in `931,584` cases.

So this is not a useful global sufficient condition.

## 5. Why the two mode inequalities alone do not force `g>=0`

Script:

- `attack2_local_tuple_obstruction_2026_02_19.py`

Command:

```bash
python3 attack2_local_tuple_obstruction_2026_02_19.py --max-val 80 --out results/attack2_local_tuple_obstruction.json
```

This searches abstract local tuples around `m` with:

- `g<0`,
- both mode inequalities `I1>=0`, `I2>=0`,
- local LC for `P` and `Q`,
- local dominance `q_{k+1} <= p_k`,
- optional `Q` shape constraints.

Found explicit witnesses in every regime. Strongest regime (`Q` shape + positive tail):

- `P`-local: `(p_{m-2},p_{m-1},p_m,p_{m+1}) = (5,4,3,1)` (so `g=-1`)
- `Q`-local: `(q_{m-2},q_{m-1},q_m,q_{m+1}) = (1,2,4,1)`
- `I1 = 0`, `I2 = 5`.

So these inequalities + local LC/dominance are not sufficient by themselves to deduce `g>=0`; extra tree-structural information is required.

## 6. Takeaway for the direct attack

The direct coefficient route gives a clean reduction:

- proving `mode(P) >= m-1` is equivalent to forbidding `Q`-compensation patterns that can offset negative `P` slope in `I1`.

Empirically (full `n<=23`), the target holds with margin (`g>=1` always), while `I1` can be near-tight and often uses substantial `Q` compensation. Therefore a symbolic closure likely needs one more structure lemma linking:

- two-step rise `q_m-q_{m-2}` (or `P'` local behavior),
- to `P` local slope near `m-1`,
- beyond generic LC + coefficientwise domination.

