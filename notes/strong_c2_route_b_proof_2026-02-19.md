# STRONG C2 Route B: Algebraic Attempt on `R = cross + mismatch`

Date: 2026-02-19

## 1. Setup

For `T` with `d_leaf <= 1`, choose canonical leaf `l` with degree-2 support `s`, and let `u` be the other neighbor of `s`.

Define:
- `B = T - {l, s}`
- `P = dp_B[u][0]`
- `Q = dp_B[u][1] = x G`, where `G = prod_c g_c` and `g_c = dp0[c]` over children `c` of `u` in `B`
- `P = prod_c (g_c + h_c)` with `h_c = dp1[c]`

At `m = mode(I(T))`, set `r = m-2` and
- `p0 = p_r`, `p1 = p_{r+1}`, `pm = p_{r+2}`
- `q0 = q_r = G_{r-1}`, `q1 = q_{r+1} = G_r`, `qm = q_{r+2} = G_{r+1}`

Route B residual:

`R = 2*p1*q1 - pm*q0 - p0*qm + p0*q1 - p1*q0`.

Equivalent split:
- `D1 = p1*q1 - pm*q0`
- `D2 = p1*q1 - p0*qm`
- `D3 = p0*q1 - p1*q0 = mismatch`
- `R = D1 + D2 + D3`
- `cross = D1 + D2`

---

## 2. Exact product-linearity identity

For fixed `G` and `r`, define the linear functional in `A`:

`L(A) = 2*a_{r+1}*G_r - a_{r+2}*G_{r-1} - a_r*G_{r+1} + a_r*G_r - a_{r+1}*G_{r-1}`.

Then `R = L(P)` exactly.

Since
`P = sum_{S subseteq children(u)} H_S`, with
`H_S = (prod_{i in S} h_i) (prod_{i notin S} g_i)`,
we get:

`R = sum_S L(H_S)`.

This is the clean algebraic reduction from the product structure.

---

## 3. What worked and what failed

### 3.1 Term-level behavior (n<=23, 931,595 trees)

Using `verify_strong_c2_route_b_pair_bounds_2026_02_19.py`:

- `R < 0`: `0` failures, `min(R)=4`.
- `cross < 0`: `0` failures, `min(cross)=3`.
- `D2 + D3 < 0`: `0` failures, `min(D2+D3)=2`.
- `cross - |mismatch| < 0`: `0` failures, `min=2`.

But:
- `D1 < 0`: `1` witness.
- `D2 < 0`: many (`14,283`).
- `D3 < 0`: sparse (`129`).
- `D1 + D3 < 0`: `7` witnesses.

So the earlier heuristic “Term A (`D1`) is always nonnegative” is false on this frontier.

Explicit `D1<0` witness:
- `n=22`, `m=7`, `r=5`
- `g6 = U?????????O?O?G?A??O?@??A??A??@???O?~}??`
- `(p0,p1,pm)=(6048,9408,9984)`
- `(q0,q1,qm)=(126,126,84)`
- `D1=-72576`, `D2=677376`, `D3=-423360`, `R=181440`

### 3.2 Subset expansion for full `R`

Using `verify_strong_c2_route_b_componentwise_2026_02_19.py`:

- `R` linearity identities: `0` failures.
- Negative subset contributions `L(H_S)` do occur (`328` trees in explicit subset checks, `min=-211320`).
- Negative `|S|`-layer contributions also occur (`331` trees, layer min `-396128`).

So `R = sum_S L(H_S)` is exact but **not** termwise-positive.

---

## 4. Strongest new candidate inequality

Define

`M := D2 + D3 = p1*(q1-q0) + p0*(q1-qm)`.

Empirically (full n<=23):
- `M >= 2` always (`M_neg = 0`, `M_min = 2`).

This was confirmed independently in:
- `verify_strong_c2_route_b_pair_bounds_2026_02_19.py`
- `verify_strong_c2_route_b_subset_M_2026_02_19.py`

### 4.1 Algebraic form via `P = G + E`

Write `E = P - G` (`E_k >= 0`). Then:

`M = (G_r^2 - G_{r-1}G_{r+1}) + E_{r+1}(G_r-G_{r-1}) + E_r(G_r-G_{r+1})`.

So `M` is LC-surplus of `G` plus an `E`-correction.

This is structurally useful, but I do not yet have a general sign proof for the correction term.

### 4.2 Subset test for `M`

For `M_A := a_{r+1}(G_r-G_{r-1}) + a_r(G_r-G_{r+1})`, we have
`M = sum_S M_{H_S}`.

Result (explicit subset checks, n<=23):
- almost all subset contributions are nonnegative,
- but there is one negative subset contribution witness:
  - `n=23`, `m=7`, `r=5`, `deg_u(B)=1`, `mask=1`
  - `g6 = V?????????????????E??K??K??E??@_??K???oDTT??`
  - contribution `-27984`
  - total `M=7860000` (still strongly positive).

So even for `M`, full termwise subset positivity is false.

---

## 5. Cauchy-Schwarz / rearrangement status

### 5.1 AM-GM/Cauchy route on `cross`

The standard LC bounds give
`p1^2 >= p0*pm`, `q1^2 >= q0*qm`, hence
`p1*q1 >= sqrt(p0*pm*q0*qm)`.

This does not imply
`2*p1*q1 >= pm*q0 + p0*qm`.
So it does not directly prove `cross >= 0` or `R >= 0`.

### 5.2 Rearrangement form of `M`

`M >= 0` is equivalent to

`q1 >= (p1*q0 + p0*qm)/(p0+p1)`.

This looks like a weighted midpoint condition for `q`, but LC alone (`q1^2 >= q0*qm`) does not force this weighted average bound for arbitrary `(p0,p1)`.

No clean rearrangement closure proof found yet.

---

## 6. Current best algebraic reduction

We now have:

`R = D1 + M`, where
- `D1 = p1*q1 - pm*q0`,
- `M = p1*(q1-q0) + p0*(q1-qm)`.

Empirically:
- `M >= 2` always,
- `D1` can be negative (rarely),
- but `R >= 4` always.

So the unresolved symbolic gap is: a general lower bound coupling `D1` and `M` (or directly `cross` vs `|mismatch|`) from the product structure.

---

## 7. New scripts and artifacts (new files only)

Scripts added:
- `verify_strong_c2_route_b_componentwise_2026_02_19.py`
- `verify_strong_c2_route_b_pair_bounds_2026_02_19.py`
- `verify_strong_c2_route_b_subset_M_2026_02_19.py`
- `find_negative_subset_M_witness_2026_02_19.py`

Main outputs:
- `results/verify_strong_c2_route_b_componentwise_2026_02_19_n23.json`
- `results/verify_strong_c2_route_b_pair_bounds_2026_02_19_n23.json`
- `results/verify_strong_c2_route_b_subset_M_2026_02_19_n23.json`

Commands used:

```bash
python3 verify_strong_c2_route_b_componentwise_2026_02_19.py \
  --min-n 4 --max-n 23 \
  --out results/verify_strong_c2_route_b_componentwise_2026_02_19_n23.json

python3 verify_strong_c2_route_b_pair_bounds_2026_02_19.py \
  --min-n 4 --max-n 23 \
  --out results/verify_strong_c2_route_b_pair_bounds_2026_02_19_n23.json

python3 verify_strong_c2_route_b_subset_M_2026_02_19.py \
  --min-n 4 --max-n 23 \
  --out results/verify_strong_c2_route_b_subset_M_2026_02_19_n23.json

python3 find_negative_subset_M_witness_2026_02_19.py \
  --min-n 23 --max-n 23
```

---

## 8. Bottom line

I do **not** have a complete symbolic proof of `R >= 0` yet.

What is new and solid:
- exact product-linearity decomposition for `R`;
- correction of a false intermediate claim (`D1` not always nonnegative);
- strong universal inequality candidate `M = D2 + D3 >= 0` (empirical min `2`), with an explicit algebraic expression in `G` and `E=P-G`.

The most promising next symbolic target is proving `M >= 0` from product structure and then deriving a coupled lower bound for `D1` in the rare `D1<0` regime.
