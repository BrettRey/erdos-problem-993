# Sub-claim A, Lane #2: envelope closure attempt (2026-02-18)

Goal in lane #2:

- prove `lambda(k,j+2) >= lambda(k,j)` for mixed spiders `S(2^k,1^j)`.

Using prior notation (`m = mode(k,j)`):

- `A=f_m, B=f_{m-1}, C=f_{m-2}, D=f_{m+1}`,
- `G=g_m, H=g_{m-1}`,
- `f_t=[x^t](1+2x)^k(1+x)^j`,
- `g_t=[x^t]x(1+x)^k`.

Equivalent inequality:

- `F := (A+2B+C)(A+G) - (D+2A+B)(B+H) >= 0`.

Write:

- `F = T1 + T2`,
- `T1 = (A^2-BD) + (AC-B^2)`,
- `T2 = G(A+2B+C) - H(2A+B+D)`.

## 1) Envelope idea for the hard sign (`T2<0`)

Let `s := G/A`. Then

- `F/A^2 = T1/A^2 + s * (T2/(A G))` when `G>0`.

If `T2<0`, this is decreasing in `s`, so any explicit upper bound `s <= s_hat`
gives

- `F/A^2 >= T1/A^2 + s_hat * (T2/(A G))`.

## 2) Explicit combinatorial upper bound on `s=G/A`

From `A=[x^m](1+2x)^k(1+x)^j`:

1. If `m <= k`, use the `a=m` and `a=m-1` terms:
   - `A >= 2^m C(k,m) + j 2^(m-1) C(k,m-1)`,
   - hence
     `s = C(k,m-1)/A <= m / (2^m (k-m+1) + j 2^(m-1) m)`.

2. If `m = k+1`, use `a=k` and `a=k-1`:
   - `A >= 2^k j + 2^(k-1) k C(j,2)`,
   - hence
     `s <= 1 / (2^k j + 2^(k-1) k C(j,2))`.

(`m>k+1` implies `G=0`, so this branch is irrelevant.)

## 3) New verifier

Script:

- `prove_subclaim_A_lane2_envelope.py`

What it checks:

1. direct `F>=0`,
2. mode-step `m(k,j+2)=m(k,j)+1`,
3. `mode(I_{k,j}) = mode((1+2x)^k(1+x)^j)`,
4. correctness of the explicit `s` upper bound above,
5. envelope inequality in `T2<0` regime:
   - `LB := T1/A^2 + s_hat * T2/(AG) >= 0`.

## 4) Runs

### Run A

`python3 prove_subclaim_A_lane2_envelope.py --k-max 3000 --j-max 120 --out results/subclaimA_lane2_envelope_k3000_j120.json`

Results:

- checked pairs: `362,395`,
- direct `F>=0` failures: `0`,
- envelope failures (`LB<0`): `0`,
- `s`-upper-bound failures: `0`,
- mode-step failures: `0`,
- mode mismatch (`I` vs product mode): `0`,
- `T2<0` pairs: `352,648`,
- `T2>=0` pairs: `0`.

Tight witnesses:

- `min F/A^2 = 1.3129329744128296e-06` at `(k,j,m)=(3000,119,2060)`,
- `min envelope LB = 1.3129329744128296e-06` at the same witness,
- `max (s/s_hat) = 1` at `(k,j,m)=(6,0,4)` (bound is tight there).

### Run B

`python3 prove_subclaim_A_lane2_envelope.py --k-max 1200 --j-max 160 --out results/subclaimA_lane2_envelope_k1200_j160.json`

Results:

- checked pairs: `192,395`,
- same zero-failure pattern as Run A,
- `min F/A^2 = min envelope LB = 6.1487527760063325e-06` at `(1200,159,880)`.

### Run C (patched script with explicit ratio output)

`python3 prove_subclaim_A_lane2_envelope.py --k-max 3000 --j-max 120 --out results/subclaimA_lane2_envelope_k3000_j120_v2.json`

Results:

- same zero-failure pattern as Run A,
- explicit ratio minimum:
  - `min envelope ratio = 1.04020450460594072` at `(k,j,m)=(6,3,6)`.

## 5) Status

This gives a stronger algebraic lane than raw brute-force `F>=0` checking:

- the difficult `T2<0` regime is controlled by an explicit closed `s` envelope;
- the envelope certificate passed all large scans above with zero failures.

Remaining closure task is to convert this into a full symbolic proof (without
range scans) that the envelope lower bound is nonnegative for all `k>=6, j>=0`.

## 6) Ratio form of the envelope condition

Equivalent ratio (in the `T2<0` branch):

- `R := (T1 * A_lb) / (|T2| * A)`,
  where `A_lb` is the explicit two-term lower bound on `A`.
- Envelope condition is `R >= 1`.

Full-range computation (`k=6..3000`, `j=0..120`) gives:

- minimum `R = 1.0402045046059407` at `(k,j,m)=(6,3,6)`,
- checked pairs in `T2<0` branch: `352,648`.

So the numeric margin in this ratio form is about `+4.02%` over the tested range.
The remaining task is symbolic closure of `R>=1`.
