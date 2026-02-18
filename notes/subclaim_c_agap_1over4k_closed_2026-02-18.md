# Sub-claim C: closed-form proof of `A_gap >= 1/(4k)` (2026-02-18)

Goal (remaining algebraic gap in Sub-claim C):

- For `k >= 6`, `k ≡ 0,2 (mod 3)`, show
  `A_gap := A_0 - A_1 >= 1/(4k)`.

Here:

- `A_0 = 2k*lambda_0/(1+2lambda_0) - (m_0-1)`,
- `A_1 = 2k*lambda_1/(1+2lambda_1) + lambda_1/(1+lambda_1) - (m_1-1)`,
- `m_0=(2k+1)//3`, `m_1=(2k+3)//3`,
- `lambda_0, lambda_1` are the exact tie-fugacities used in `prove_subclaim_c_algebra.py`.

---

## 1) Exact closed forms for `A_gap - 1/(4k)`

Set `p = 4^t`.

For `k=3t` (`t>=2`):

`A_gap - 1/(4k) = N_0(t,p) / D_0(t,p)`,

with

- `D_0 = 12*t*(3pt^2+7pt+2p+10t^2)*(8pt^2+9pt+2p+6t^2+5t+1)*(12pt^2+13pt+3p+10t^2+7t+1)`,
- `N_0 = p^3*C3_0(t) + p^2*C2_0(t) + p*C1_0(t) - C0_0(t)`,
- `C1_0(t) = (2t+1)R_0(t)`,
- `R_0(t) = 1296t^7-2928t^6-6022t^5-2387t^4+15t^3+91t^2-3t-2`,
- `C0_0(t) = 2t^2(2t+1)^2(348t^3+147t^2+52t+5)`.

For `k=3t+2` (`t>=1`):

`A_gap - 1/(4k) = N_2(t,p) / D_2(t,p)`,

with

- `D_2 = 4*(3t+2)*(8pt+13p+3t+3)*(12pt+18p+5t+4)*(6pt^2+24pt+18p+10t^2+11t+3)`,
- `N_2 = p^3*C3_2(t) + p^2*C2_2(t) + p*C1_2(t) - C0_2(t)`,
- `C1_2(t) = 3R_2(t)`,
- `R_2(t) = 432t^6+608t^5-2542t^4-6900t^3-6301t^2-2455t-338`,
- `C0_2(t) = (2t+1)(348t^4+1015t^3+1180t^2+637t+132)`.

All denominator factors are positive for `t>=1`, `p>0`, so sign is controlled by `N_0, N_2`.

---

## 2) Tail positivity for `k=3t` (`t>=4`)

Because `N_0 = p^3*C3_0 + p^2*C2_0 + p*C1_0 - C0_0`, it is enough to prove:

`p*C1_0 - C0_0 > 0`.

First bound `C1_0` from below.

`R_0 - 150t^7 = t^4*Q(t) + (15t^3+91t^2-3t-2)`,
where `Q(t)=1146t^3-2928t^2-6022t-2387`.

- `Q''(t)=6876t-5856 > 0` for `t>=1`, so `Q'` is increasing.
- `Q'(4)=25562>0`, so `Q' > 0` for `t>=4`.
- `Q(4)=9621>0`, hence `Q(t)>0` for `t>=4`.

Therefore `R_0 >= 150t^7` for `t>=4`, and

`C1_0 = (2t+1)R_0 >= 300t^8`.

Now bound `C0_0` above:

- `(2t+1)^2 <= (3t)^2`,
- `348t^3+147t^2+52t+5 <= 552t^3` for `t>=1`,

so `C0_0 <= 9936t^7`.

Hence:

`N_0 >= p*C1_0 - C0_0 >= 4^t*300t^8 - 9936t^7 = t^7(300t4^t-9936) > 0`

for `t>=4` (already positive at `t=4`).

So `A_gap - 1/(4k) > 0` for `k=3t`, `t>=4` (`k>=12`).

---

## 3) Tail positivity for `k=3t+2` (`t>=4`)

Similarly, enough to prove:

`p*C1_2 - C0_2 > 0`.

For `t>=4`, absorb lower-order negatives in `R_2`:

- `-6900t^3 >= -1725t^4`,
- `-6301t^2 >= -394t^4`,
- `-2455t >= -39t^4`,
- `-338 >= -2t^4`.

Thus

`R_2 >= 432t^6+608t^5-4702t^4 = t^4(432t^2+608t-4702)`.

Since `332t^2+608t-4702 > 0` for `t>=4`, we get

`R_2 >= 100t^6`, hence `C1_2 = 3R_2 >= 300t^6`.

Also

`C0_2 <= (3t)*(3312t^4) = 9936t^5`.

So:

`N_2 >= p*C1_2 - C0_2 >= 4^t*300t^6 - 9936t^5 = t^5(300t4^t-9936) > 0`

for `t>=4`.

Therefore `A_gap - 1/(4k) > 0` for `k=3t+2`, `t>=4` (`k>=14`).

---

## 4) Finite bases

The only `k>=6`, `k ≡ 0,2 (mod 3)` not in the two tails above are:

- `k in {6, 8, 9, 11}`.

Exact values:

- `k=6`: `A_gap=0.04210996899181141`, `1/(4k)=0.041666666666666664`.
- `k=8`: `A_gap=0.046869125179519106`, `1/(4k)=0.03125`.
- `k=9`: `A_gap=0.0400014840066271`, `1/(4k)=0.027777777777777776`.
- `k=11`: `A_gap=0.03657205192274994`, `1/(4k)=0.022727272727272728`.

All satisfy `A_gap >= 1/(4k)`.

Hence the bound holds for all `k>=6`, `k ≡ 0,2 (mod 3)`.

---

## 5) Reproducibility

New script:

- `prove_subclaim_c_agap_bound.py`

Run:

`python3 prove_subclaim_c_agap_bound.py --t-max-id 60 --t-max-lane 200 --k-max-scan 600`

Output in this run:

- closed-form `num/den` identities: PASS,
- finite bases `{6,8,9,11}`: PASS,
- large-`t` lane checks: PASS,
- direct exact scan (`k<=600`): `0` failures.
