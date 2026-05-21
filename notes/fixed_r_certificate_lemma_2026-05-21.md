# Fixed-`r` Route-2 Certificate Lemma

## Purpose

This note extracts the current fixed-`r` computation into a theorem schema.
It does not prove the arbitrary fixed-`r` theorem by itself.  It proves the
conditional implication:

```text
if the finite, mode, hub-off, and hub-on certificates pass for a lane
S(2^a,r), then Route-2 holds for that lane for every a.
```

This is the bridge between the exact certificate scripts and a manuscript
proof.

## Setup

Let `P_j(x)` be the independence polynomial of the path on `j` vertices.  For
the spider

```text
T_{a,r} = S(2^a,r),
```

write

```text
I(T_{a,r};x) = F_a(x) + G_a(x),
F_a(x) = P_r(x)(1+2x)^a,
G_a(x) = xP_{r-1}(x)(1+x)^a.
```

After removing one length-2 arm, the bridge tree is

```text
B_{a,r} = S(2^(a-1),r),
```

with polynomial

```text
I(B_{a,r};x) = F^-_a(x) + G^-_a(x),
F^-_a(x) = P_r(x)(1+2x)^(a-1),
G^-_a(x) = xP_{r-1}(x)(1+x)^(a-1).
```

For a polynomial `H(x)=sum h_k x^k`, define

```text
mu_H(lambda) = lambda H'(lambda) / H(lambda).
```

The Route-2 target is:

```text
mu_{B_{a,r}}(lambda) >= m - 3/2,
lambda = i_{m-1}(T_{a,r}) / i_m(T_{a,r}),
```

where `m` is the leftmost mode of `I(T_{a,r})`.

## Certificate Data

Fix `r` and a threshold `A`.  For each residue class `q in {0,1,2}`, the
certificate supplies an integer shift `D_q` and sets

```text
a = 3t + q,
m = (2a + D_q)/3.
```

The certificate also supplies positive denominators:

```text
M = mode-margin denominator,
R = hub-off reserve denominator.
```

Current runs usually use `M=R=1000`, except the early `r=8` certificate where
stronger denominators were also checked.

## Lemma

For a fixed lane `r`, suppose the following four certificates hold.

### C0. Finite Certificate

For every `1 <= a < A`, the exact finite check verifies

```text
mu_{B_{a,r}}(lambda) >= m - 3/2,
```

where `m` is the actual leftmost mode of `I(T_{a,r})` and
`lambda=i_{m-1}(T_{a,r})/i_m(T_{a,r})`.

### C1. Mode Certificate

For every `a >= A`, with `a == q mod 3`, the integer `m=(2a+D_q)/3` is the
leftmost mode of `I(T_{a,r})`.

A sufficient way to prove this is:

```text
F_m - F_k >= F_m/(M a)       for every k != m,
2 max_k G_k < F_m/(M a).
```

Indeed,

```text
(F_m+G_m) - (F_k+G_k)
  >= (F_m-F_k) - |G_m-G_k|
  >= F_m/(M a) - 2 max_j G_j
  > 0.
```

In the current scripts, the global `F`-mode condition is reduced to adjacent
margin checks plus the real-rooted/log-concave structure of
`P_r(x)(1+2x)^a`.  This reduction should be stated explicitly in the paper.

### C2. Hub-Off Reserve Certificate

Let

```text
lambda_0 = F_{m-1}/F_m.
```

For every `a >= A`, the certificate proves:

```text
3/4 <= lambda_0 <= 2
```

and

```text
mu_{F^-_a}(lambda_0) >= m - 4/3 + 1/(R a).
```

### C3. Hub-On Perturbation Certificate

For every `a >= A`,

```text
|mu_{F^-_a+G^-_a}(lambda) - mu_{F^-_a}(lambda_0)| <= 1/(R a).
```

The current implementation proves this by bounding two errors:

```text
1. fugacity shift: lambda versus lambda_0;
2. mixture shift: replacing F^-_a by F^-_a+G^-_a.
```

The proof uses the compact interval for `lambda_0`, the mode-domination bound
on `G/F`, and a conservative variance estimate.

### Conclusion

Under C0-C3,

```text
mu_{B_{a,r}}(lambda) >= m - 3/2
```

for every `a >= 1`.

## Proof

For `a < A`, this is exactly C0.

Now take `a >= A`.  By C1, the certificate's `m` is the actual leftmost mode of
`I(T_{a,r})`, so the Route-2 fugacity is

```text
lambda = i_{m-1}(T_{a,r})/i_m(T_{a,r}).
```

By C2 and C3,

```text
mu_{B_{a,r}}(lambda)
  = mu_{F^-_a+G^-_a}(lambda)
  >= mu_{F^-_a}(lambda_0) - 1/(R a)
  >= m - 4/3.
```

Since

```text
m - 4/3 > m - 3/2,
```

the Route-2 target follows.

## Script Mapping

### C0 finite certificate

Script:

```bash
python3 gpt_attack/fixed_r_finite_route2_check.py
```

Exact sign checked:

```text
sum_k (2k - 2m + 3) b_k L^k M^(d-k) >= 0
```

where `lambda=L/M`, `I(B)=sum b_k x^k`, and `d=deg B`.

### C1 mode certificate

Scripts:

```bash
python3 gpt_attack/fixed_r_huboff_certificate.py
python3 gpt_attack/fixed_r_hubon_mode_certificate.py
```

The hub-off script checks the `F` boundary margins and the hub-on mode script
checks:

```text
2 max_k G_k / F_m < 1/(M a)
```

at the threshold, plus monotonicity as `a -> a+3` inside each residue class.

Paper obligation:

```text
state why the checked F boundary margins imply the needed global F margin.
```

The likely route is real-rootedness/log-concavity of `P_r(x)(1+2x)^a`.

### C2 hub-off reserve certificate

Scripts:

```bash
python3 gpt_attack/fixed_r_huboff_certificate.py
python3 gpt_attack/fixed_r_huboff_cpp_certificate.py
```

The Python path works for smaller `r`; the C++/GMP helper handles dense
high-degree positivity for larger `r`.

The core exact recurrence is:

```text
E_n = Z^deg(P_n) P_n(L/Z),
N_n = Z^deg(P_n) (L/Z)P'_n(L/Z),
```

where `lambda_0=L/Z`.

The shifted numerator and denominator polynomials are checked to have positive
coefficients after the threshold shift in each residue class.

### C3 hub-on perturbation certificate

Script:

```bash
python3 gpt_attack/fixed_r_hubon_route2_perturbation.py
```

It proves that the sum of the fugacity-shift error and mixture error is below
`1/(R a)` at the threshold and remains decreasing in each residue class.

## Current Instantiations

Complete lanes under this lemma:

```text
sampled set: {4,8,12,16,20,24,32,40,60,80}, A=200
r=100, A=200
r=120, A=250
r=160, A=300
r=200, A=350
r=240, A=400
r=280, A=450
```

Partial next lane:

```text
r=320, A=500:
  C0 finite certificate passes for a <= 499;
  C3 hub-on perturbation and hub-on mode certificates pass for a >= 500;
  C2 hub-off reserve has not yet been run.
```

## Load-Bearing Assumptions and Failure Conditions

### A1. Boundary margins imply global F margins

Failure condition:

```text
F_m beats its neighbors but some farther F_k is within F_m/(M a) or larger.
```

Why this is unlikely:

```text
F=P_r(1+2x)^a has real nonpositive roots, hence strong coefficient
log-concavity/unimodality.  But the quantitative global margin still needs a
clean statement.
```

Status:

```text
Drafted in notes/fixed_r_global_f_mode_margin_sublemma_2026-05-21.md.
```

### A2. C3 perturbation is stated for the actual lambda

Failure condition:

```text
mode localization identifies m, but the bound on |lambda-lambda_0| is too weak
or depends on a hidden stronger margin than C1 supplies.
```

Status:

```text
Drafted in notes/fixed_r_fugacity_shift_sublemma_2026-05-21.md.
```

### A3. Positive shifted coefficients are accepted as exact proof

Failure condition:

```text
the certificate proves positivity only for the sampled integer t values, not
for all t >= threshold.
```

Why current checks avoid this:

```text
they shift t by the threshold and check every coefficient is positive, so the
polynomial is positive for every nonnegative integer t.
```

Status:

```text
Drafted in notes/fixed_r_shifted_coefficient_positivity_2026-05-21.md.
```

## Next Proof-Writing Step

The fugacity-shift sublemma has been drafted in:

```text
notes/fixed_r_fugacity_shift_sublemma_2026-05-21.md
```

It proves the two bounds used in `fixed_r_hubon_route2_perturbation.py`:

```text
|lambda-lambda_0| <= 4T,
|mu_{F^-}(lambda)-mu_{F^-}(lambda_0)| <= 2N^2T.
```

The next remaining proof-writing target is the global `F`-mode margin:

```text
show that the checked adjacent margins for F=P_r(1+2x)^a imply the global
margin F_m-F_k >= F_m/(M'a) needed for mode domination.
```

This has now been drafted in:

```text
notes/fixed_r_global_f_mode_margin_sublemma_2026-05-21.md
```

It uses real-rootedness of `P_r(x)(1+2x)^a` and unimodality to reduce the
global margin to the checked adjacent margins.

The shifted-coefficient positivity principle and the hub-off reserve recurrence
certificate have also been drafted in:

```text
notes/fixed_r_shifted_coefficient_positivity_2026-05-21.md
notes/fixed_r_huboff_reserve_recurrence_certificate_2026-05-21.md
```

At this point, the local certificate-composition proof is written as notes.
The remaining gap is global: producing the certificate data for arbitrary fixed
`r`, or replacing the per-lane hub-off reserve certificate by an analytic
lower bound.
