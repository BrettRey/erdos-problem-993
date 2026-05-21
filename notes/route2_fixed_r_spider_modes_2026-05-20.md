# Fixed-`r` Spider Mode Probe (2026-05-20)

## Purpose

After the pure-spider proof, the next natural family is:

```text
S(2^a,r)
```

with `r` fixed and `a -> infinity`, removing a length-2 arm for the Route-2
bridge.

The independence polynomial has the exact form:

```text
I_{a,r}(x) = P_r(x)(1+2x)^a + x P_{r-1}(x)(1+x)^a,
```

where `P_r` is the independence polynomial of the path on `r` vertices.  The
second term is again exponentially small near the large-`a` saddle.  Thus the
mode should be governed by the fixed convolution `P_r(x)(1+2x)^a`, with an
eventually constant residue-class shift:

```text
3 m_{a,r} - 2a = D_{r,q},     q = a mod 3.
```

## Probe

Script:

```bash
python3 gpt_attack/fixed_r_spider_mode_probe.py \
  --r-max 80 --a-start 200 --a-end 260 \
  --out results/fixed_r_spider_mode_shifts_r80_a200_260.json
```

Output:

```text
checked r=2..80, a=200..260
unstable records in window: 0

selected shifts D_q = 3m - 2a
r= 2: q0=  3 q1=  1 q2=  2
r= 3: q0=  3 q1=  4 q2=  2
r= 4: q0=  3 q1=  4 q2=  5
r= 5: q0=  6 q1=  4 q2=  5
r= 6: q0=  6 q1=  7 q2=  5
r= 7: q0=  6 q1=  7 q2=  8
r= 8: q0=  9 q1=  7 q2=  8
r= 9: q0=  9 q1=  7 q2=  8
r=10: q0=  9 q1= 10 q2=  8
r=12: q0= 12 q1= 10 q2= 11
r=16: q0= 15 q1= 13 q2= 14
r=20: q0= 18 q1= 19 q2= 17
r=40: q0= 33 q1= 34 q2= 35
r=80: q0= 66 q1= 67 q2= 68
```

A wider window starting at `a=50` was also run:

```bash
python3 gpt_attack/fixed_r_spider_mode_probe.py \
  --r-max 80 --a-start 50 --a-end 260
```

Only three `r` values had non-stabilized shifts over this wider window:

```text
r=26: q0 in {21,24}
r=55: q0 in {45,48}
r=61: q2 in {50,53}
```

By `a >= 200`, all tested `r <= 80` are stabilized.

## Interpretation

This is not a theorem yet, but it makes the next proof target concrete:

```text
For each fixed r, prove eventual mode localization
m_{a,r} = (2a + D_{r,a mod 3})/3
```

with an effectively computable finite threshold.  Once that is known, the
large-`a` Route-2 asymptotic should follow the same pattern as the pure
spider:

```text
Delta_{a,r} = Delta^0_{a,r} + O_r(a beta_r^a),
```

where `Delta^0_{a,r}` is the hub-off/binomial-convolution slack and the
hub-on branch contributes the exponentially small error.

The finite endpoint dip at `S(2^12,8)` is therefore probably a small-`a`
phenomenon, not the start of a worse asymptotic family.

## `r=8` Symbolic Certificate Pieces

The first serious fixed-`r` target should be `r=8`, since the endpoint scan's
finite worst witness is `S(2^12,8)`.

Script:

```bash
python3 gpt_attack/fixed_r8_symbolic_certificate.py
```

For `S(2^a,8)`, the polynomial is:

```text
I_{a,8}(x) = P_8(x)(1+2x)^a + xP_7(x)(1+x)^a
```

with:

```text
P_8(x) = 1 + 8x + 21x^2 + 20x^3 + 5x^4,
P_7(x) = 1 + 7x + 15x^2 + 10x^3 + x^4.
```

The stabilized shifts are:

```text
a = 3t:     m = 2t + 3,     3m - 2a = 9,
a = 3t+1:   m = 2t + 3,     3m - 2a = 7,
a = 3t+2:   m = 2t + 4,     3m - 2a = 8.
```

These shifts are asymptotic/stabilized shifts, not small-`a` identities.  The
tiny cases `a=6` and `a=9` have the actual full mode one step left of this
formula.

The script derives exact formulas for:

- the full-polynomial boundary differences `c_m-c_{m-1}` and `c_{m+1}-c_m`;
- the hub-off fugacity `lambda_0`;
- the hub-off reserve `Delta^0_{a,8} - 1/6`.

For all three residue classes, the hub-off reserve is a rational function in
`t` whose numerator and denominator have strictly positive coefficients.  This
proves:

```text
Delta^0_{a,8} > 1/6
```

for the stabilized hub-off lane.

Exact Route-2 check for the complete `S(2^a,8)` lane:

```bash
python3 - <<'PY'
from gpt_attack.route2_spider_lane_scan import path_polys, route2_for_removed_arm
paths = path_polys(8)
mins = None
fails = []
for a in range(1, 81):
    rec = route2_for_removed_arm({2: a, 8: 1}, 2, paths)
    if rec["route2_slack"] <= 0:
        fails.append(a)
    if mins is None or rec["route2_slack"] < mins[1]:
        mins = (a, rec["route2_slack"], rec["m"])
print(len(fails), mins[0], float(mins[1]), mins[2])
PY
```

Output:

```text
checked a=1..80 exactly
failures: 0
minimum slack: a=78, 0.1710112823235669, m=55
```

## Remaining Proof Obligations for `r=8`

The `r=8` lane is not fully proved yet.  What remains:

1. Convert the boundary-difference formulas into a full mode-localization
   proof for the complete polynomial.  Boundary signs are in hand; the missing
   step is a clean global-unimodality/log-concavity or domination argument.
2. Bound the exponentially small perturbation from
   `xP_7(x)(1+x)^a`, showing the full fugacity and full mean stay within the
   positive hub-off reserve.
3. Use exact finite checking only up to the resulting threshold.

## `r=8` Domination Probe

The current best mode-proof mechanism is domination:

- `F_a=P_8(x)(1+2x)^a` is the hub-off term.
- `G_a=xP_7(x)(1+x)^a` is the hub-on term.
- `G_a` peaks well before the target mode.
- On the interval from the `G_a` peak to the target mode, `G_a/F_a` is
  exponentially small, so the increasing margin of `F_a` should dominate any
  negative difference from `G_a`.

Script:

```bash
python3 gpt_attack/fixed_r8_domination_probe.py \
  --a-max 250 \
  --out results/fixed_r8_domination_probe_a250.json
```

Output summary:

```text
checked a=1..250
diagnostic failures: 9
min Route-2 slack: a=249 slack=0.16810129314
```

All diagnostic failures are small (`a <= 12`).  From `a=30` onward, the
domination condition is very comfortable:

```text
a=  30 m=  23 gmode=  18 tailE=6.355e-06 deltaL=2.781e-03 margin=2.768e-03
a=  81 m=  57 gmode=  44 tailE=8.422e-14 deltaL=1.403e-03 margin=1.403e-03
a= 200 m= 136 gmode= 103 tailE=1.337e-31 deltaL=7.921e-03 margin=7.921e-03
```

This suggests the right proof split for `r=8`:

```text
exact finite check: a <= 29,
symbolic domination: a >= 30.
```

After trying to make the constants proof-ready, the cleaner conservative split
is:

```text
exact finite Route-2 check: a <= 200,
asymptotic proof:          a >= 200.
```

Script:

```bash
python3 gpt_attack/fixed_r8_asymptotic_bounds.py
```

Output:

```text
finite exact Route-2 range: a=1..200
finite failures: 0
minimum finite slack: a=198 slack=0.168459434752 m=135

mode localization crude bound at a=200:
  918(2/3)^66 ~= 2.19189418012e-09
  2*918(2/3)^66 ~= 4.38378836024e-09
  comparison reserve 1/(10a) = 0.0005

perturbation crude bound at a=200:
  10000 a^8 / 2^(a/2) ~= 2.01948391737e-08
  10000 a^3 (3/4)^(a-1) ~= 1.09718889151e-14
  total perturbation bound ~= 2.01948501455e-08
  comparison reserve 1/(4a) = 0.00125
```

All exponential comparison terms decrease after `a=200`; the successive-ratio
bounds at `a=200` are below `1` in the certificate script.

The formal ingredients for the `a >= 200` proof are now:

1. Hub-off relative margin:

   ```text
   (F_m-F_{m-1})/F_m >= 1/(10a).
   ```

   The residue-class formulas show this; the only slightly delicate class is
   `a=3t`, where the numerator of
   `delta_F - 1/(10a)` is

   ```text
   40t^4 - 572t^3 - 947t^2 - 441t - 60,
   ```

   positive for `t >= 20`.  The companion script
   `gpt_attack/fixed_r8_symbolic_inequality_certificate.py` checks this and
   the right margin by shifting `t` to the threshold and verifying positive
   coefficients.
2. Hub-off reserve:

   ```text
   Delta^0_{a,8} - 1/6 >= 1/(4a).
   ```

   The same symbolic certificate verifies this in all three residue classes
   for `a >= 200`, along with `3/4 <= lambda_0 <= 2`.
3. Mode localization:

   Let `F_k=[x^k]P_8(x)(1+2x)^a` and
   `G_k=[x^k]xP_7(x)(1+x)^a`.  Since `P_8` is real-rooted with negative roots,
   `F` is log-concave and unimodal.  The margins in item 1 imply

   ```text
   F_m - F_k >= F_m/(10a)       for every k != m.
   ```

   Meanwhile

   ```text
   G_k <= 34 * 2^a,
   F_m >= 2^m binom(a,m) >= 2^(2a/3) 3^(a/3-3),
   ```

   so

   ```text
   G_k/F_m <= 918(2/3)^(a/3).
   ```

   Comparing full coefficients requires the two-sided bound
   `|G_m-G_k| <= G_m+G_k <= 2G_max`, so the certificate checks

   ```text
   2*918(2/3)^(a/3) < 1/(10a)
   ```

   for `a >= 200`.  Hence the full coefficient `C_k=F_k+G_k` has its unique
   mode at the stabilized `m`.
4. Full fugacity/mean perturbation:

   At the mode indices, the ratios `G_k/F_k` are bounded by
   `34a^5/2^(a/2)`.  This gives a fugacity shift contribution bounded by
   `10000a^8/2^(a/2)`.  The residual hub-on mixture in `B=S(2^(a-1),8)` is
   bounded by `10000a^3(3/4)^(a-1)`.  The asymptotic-bounds script verifies
   that their sum is below the reserve `1/(4a)` for all `a >= 200`.

Together with the exact finite check for `a <= 200`, this gives a plausible
complete proof of Route-2 for the `S(2^a,8)` lane.  It still needs an
independent audit before being treated as manuscript-grade.

## General Fixed-`r` Hub-Off Certificate

I added a parameterized hub-off certificate:

```bash
python3 gpt_attack/fixed_r_huboff_certificate.py \
  --r 8 --threshold 200 --margin-denom 10 --reserve-denom 4
```

For a general fixed `r`, it:

1. computes the stabilized shifts `D_q = 3m-2a` at the threshold;
2. checks, residue class by residue class, the hub-off margins
   `(F_m-F_{m-1})/F_m` and `(F_m-F_{m+1})/F_m`;
3. checks the hub-off Route-2 reserve;
4. checks `3/4 <= lambda_0 <= 2`;
5. proves each inequality by shifting `t` to the threshold and verifying
   positive coefficients.

The implementation now uses normalized coefficient ratios
`F_{m+s}/(2^m binom(a,m))` and constructs the reserve numerator after shifting
`t` to the threshold.  This avoids the direct Sympy expansion bottleneck that
blocked the first `r=20` attempt.

Runs completed:

```bash
python3 gpt_attack/fixed_r_huboff_certificate.py --r 4  --threshold 200 --margin-denom 1000 --reserve-denom 1000
python3 gpt_attack/fixed_r_huboff_certificate.py --r 8  --threshold 200 --margin-denom 10   --reserve-denom 4
python3 gpt_attack/fixed_r_huboff_certificate.py --r 12 --threshold 200 --margin-denom 1000 --reserve-denom 1000
python3 gpt_attack/fixed_r_huboff_certificate.py --r 16 --threshold 200 --margin-denom 1000 --reserve-denom 1000
python3 gpt_attack/fixed_r_huboff_certificate.py --r 20 --threshold 200 --margin-denom 1000 --reserve-denom 1000
python3 gpt_attack/fixed_r_huboff_certificate.py --r 24 --threshold 200 --margin-denom 1000 --reserve-denom 1000
python3 gpt_attack/fixed_r_huboff_certificate.py --r 32 --threshold 200 --margin-denom 1000 --reserve-denom 1000
python3 gpt_attack/fixed_r_huboff_certificate.py --r 40 --threshold 200 --margin-denom 1000 --reserve-denom 1000
python3 gpt_attack/fixed_r_huboff_certificate.py --r 60 --threshold 200 --margin-denom 1000 --reserve-denom 1000
python3 gpt_attack/fixed_r_huboff_certificate.py --r 80 --threshold 200 --margin-denom 1000 --reserve-denom 1000
```

All selected runs completed with every symbolic inequality certified.  The
`r=80` run is noticeably slower than the lower lanes but completes under the
optimized exact-polynomial method.

At this checkpoint, before adding the hub-on certificates, the practical
frontier was:

```text
hub-off symbolic certificate works for selected fixed-r lanes through r=80;
the remaining fixed-r proof work is hub-on perturbation/mode-localization,
not the hub-off reserve.
```

## General Fixed-`r` Hub-On Certificates

The hub-on branch is now split into two exact certificate scripts.

Mode domination:

```bash
python3 gpt_attack/fixed_r_hubon_mode_certificate.py --r 80 --threshold 200 --margin-denom 1000
```

For each residue class, the script chooses a witness term `p_s x^s` in
`P_r(x)` and uses:

```text
G_k <= P_{r-1}(1) 2^a,
F_m >= p_s 2^(m-s) binom(a,m-s).
```

It checks at the threshold that:

```text
2 max_k G_k/F_m < 1/(D a),
```

where `D` is the hub-off margin denominator, and proves the bound decreases
under `a -> a+3` in the same residue class.  This gives full-polynomial mode
localization once the hub-off certificate supplies the margin.

Route-2 perturbation:

```bash
python3 gpt_attack/fixed_r_hubon_route2_perturbation.py --r 80 --threshold 200 --reserve-denom 1000
```

The perturbation certificate uses two coarse bounds:

```text
|lambda-lambda_0| contribution <= 2 N^2 T,
hub-on mixture in B <= N C_r (3/4)^(a-1),
```

where `T` is the certified bound on `max_k G_k/F_m`, `N=a+r`, and

```text
C_r = 2 P_{r-1}(2)/P_r(1/2).
```

The `3/4` base is justified by the strengthened interval
`lambda_0 >= 3/4` plus the tiny mode-branch error, which keeps the full
fugacity in `[1/2,2]`.

Both scripts were run successfully for:

```text
r in {4,8,12,16,20,24,32,40,60,80}
```

using the same threshold `a >= 200`, with the strong `r=8` denominators
`D=10`, `R=4` and weak `D=R=1000` for the other sampled lanes.  For `r=80`,
the perturbation certificate reports the largest sampled mixture bound as
about `8.37e-10`, still far below the weak reserve `5e-6`.

Current practical frontier:

```text
For selected fixed-r lanes through r=80, the proof template now has certified
hub-off reserve, full-mode localization, and full Route-2 perturbation control
for a >= 200.  What remains is packaging this as a clean parametric theorem
for arbitrary fixed r, plus exact finite checks below the chosen threshold.
```

## Exact Finite Checks Below the Threshold

I added an exact integer finite checker:

```bash
python3 gpt_attack/fixed_r_finite_route2_check.py \
  --r-values selected \
  --a-max 199 \
  --out results/fixed_r_finite_route2_selected_a199.json
```

It clears denominators at `lambda=i_{m-1}(T)/i_m(T)`.  For
`B(x)=sum b_k x^k` and `lambda=L/M`, Route-2 is checked by the exact integer
sign:

```text
sum_k (2k-2m+3) b_k L^k M^(d-k) >= 0.
```

The stronger threshold

```text
mu_B(lambda) >= m - 1 - lambda/(1+lambda)
```

is checked by the analogous cleared-denominator integer sign.

Run result:

```text
checked records: 1990
skipped records: 0
Route-2 failures: 0
stronger-threshold failures: 0
global min Route-2: r=8, a=198, slack ~= 0.1684594347526992
global min stronger: r=8, a=198, slack ~= 0.1682996084886373
```

Thus, for the sampled lanes

```text
r in {4,8,12,16,20,24,32,40,60,80},
```

Route-2 is now certified by exact finite checking for `a <= 199` and by the
hub-off/hub-on asymptotic certificates for `a >= 200`.
