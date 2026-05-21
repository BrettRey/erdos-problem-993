# Fixed-`r` Route-2 Theorem Plan After Reading the Unit-Distance CoT

## Purpose

The OpenAI unit-distance CoT is useful here as a proof-architecture example,
not because its algebraic number theory transfers directly to Erdos 993.

The main lesson is procedural:

```text
keep the quantifiers exact;
search for the compression mechanism;
kill insufficient routes explicitly;
turn computations into a theorem-shaped certificate;
track the remaining hidden obstructions.
```

This note applies that pattern to the fixed-`r` Route-2 spider lanes

```text
S(2^a,r),    a -> infinity, r fixed.
```

Sources:

- OpenAI announcement:
  https://openai.com/index/model-disproves-discrete-geometry-conjecture/
- Rewritten CoT summary:
  https://cdn.openai.com/pdf/1625eff6-5ac1-40d8-b1db-5d5cf925de8b/unit-distance-cot.pdf

## Target Statement

For fixed `r >= 2`, let

```text
T_{a,r} = S(2^a,r)
```

and remove a length-2 arm to get `B_{a,r}`.  Let `m=m(a,r)` be the leftmost
mode of `I(T_{a,r})`, and set

```text
lambda = i_{m-1}(T_{a,r}) / i_m(T_{a,r}).
```

The fixed-`r` Route-2 target is:

```text
mu_{B_{a,r}}(lambda) >= m - 3/2
```

for every `a >= 1`.

The theorem shape should be:

```text
For every fixed r, there is an explicit threshold A(r) such that:

1. exact finite checking proves the target for 1 <= a < A(r);
2. symbolic hub-off and hub-on certificates prove the target for a >= A(r).
```

This is deliberately not a uniform-in-`r` statement.  The evidence already
shows that a universal threshold such as `a=200` is false for the current proof
template.

## Exact Polynomial Split

The key fixed-`r` identity is:

```text
I(T_{a,r};x) = P_r(x)(1+2x)^a + xP_{r-1}(x)(1+x)^a.
```

Write:

```text
F_{a,r}(x) = P_r(x)(1+2x)^a
G_{a,r}(x) = xP_{r-1}(x)(1+x)^a.
```

The compression mechanism is that `G_{a,r}` is exponentially small near the
`F_{a,r}` saddle once `r` is fixed and `a` is large enough.  The fixed
path-polynomial convolution `P_r` shifts the mode by a residue-class constant,
but it does not change the binomial spine.

The expected mode form is:

```text
m(a,r) = (2a + D_{r,q})/3,    q = a mod 3,
```

for all `a >= A(r)`.

This is the analogue of the unit-distance CoT's "find the actual compression
mechanism" step: many directions in a low-rank additive object there; here a
fixed path-polynomial convolution around a binomial saddle.

## Certificate Components

### 1. Mode-Shift Certificate

For each residue class `q in {0,1,2}`, produce:

```text
D_{r,q}
A_mode(r)
```

such that for all `a >= A_mode(r)`, `a == q mod 3`,

```text
c_m(T_{a,r}) >= c_{m-1}(T_{a,r})
c_m(T_{a,r}) >= c_{m+1}(T_{a,r})
```

with `m=(2a+D_{r,q})/3`.

Current scripts find these shifts by exact coefficient checks at a proposed
threshold.  A theorem-ready version should derive `D_{r,q}` from the fixed
convolution `P_r(x)(1+2x)^a` and then show the exponentially small `G` term
cannot move the mode.

### 2. Hub-Off Reserve Certificate

Set

```text
lambda_0 = F_{m-1}/F_m.
```

The hub-off branch uses the post-removal approximation

```text
B^0_{a,r}(x) = P_r(x)(1+2x)^(a-1).
```

The current certificate proves the stronger reserve

```text
mu_{B^0_{a,r}}(lambda_0) >= m - 4/3 + 1/(1000a)
```

and also proves a compact fugacity interval:

```text
3/4 <= lambda_0 <= 2.
```

The reserve is checked exactly after shifting `t` by the threshold in each
residue class.  The path recurrence

```text
E_n = Z^deg(P_n) P_n(L/Z)
N_n = Z^deg(P_n) (L/Z)P'_n(L/Z)
```

keeps this from expanding powers of `L` and `Z` directly.

Current implementation status:

```text
Python symbolic recurrence: good through r=120.
GMP threaded exact recurrence: good through r=280.
Next unrun stress point: r=320, threshold=500.
```

The theorem gap is not numerical failure.  It is replacing dense polynomial
positivity checks by either:

```text
1. a readable inequality for all fixed r; or
2. a certified algorithm with a proof of termination/correctness for each r.
```

### 3. Hub-On Perturbation Certificate

The full polynomial differs from the hub-off approximation by the hub-on term.
The current proof bounds two effects:

```text
fugacity shift:     |lambda - lambda_0|
mixture shift:      G_B(lambda)/F_B(lambda)
```

using:

```text
G/F tail domination near the target mode;
lambda in [3/4,2];
degree and variance crude bounds;
an explicit threshold check plus monotonicity in a -> a+3.
```

This part scales better than hub-off.  Current grid:

```text
r=80  -> A=200
r=100 -> A=200
r=120 -> A=250
r=160 -> A=300
r=200 -> A=350
r=240 -> A=400
r=280 -> A=450
r=320 -> A=500
r=400 -> A=650
```

This suggests `A(r)` grows roughly linearly for the tested range, with
threshold jumps caused by the conservative mixture bound.

### 4. Finite Exact Certificate

The finite side now evaluates

```text
B(lambda), lambda B'(lambda)
```

directly from

```text
B(x) = P_r(x)(1+2x)^(a-1) + xP_{r-1}(x)(1+x)^(a-1),
```

clearing one common denominator.  This avoids coefficient-by-coefficient
summation at `lambda=L/M`.

Current finite coverage:

```text
selected r <= 80:       a <= 199, 1990 records, 0 failures
r in {100,120,160,200}: a <= 349, 1396 records, 0 failures
r=240:                  a <= 399, 399 records, 0 failures
r=280:                  a <= 449, 449 records, 0 failures
r=320:                  a <= 499, 499 records, 0 failures
```

The finite checker is no longer the main bottleneck for the tested fixed-`r`
lanes.

## Complete and Partial Lanes

Complete under the current certificate split:

```text
sampled set: {4,8,12,16,20,24,32,40,60,80}, A=200
r=100, A=200
r=120, A=250
r=160, A=300
r=200, A=350
r=240, A=400
r=280, A=450
```

Partial:

```text
r=320:
  finite exact check passes for a <= 499;
  hub-on certificates pass for a >= 500;
  hub-off reserve at a >= 500 has not yet been run.

r=400:
  hub-on certificates pass for a >= 650;
  finite and hub-off pieces not yet packaged.
```

## Dead or Insufficient Routes

These should not be mistaken for theorem paths:

```text
1. More finite lane checks alone.
   Useful for stress-testing, but they do not prove arbitrary fixed r.

2. A uniform threshold such as a=200.
   Already false for the conservative hub-on certificate.

3. Generic unimodality/log-concavity language.
   The proof needs the specific tie fugacity and Route-2 inequality.

4. Treating the GMP helper as the theorem.
   It is an exact certificate engine, but the paper needs a clean statement of
   what it certifies and why that certificate is enough.

5. Re-running old SCC or obsolete notes.
   Those are marked dead ends elsewhere and do not touch the current fixed-r
   proof gap.
```

## Risk Register

### R1. Stabilized shifts

Current scripts extract `D_{r,q}` from the proposed threshold.  A theorem needs
an effective rule or certificate for these shifts.

Acceptable next step:

```text
Define a shift certificate as the finite list of exact rational functions
for the two boundary differences in each residue class, shifted by A(r), with
positive coefficients.
```

Better theorem step:

```text
Derive D_{r,q} from the first nonzero terms in the local expansion of
P_r(x)(1+2x)^a around the binomial saddle.
```

### R2. Hub-off reserve positivity

Current positivity checks are exact but dense.  For `r=280`, the reserve
polynomial degrees are about `14382`.  This is fine as a certificate artifact
for fixed lanes, but not a conceptual proof for arbitrary fixed `r`.

Possible theorem paths:

```text
1. Prove an analytic lower bound for the hub-off reserve using saddle-point
   expansion, uniform for fixed r once a >= C r.

2. Prove coefficient positivity by induction through the path recurrence after
   shifting t by A(r).

3. State a certified-algorithm theorem: given r, the algorithm constructs
   rational functions whose shifted numerators/denominators have positive
   coefficients, and this proves the lane.  This still needs a termination
   theorem or an explicit proof for each lane used in the paper.
```

### R3. Hub-on constants

Hub-on is conservative but stable.  The key risk is that crude constants make
`A(r)` larger than necessary, increasing finite burden and hub-off polynomial
degree.

This is engineering risk, not current mathematical failure.

### R4. Paper scope

There are two different publishable levels:

```text
Level 1: certified finite family of fixed-r lanes, including r up to 280.
Level 2: theorem for every fixed r.
```

Do not blur these.  The current evidence supports Level 2, but the written
proof is not there yet.

## Next Moves

### Best immediate proof move

Write a formal "certificate lemma" for fixed `r`:

```text
If the following four exact certificate files/checks pass for a given r and A:
  finite exact range;
  mode shift;
  hub-off reserve;
  hub-on perturbation;
then Route-2 holds for S(2^a,r) for all a.
```

This converts the computation into a checkable theorem schema.

This has been drafted in:

```text
notes/fixed_r_certificate_lemma_2026-05-21.md
```

The draft identifies the fugacity-shift sublemma as the least formalized link
in the current certificate chain.

That link and the adjacent-to-global mode-margin link have now been separated
into:

```text
notes/fixed_r_fugacity_shift_sublemma_2026-05-21.md
notes/fixed_r_global_f_mode_margin_sublemma_2026-05-21.md
```

The remaining local certificate justifications have also been drafted:

```text
notes/fixed_r_shifted_coefficient_positivity_2026-05-21.md
notes/fixed_r_huboff_reserve_recurrence_certificate_2026-05-21.md
```

The fixed-`r` mode-shift rule has also been separated into:

```text
notes/fixed_r_shift_rule_from_path_moments_2026-05-21.md
```

It derives the first-order candidate interval
`3M_1-1 <= D <= 3M_1+2`; computation through `r=2000` found only the
boundary cases `r=2,3`.  The note now proves the boundary classification:
using `P_r(1)=F_{r+2}` and
`P'_r(1)=(rL_{r+1}+2F_r)/5`, boundary would force
`F_{r+2} | 3(r+2)`, impossible for `r>=7`, with `r=4,5,6` checked directly.

The positivity target for `C_{r,q}` has been separated into:

```text
notes/fixed_r_c_positivity_reduction_2026-05-21.md
```

It shows that, after the shift interval, it is enough to prove the path moment
inequality `V+3K_3>=0`, equivalently a finite trigonometric sum over the roots
of `P_r`.  The note now proves that inequality by solving the raw path-moment
recurrences and reducing the result to a Fibonacci/Lucas expression.

### Best analytic move

Derive the hub-off reserve asymptotic:

```text
mu_{B^0}(lambda_0) - (m - 4/3)
```

for

```text
F=P_r(x)(1+2x)^a,
m=(2a+D_{r,q})/3.
```

Goal:

```text
mu_{B^0}(lambda_0) - (m - 4/3) >= c_r/a
```

or preferably a lower bound with `c_r` controlled well enough that the hub-on
perturbation can be made smaller by taking `A(r)`.

This has now been started in:

```text
gpt_attack/fixed_r_huboff_asymptotic.py
notes/fixed_r_huboff_asymptotic_reserve_2026-05-21.md
```

The expansion finds:

```text
mu_{B^0}(lambda_0) - (m - 4/3)
  = C_{r,q}/a + O_r(a^-2)
```

with positive tested constants through `r=400`.  The arbitrary fixed-`r`
reserve problem is now to prove `C_{r,q}>0` for the stabilized shifts and make
the `O_r(a^-2)` remainder effective.

The leading constant is now identified in closed form.  If
`P_r(x)=sum_s p_s x^s`, `Pr(S=s)=p_s/P_r(1)`, `M_j=E[S^j]`, and
`delta=D_{r,q}/3`, then

```text
C_{r,q}
 = (-24 delta + 54 M_1^3 - 9 M_1^2 - 81 M_1 M_2
    + 24 M_1 + 9 M_2 + 27 M_3 + 16)/12
 = (9 Var(S) + 27 E[(S-M_1)^3] + 24(M_1-delta) + 16)/12.
```

So the analytic subproblem has tightened again: handle the effective
`O_r(a^-2)` remainder, write the adjacent-ratio remainder that turns the
first-order shift rule into an eventual mode-shift theorem, and then connect
the already-separated hub-on domination lemma.

The eventual fixed-`r` theorem schema is now written in:

```text
notes/fixed_r_eventual_route2_theorem_schema_2026-05-21.md
```

It packages the upgraded state as a computable-threshold theorem: for every
fixed `r`, rational-function root isolation plus exponential hub-on
domination should give an `A(r)` after which Route-2 holds, with exact finite
checking below `A(r)`.

The corresponding threshold probe is now implemented in:

```text
gpt_attack/fixed_r_effective_threshold_probe.py
```

Exact root isolation gives `A=15` for `r=4`, `A=27` for `r=8`, and `A=97` for
`r=20`.  The probe confirms the theorem schema is operational, but also shows
that naive exact isolation slows by `r=24`; a Cauchy fallback proves
computability but gives impractically huge thresholds at `r=40`.

The threshold probe has now been improved with shifted-coefficient positivity.
This method gives practical hub-off/mode/lambda thresholds:

```text
r=4 -> 12, r=8 -> 24, r=20 -> 94, r=24 -> 22,
r=40 -> 26, r=80 -> 44, r=120 -> 62.
```

The effective hub-on probe is:

```text
gpt_attack/fixed_r_hubon_effective_threshold_probe.py
```

Using half of the actual first-order mode margins and half of `C_{r,q}`, with
corrected `a * perturbation` monotonicity, the first passing thresholds in a
step-5 grid are:

```text
r=80 -> 145, r=120 -> 200, r=160 -> 255, r=200 -> 305,
r=240 -> 360, r=280 -> 415, r=320 -> 470, r=400 -> 575.
```

This identifies hub-on perturbation, not hub-off reserve, as the practical
threshold driver in the tested range.

### Best stress test

Run the remaining hub-off certificate:

```bash
python3 gpt_attack/fixed_r_huboff_cpp_certificate.py \
  --r 320 --threshold 500 --margin-denom 1000 --reserve-denom 1000 \
  --reserve-threads 8
```

This is not the main theorem path, but it tests whether the certificate engine
continues to scale past the current frontier.

## Current Call

Push on, but pivot upward.

The computation has done its job: it found a robust fixed-`r` pattern and
removed several false bottlenecks.  The next deliverable should be a formal
certificate lemma plus a hub-off asymptotic/reserve lemma.  More lane runs are
useful only insofar as they stress-test those lemmas.
