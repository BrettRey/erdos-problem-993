# Broom unimodality proof attempt

Goal: prove unimodality for all brooms B(p,s) with p >= 2, s >= 0.

## 1. Closed form and coefficient notation

Let
  A(x) = I(P_{p-1}; x) = sum_{t>=0} a_t x^t,
  B(x) = I(P_{p-2}; x) = sum_{t>=0} b_t x^t,
where
  a_t = C(p - t, t),
  b_t = C(p - 1 - t, t),
with support t <= floor(p/2) and t <= floor((p-1)/2) respectively.

For the broom B(p,s),
  I(B(p,s); x) = (1+x)^s A(x) + x B(x).

Write c_k = [x^k] I(B(p,s); x). Then
  c_k = d_k + e_k,
  d_k = sum_{j=0}^{min(s,k)} C(s,j) a_{k-j},
  e_k = b_{k-1} (with b_{-1} = 0).

So d_k is the binomial convolution of a_t with C(s,j), and e_k is a short
positive tail supported on k = 1..S where S = floor((p-1)/2) + 1.

## 2. Structural facts

F1. A(x) and B(x) are Fibonacci polynomials. Their zeros are real and
    negative, so (a_t) and (b_t) are PF-infinity sequences (in particular
    log-concave with no internal zeros).

F2. Binomial coefficients form a PF-infinity sequence. By Hoggar's theorem,
    the convolution of PF-infinity sequences is PF-infinity. Therefore d_k is
    log-concave and unimodal.

F3. e_k is log-concave and unimodal (it is just b_{k-1}).

These facts do NOT imply c_k is log-concave or unimodal in general, because
sums of log-concave sequences can fail to be unimodal.

## 3. A sufficient monotonicity inequality (case s >= p)

Define S = floor((p-1)/2) + 1, the last index where e_k is nonzero.
If we can show

  (I)  d_{k+1} - d_k >= e_k - e_{k+1}  for all 0 <= k <= S - 1,

then c_k = d_k + e_k is nondecreasing up to k = S. Since e_k = 0 for k > S,
we then have c_k = d_k for all k >= S. Because d_k is unimodal, c_k is
unimodal.

Empirical check: (I) holds for all p <= 60, s <= 120 with s >= p.

### Proof sketch for s >= p (fully symbolic)

Fix p >= 2 and s >= p. Let S = floor((p-1)/2)+1, so e_k = 0 for k >= S+1.
We show c_{k+1} >= c_k for all 0 <= k <= S-1.

Rewrite the d-difference by changing variables t = k+1-j:
  d_{k+1} - d_k = sum_t a_t [C(s, k+1-t) - C(s, k-t)].

For k <= S-1 we have k+1 <= S = ceil(p/2) <= ceil(s/2), hence for every t in
the summation, r = k+1-t satisfies r <= k+1 <= ceil(s/2). The binomial row
is nondecreasing up to r = ceil(s/2), so
  C(s, r) - C(s, r-1) >= 0.
Thus all terms are nonnegative.

In particular, keeping only the t=1 term gives
  d_{k+1} - d_k >= a_1 [C(s,k) - C(s,k-1)] = (p-1)[C(s,k) - C(s,k-1)].

Since s >= p and k <= floor(p/2), the function
  f(s) = C(s,k) - C(s,k-1) = C(s,k-1) * (s-2k+1)/k
is increasing in s on s >= 2k-1, hence
  C(s,k) - C(s,k-1) >= C(p,k) - C(p,k-1).

Also
  b_{k-1} - b_k = C(p-k, k-1) - C(p-k-1, k) <= C(p-k, k-1) <= C(p,k-1).

Therefore it suffices to show
  (p-1)[C(p,k) - C(p,k-1)] >= C(p,k-1).
Using C(p,k) = C(p,k-1) * (p-k+1)/k, this reduces to
  (p-1)(p-2k+1) >= k,
which holds for k <= floor(p/2).

Conclusion: c_{k+1} >= c_k for all 0 <= k <= S-1 when s >= p.
Since e_k = 0 for k >= S+1, the tail equals d_k, which is unimodal
(convolution of unimodal sequences). Hence c_k is unimodal for all s >= p.

## 4. Literature closure (stronger than needed)

Theorem. For every p >= 2 and s >= 0, the broom B(p,s) has a log-concave
independence sequence; in particular it is unimodal.

Proof.
1) Log-concavity implies unimodality.
   If a_k^2 >= a_{k-1} a_{k+1} for all k, then r_k = a_k / a_{k-1} is
   nonincreasing, hence crosses 1 at most once. So the sequence increases
   and then decreases.

2) Reduce brooms to paths or spiders.
   If s = 0, then B(p,0) = P_p (a path).
   If s = 1, then B(p,1) = P_{p+1}.
   If s >= 2, then B(p,s) has a unique vertex of degree >= 3 (the attachment
   vertex), with one long leg and s legs of length 1, i.e., a spider.

3) Spiders are strongly log-concave.
   Li–Li–Yang–Zhang (2025, arXiv:2501.04245) prove that every spider has a
   strongly log-concave independence polynomial (their Theorem 3.1).
   Strong log-concavity implies a_k^2 >= a_{k-1} a_{k+1} by taking the
   adjacent-index specialization.

4) Apply to brooms.
   For s >= 2, B(p,s) is a spider, hence log-concave.
   For s <= 1, B(p,s) is a path; path independence polynomials are real-rooted
   and therefore log-concave.
   Thus B(p,s) is log-concave for all s >= 0, hence unimodal.  □

The remainder of this section is a superseded attempt at an elementary
closure for s < p; it is retained for historical context.

## 5. Remaining case s < p (superseded attempt)

The inequality (I) above need not hold for small s. A second argument is
needed for s < p.

Two possible routes:

(A) Induction in s with a unimodality-preserving lemma. The recurrence
    I(B(p,s+1)) = I(B(p,s)) + x (1+x)^s A(x)
    adds a log-concave sequence to c. If one can show the added sequence has
    its mode >= the mode of I(B(p,s)) and is sufficiently "flat" near the
    peak, unimodality is preserved.

(B) Direct ratio control. Define r_k = c_{k+1}/c_k. Show r_k crosses 1 at
    most once by bounding r_k between ratios for d_k and e_k. This likely
    needs sharp inequalities for c_{k+1}c_{k-1} - c_k^2.

Empirical check: no log-concavity failures for p <= 60, s <= 120.
This suggests the stronger statement "all brooms are log-concave" may be true.

Additional empirical check: inequality (I) holds for all p <= 60, s <= 120
whenever s >= ceil(p/2). This suggests the s >= p proof can likely be
extended to s >= p/2 with a sharper lower bound (using a different t in the
sum for d_{k+1}-d_k).

## 5b. Stronger sufficient inequality (tested for s >= p/2)

Define A_k = e_k (independent sets of size k that include the hub) and
B_k = d_k (independent sets of size k that exclude the hub). Then
c_k = A_k + B_k and

  c_{k+1} - c_k = (B_{k+1} - B_k) - (A_k - A_{k+1}).

A sufficient (stronger) condition is

  (III)  B_{k+1} - B_k >= A_k,  for all 0 <= k <= S-1.

If (III) holds, then c_{k+1} - c_k >= A_{k+1} >= 0, so the sequence is
nondecreasing up to k = S. (The tail equals B_k and is unimodal.)

Empirical check: (III) holds for all p <= 60, s <= 120 whenever
s >= ceil(p/2). In particular, this is stronger than (I) in that regime.

Potential path: use the form
  B_{k+1} - B_k = sum_j C(s,j) (a_{k+1-j} - a_{k-j})
and show the sum is >= b_{k-1} for k <= S-1 when s >= p/2.
The positive term with j = k-1 contributes
  C(s, k-1) (a_2 - a_1),
which is large for s >= p/2. The negative terms correspond to small j and
can potentially be bounded by a polynomial in p, while C(s,k-1) grows
exponentially in s once k is in the central range.

Combinatorial injection idea: if we can inject all size-k independent sets
into size-(k+1) independent sets in a way that treats hub-containing sets
separately, then (III) follows directly. This likely requires at least three
reserved leaves and a collision-free mapping (not yet formalized).

## 5c. Attempted analytic split (r <= s/2 vs r > s/2)

Write
  d_{k+1} - d_k = sum_{r=0}^{k+1} a_{k+1-r} * DeltaC(s,r),
where DeltaC(s,r) = C(s,r) - C(s,r-1).

Let r* = floor(s/2) and t0 = k+1 - r*. For k <= S-1 and s >= p/2,
we have t0 <= p/4, which is below the mode of a_t. Thus a_t is increasing
on [0, t0].

Split the sum:
  P = sum_{r=0}^{r*} a_{k+1-r} * DeltaC(s,r)   (nonnegative terms),
  N = sum_{r=r*+1}^{k+1} a_{k+1-r} * DeltaC(s,r) (negative terms).

If a_t is increasing on [0,t0], then for r <= r* we have t >= t0 and
a_{k+1-r} >= a_{t0}, so
  P >= a_{t0} * C(s, r*).

For r > r* we have t < t0, and since a_t >= 1,
  N >= 1 * (C(s, k+1) - C(s, r*))   (telescoping sum of DeltaC).

Therefore
  d_{k+1} - d_k >= C(s, k+1) + (a_{t0} - 1) C(s, r*).

This bound is too weak when t0 = 0 (k+1 <= r*). Empirically, the inequality
d_{k+1} - d_k >= b_{k-1} still holds for s >= p/2, but a tighter lower bound
is required in the regime k <= s/2 - 1. The likely fix is to keep several
positive terms (t = 0,1,2,...) rather than bounding by a single term.

## 5d. Multi-term lower bound for the early regime (k <= s/2 - 1)

When k <= s/2 - 1, all DeltaC(s,r) for r <= k+1 are nonnegative, so we can
drop all but a few positive terms:

  d_{k+1} - d_k
    = sum_t a_t DeltaC(s, k+1-t)
    >= a_0 DeltaC(s,k+1) + a_1 DeltaC(s,k) + a_2 DeltaC(s,k-1).

Since
  a_0 = 1,
  a_1 = p-1,
  a_2 = C(p-2,2),
this yields the explicit bound

  d_{k+1} - d_k >= C(s,k+1)
    + (p-1)(C(s,k) - C(s,k-1))
    + C(p-2,2)(C(s,k-1) - C(s,k-2)).

Empirically, for all p <= 120 and s <= 240 with s >= ceil(p/2), this lower
bound is >= b_{k-1} - b_k throughout the early regime k <= s/2 - 1.

If this inequality can be proved symbolically (for s >= p/2, k <= s/2 - 1),
then monotonicity follows for all k in that early regime. The remaining
regime is k >= s/2, where DeltaC changes sign and the split argument of 4c
must be tightened (or replaced).

## 5e. New empirical pattern: the minimum of F_k occurs at k = 0

Define
  F_k(p,s) = (d_{k+1} - d_k) - b_{k-1}.
Then the sufficient inequality (III) is F_k(p,s) >= 0 for all k <= S-1.

Empirically (checked by brute force for all p <= 200 and all s in [ceil(p/2), p]),
the minimum of F_k(p,s) over k in [0, S-1] always occurs at k = 0.

At k = 0,
  F_0 = d_1 - d_0 = (p-1) + s - 1 = p + s - 2.
So if one can prove
  F_k(p,s) >= F_0(p,s)  for all k <= S-1 and s >= p/2,
then (III) follows immediately.

This suggests a promising proof route:
  1) Express F_k using increments of a_t:
       d_{k+1} - d_k = sum_{j=0}^s C(s,j) (a_{k+1-j} - a_{k-j})
     where a_t = C(p-t, t) and the increment sequence
       Delta a_t = a_{t+1} - a_t
     is positive for t <= floor((p-1)/3) and negative afterward
     (the mode of a_t is around p/3).
  2) For s >= p/2 and k <= S, the binomial weights C(s,j) in the sum are
     largest near j = s/2, so the effective t = k - j concentrates in
     the range t <= p/4, where Delta a_t is strictly positive.
     The negative Delta a_t terms correspond to small j (tail of the binomial),
     which carry exponentially smaller total weight.
  3) Formalize this with a binomial-tail bound:
       sum_{j=0}^J C(s,j) << C(s, floor(s/2))
     for J ≈ p/6, while the positive contributions from the central window
     exceed p+s-2.

So far this is only an empirical pattern + proof sketch. It is strong enough
to justify a targeted analytic effort: show F_k is bounded below by F_0 via
central-binomial dominance once s >= p/2, then check finitely many small p
to close the remaining cases.

## 5f. Useful identity and a clean early-regime lemma

The path coefficients satisfy
  a_t = b_t + b_{t-1},
so the increment is
  Delta a_t = a_t - a_{t-1} = b_t - b_{t-2}.

Using the convolution form
  d_{k+1} - d_k = sum_t Delta a_t * C(s, k+1 - t),
we obtain the cleaner identity

  d_{k+1} - d_k
    = sum_t b_t [ C(s, k+1 - t) - C(s, k-1 - t) ].

Define w_r = C(s,r) - C(s,r-2), with the convention C(s,r)=0 for r<0.
Then w_r >= 0 for all r <= floor(s/2), since C(s,r) is increasing there.

Lemma (early regime): If k <= s/2 - 1 and k >= 1, then
  d_{k+1} - d_k >= b_{k-1} * (C(s,2) - 1).
In particular, for s >= 3 this implies
  d_{k+1} - d_k >= b_{k-1},
so the strong inequality (III) holds throughout the early regime.

Proof sketch: for k <= s/2 - 1, we have k+1 <= s/2, hence for every t >= 0
the index r = k+1 - t satisfies r <= s/2 (or r < 0, where w_r = 0). Thus
all terms in the sum are nonnegative. The term t = k-1 yields r = 2 and
contributes b_{k-1} (C(s,2) - C(s,0)) = b_{k-1}(C(s,2) - 1), giving the bound.

This reduces the remaining analytic work for s >= p/2 to the late regime
k >= s/2, where w_r changes sign and cancellation must be controlled.

## 5g. Late-regime reduction via the paired w_r identity

Recall
  d_{k+1} - d_k = sum_t b_t [C(s, k+1 - t) - C(s, k-1 - t)]
and define w_r = C(s,r) - C(s,r-2), with C(s,r)=0 for r<0.
Then
  d_{k+1} - d_k = sum_r w_r * b_{k+1-r},
and w_r satisfies the antisymmetry
  w_r = - w_{s-r+2}.

Let r* = floor(s/2) + 1. Pairing r with s-r+2 gives
  d_{k+1} - d_k
    = sum_{r=0}^{r*} w_r [ b_{k+1-r} - b_{k-s-1+r} ].

Late regime setup: assume s >= ceil(p/2) and k >= s/2.
Then k <= S-1 <= floor((p-1)/2), which implies k <= s-1 (since s > S-1).
Hence k-s+1 <= 0 and, in particular, b_{k-s+1} in {0,1}.

For each r in [0, r*], define
  t  = k+1-r,
  t' = k-s-1+r.
Then t >= t' and
  t + t' = 2k - s <= (p-1)/2 - 1.

Empirically (checked for all p <= 200, s in [ceil(p/2), p], and all admissible k,r),
we have
  b_t >= b_{t'}  whenever t >= t' and t + t' <= (p-1)/2.

If this binomial-coefficient inequality can be proved in general, then each
paired term is nonnegative and we may keep only r=2:

  d_{k+1} - d_k >= w_2 [b_{k-1} - b_{k-s+1}].

Since k <= s-1, we have b_{k-s+1} in {0,1}, and w_2 = C(s,2) - 1 >= 1
for s >= 3. This yields d_{k+1} - d_k >= b_{k-1}, i.e. (III), throughout
the late regime.

Thus the remaining analytic task is to prove:

  Lemma (binomial inequality).
  For b_t = C(p-1-t, t), if t >= t' and t + t' <= (p-1)/2, then b_t >= b_{t'}.

This is now a purely binomial-coefficient statement. A possible route is to
analyze the ratio
  r_i = b_i / b_{i-1} = (p-2i)(p+1-2i) / (i(p-i)),
and show that for fixed M = t+t' <= (p-1)/2, the product
  prod_{i=t'+1}^t r_i >= 1.
The numerical evidence for this inequality is strong; formalizing it would
close the late-regime proof.

## 6. Next steps (historical)

1. Formalize (I) for s >= p.
   - Prove d_{k+1} - d_k is positive for k <= S and obtain a lower bound.
   - Prove b_{k-1} - b_k is bounded above by that lower bound.

2. Explore s < p.
   - Attempt to show log-concavity directly from c_k^2 - c_{k-1}c_{k+1}.
   - If (I) can be proved for s >= p, check whether a symmetric inequality
     in the regime p >= 2s can also be shown (would reduce to a finite band).

3. If a clean inequality emerges, integrate into paper as a theorem.

4. Strengthen the s-regime.
   - Empirically, (I) holds for s >= p/2. Try bounding d_{k+1}-d_k by a
     term with r = floor(s/2), where binomial differences are largest, and
     compare to b_{k-1}-b_k.
   - If this works, only the regime s < p/2 remains.

5. Explore (III) as the main target.
   - Attempt to prove (III) for s >= p/2 by bounding negative terms in the
     sum for B_{k+1}-B_k and comparing to b_{k-1}.
   - If successful, this gives a clean monotonicity proof up to k = S with
     a simpler inequality than (I).
