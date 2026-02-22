You do not have file access. Treat everything below as complete ground truth.
Do not run code or discuss unrelated topics. Stay only in this notation.

Goal for this round: close the exact symbolic gap for Route 1, or reduce it to one minimal new lemma with fully explicit assumptions.

=== Setup and notation ===
- Tree polynomial: I(T;x) = sum_k i_k(T) x^k.
- m = leftmost mode index of I(T).
- Canonical bridge setup:
  - choose leaf l with minimum support degree;
  - support s has deg(s)=2;
  - u is the other neighbor of s;
  - B = T - {l,s};
  - P = dp_B[u][0], Q = dp_B[u][1], I(B)=P+Q.
- Index shorthand:
  - p0=p_{m-2}, p1=p_{m-1}, pm=p_m;
  - q0=q_{m-2}, q1=q_{m-1}, qm=q_m;
  - b0=b_{m-2}, b1=b_{m-1}, b2=b_m.
- Fugacity and means:
  - lambda = lambda_m(T) = i_{m-1}(T)/i_m(T), a=lambda/(1+lambda).
  - mu_P(lambda), mu_B(lambda).
- Transfer quantities:
  - D = mu_B(lambda)-mu_P(lambda)
  - exact_slack_B = mu_B(lambda) - (m-1-a)
  - exact_excess_D = D - (1-a)

=== Trusted identities ===
1) mu_P - (m-2) = exact_slack_B - exact_excess_D.
2) STRONG C2 split:
   a_{m-1} b_{m-1} - a_m b_{m-2} = lc_surplus + mismatch,
   mismatch = p0*b1 - p1*b0,
   lc_surplus = b1^2 - b2*b0.
3) rise-compensation:
   rise - neg = p1*db + b1*dq,
   where db=b1-b0, dq=q1-q0, neg=-mismatch, rise=b1*(b1-b0).
4) combined decomposition:
   lc_surplus + mismatch = (rise-neg) + b0*(b1-b2).
5) hard-ratio equivalence:
   if dq<0 and db>0, then
   rise-neg >= 0  iff  (-dq/db) <= p1/b1.
6) alternative decomposition:
   combined = lc_P + lc_Q + cross + mismatch,
   lc_P = p1^2 - pm*p0,
   lc_Q = q1^2 - qm*q0,
   cross = 2*p1*q1 - pm*q0 - p0*qm,
   mismatch = p0*q1 - p1*q0.

=== Newly provided bridge equations (this was missing before) ===
A) Exact polynomial identities in canonical bridge setup:
   I(T) = (1+2x)P + (1+x)Q,
   I(A) = (1+x)P + Q  where A=T-l,
   I(B) = P + Q.

B) Therefore coefficient identities (for all k, with p_{-1}=q_{-1}=0):
   i_k(T) = p_k + 2*p_{k-1} + q_k + q_{k-1}.
   Equivalently: i_k(T) = b_k + b_{k-1} + p_{k-1}.

C) At indices m-1,m:
   i_{m-1}(T) = b1 + b0 + p0,
   i_m(T)     = b2 + b1 + p1,
   lambda = (b1+b0+p0)/(b2+b1+p1).

D) Rooted DP recurrences used to define P,Q:
   For each rooted vertex v with children c:
   dp0[v] = product_c (dp0[c] + dp1[c]),
   dp1[v] = x * product_c dp0[c].
   At root u:
   P = dp0[u], Q = dp1[u].

=== Dead ends (forbidden) ===
- Do not use mode(P') >= m-1.
- Do not use TP2 multiplicative closure for >=3 factors.
- Do not use universal exact transfer D <= 1/(1+lambda).
- Do not claim D1:=p1*q1-pm*q0 is always nonnegative.

=== Task ===
Produce exactly 5 sections:

1) Route-1 closure attempt (primary):
   Using the new equations A-D, try to prove mu_P(lambda_m(T)) >= m-2.
   Separate:
   - exact identities,
   - inequalities,
   - sign/domain assumptions.

2) If full closure fails, isolate one minimal extra lemma E:
   - state E in exact symbolic form,
   - prove that E => mu_P(lambda_m(T)) >= m-2,
   - prove all intermediate reductions fully.

3) Hard-regime closure refinement:
   Use trusted identities only to reduce (-dq/db)<=p1/b1 to the smallest missing sign lemma.
   State that sign lemma exactly (no prose placeholders).

4) STRONG C2 stitching:
   Provide the exact implication chain from
   {lc_P>=0, lc_Q>=0, cross>=|mismatch|}
   to combined>=0, then to STRONG C2.
   Keep every step algebraic.

5) Final checklist:
   List each currently missing lemma with one-line status:
   - proved now,
   - reduced to E,
   - still open and why.

Output constraints:
- Use only symbols defined above.
- No renaming core quantities.
- No brute-force or computational suggestions.
- No branching route tree; give one primary route.
