You do not have file access. Treat this prompt as complete ground truth.
Stay strictly in packet symbols.

Project context (Erdos 993 proof closure):
- We are in the canonical degree-2 bridge setup for trees.
- Hard-side work is parked; this round is Route-1 only.
- Remaining target for this round: prove E_route1.

Canonical setup:
- T is a tree, I(T;x)=sum_k i_k(T)x^k.
- m is the leftmost mode index of I(T).
- Choose canonical leaf l with support s, deg(s)=2.
- Let u be the other neighbor of s.
- B := T-{l,s}.
- P := dp_B[u][0], Q := dp_B[u][1], I(B)=P+Q.

Local coefficient notation at indices (m-2,m-1,m):
- p0:=p_{m-2}, p1:=p_{m-1}, pm:=p_m.
- q0:=q_{m-2}, q1:=q_{m-1}, qm:=q_m.
- b0:=b_{m-2}, b1:=b_{m-1}, b2:=b_m.

Mode fugacity and calibrated quantities:
- lambda := lambda_m(T) = i_{m-1}(T)/i_m(T).
- a := lambda/(1+lambda).
- mu_P(lambda) is the coefficient-tilted mean for P.
- D := mu_B(lambda)-mu_P(lambda).
- exact_slack_B := mu_B(lambda) - (m-1-a).
- exact_excess_D := D - (1-a).

Target for this round:
- E_route1: exact_excess_D <= exact_slack_B.

Trusted identities:
1) mu_P-(m-2)=exact_slack_B-exact_excess_D.
2) I(T)=(1+2x)P+(1+x)Q.
3) I(B)=P+Q.
4) i_{m-1}(T)=b1+b0+p0.
5) i_m(T)=b2+b1+p1.
6) lambda=(b1+b0+p0)/(b2+b1+p1).

Equivalent restatement you may use:
- E_route1 <=> mu_P(lambda_m(T)) >= m-2, via Trusted (1).

What is already rejected (do not reuse):
- Any replacement lemma algebraically equivalent to mu_P(lambda)>=m-2 (circular).
- Any replacement that uses global geometric tail bound pattern with 1/(1-lambda) (R1* pattern).
- TP2 multiplicative closure for >=3 factors.
- mode(P') route.
- exact-transfer D<=1/(1+lambda).

Concrete failed pattern (forbid):
- R1* form:
  (m-2)p0 lambda^(m-2) + (m-1)p1 lambda^(m-1) + m pm lambda^m >= (m-2)(b2+b1+p1)/(1-lambda)
  is false in canonical regime.
  Counterexample: n=8, g6='G?`@F_', m=3, lambda=21/23,
  p0=5, p1=8, pm=4, b1=10, b2=5,
  LHS=328965/12167 < RHS=529/2.

Important restriction for R1_new:
- R1_new must be written only in symbols:
  {p0,p1,pm,q0,q1,qm,b0,b1,b2,lambda,m}
- R1_new must NOT mention:
  mu_P, exact_slack_B, exact_excess_D, P', global tails, or undefined objects.

Task (exactly 4 sections):

1) One candidate R1_new:
   - single explicit statement in allowed symbols only.

2) Proof chain:
   - explicit algebraic proof that R1_new => E_route1.
   - every implication must be shown; no “by analogy”.

3) Minimal primitive-step obligations (max 4):
   - list local inequalities needed to prove R1_new from one-child DP update steps.
   - each must be falsifiable and checkable at indices (m-2,m-1,m).

4) Binary verdict:
   - SUCCESS: complete non-circular route with explicit chain to E_route1, or
   - BLOCKED: precise reason under this symbol restriction why no such R1_new can be completed.

Output rules:
- No metaphors.
- No alternate route trees.
- No computation suggestions.
- No references to external files.
- Be implementation-ready for Lean lemma encoding.
