You do not have file access. Treat this as complete ground truth.
Stay strictly in packet symbols.

Critical correction:
Your last “convex flux” winner is circular.
Why: assuming existence of convex f with secant `-dq/db` and right derivative `p1/b1` is effectively assuming the target secant<=tangent inequality you need to prove.
So this class of assumptions is now forbidden.

Forbidden circular forms (new):
- Existence of any potential/flux/cost/transfer function whose endpoint slope is fixed to `-dq/db` and endpoint derivative fixed to `p1/b1`.
- Any replacement lemma algebraically equivalent to `mu_P(lambda_m(T))>=m-2` or `exact_excess_D<=exact_slack_B`.
- Any convexity route is admissible only if you first define the potential explicitly as a canonical expression from packet objects (`P,Q,b,q,lambda`) and then prove its convexity from packet identities/local 2x2 inequalities; postulating convexity is forbidden.

Live goals:
A) E_hard: in hard regime (`dq<0`,`db>0`,`b1>0`), prove `p1*db + b1*dq >= 0`.
B) E_route1: prove `exact_excess_D <= exact_slack_B`.

Trusted:
1) `mu_P-(m-2)=exact_slack_B-exact_excess_D`.
2) `I(T)=(1+2x)P+(1+x)Q`, `I(B)=P+Q`.
3) `i_{m-1}(T)=b1+b0+p0`, `i_m(T)=b2+b1+p1`.
4) `lambda=(b1+b0+p0)/(b2+b1+p1)`.
5) `rise-neg = p1*db + b1*dq`.
6) if `db>0,b1>0`: `rise-neg>=0` iff `(-dq/db)<=p1/b1`.

Allowed candidate from prior round:
- If you can prove from packet structure that in hard regime
  `p1>=q1` and `db>=-2dq`,
  then `p1*db+b1*dq>=0` follows (this implication is accepted).

Task (exactly 4 sections):
1) E_hard completion (non-circular):
   - derive (or reduce to one minimal local lemma each) for:
     H1: hard regime => `p1>=q1`
     H2: hard regime => `db>=-2dq`
   - then derive E_hard.

2) E_route1 non-circular replacement R1*:
   - R1* must be pure coefficient/bridge inequalities in
     `{p0,p1,pm,q0,q1,qm,b0,b1,b2,lambda}`.
   - R1* must NOT mention `mu_P`, `exact_slack_B`, `exact_excess_D`, or `P'`.
   - prove explicitly: R1* => E_route1.

3) Minimal local lemma list:
   - list exact local inequalities needed to prove H1/H2 and R1* from DP recurrences.
   - each line must be falsifiable and checkable at `(m-2,m-1,m)`.

4) Final closure chain:
   - hard side: H1+H2 => E_hard => `(-dq/db)<=p1/b1`.
   - route1 side: R1* => E_route1 => `mu_P(lambda_m(T))>=m-2`.

Constraints:
- no brute force
- no alternative route trees
- no TP2 closure for >=3 factors
- no mode(P') claim
- no exact-transfer `D<=1/(1+lambda)`
