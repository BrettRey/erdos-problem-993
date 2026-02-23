You do not have file access. Treat this as complete ground truth.
Stay strictly in packet symbols.

Audit result of prior round:
- Hard side reduction is accepted:
  H1_local: `p0>=q0` and H2_local: `p1-p0 >= 3*(q0-q1)` imply E_hard.
- Route-1 replacement R1* is rejected as false.

Concrete counterexample to R1* (inside canonical regime):
- n=8, g6='G?`@F_', m=3, lambda=21/23.
- p0=5, p1=8, pm=4, b1=10, b2=5.
- R1* LHS = `(m-2)p0 lambda^(m-2) + (m-1)p1 lambda^(m-1) + m pm lambda^m`
         = 328965/12167.
- R1* RHS = `(m-2)(b2+b1+p1)/(1-lambda)` = 529/2.
- LHS < RHS, so R1* fails.
- Yet Route-1 target holds on same witness: `mu_P(lambda)-(m-2)=0.769580... > 0`.

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

Task (exactly 4 sections):

1) Hard side completion:
   - Prove or reduce to primitive-step local lemmas for:
     H1_local: `p0>=q0`
     H2_local: `p1-p0 >= 3*(q0-q1)`
   under hard regime.
   - Then give exact derivation H1_local+H2_local => E_hard.

2) Route-1 non-circular replacement R1_new:
   - Must be pure coefficient/bridge inequalities in `{p0,p1,pm,q0,q1,qm,b0,b1,b2,lambda}`.
   - Must NOT mention `mu_P`, `exact_slack_B`, `exact_excess_D`, `P'`.
   - Must be strictly different from the rejected global-geometric-bound pattern used in R1*.
   - Prove explicitly `R1_new => E_route1`.

3) Primitive local obligations:
   - List minimal local inequalities (max 4 total) needed to prove H1_local, H2_local, and R1_new from DP recurrences.
   - Each obligation must be checkable at `(m-2,m-1,m)` in one primitive constructor step.

4) Final closure chain:
   - hard side: H1_local+H2_local => E_hard => `(-dq/db)<=p1/b1`.
   - route1 side: R1_new => E_route1 => `mu_P(lambda_m(T))>=m-2`.

Constraints:
- No circular potentials/convexity postulates.
- No brute force.
- No alternative route trees.
- No TP2 closure for >=3 factors.
- No mode(P') claim.
- No exact-transfer `D<=1/(1+lambda)`.
