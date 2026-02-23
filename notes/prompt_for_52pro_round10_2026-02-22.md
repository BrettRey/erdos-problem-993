You do not have file access. Treat this as complete ground truth.
Stay strictly in packet symbols.

Audit of prior constructor output:
- Kept: hard-side path via local lemmas
  H1_local: `p0>=q0`
  H2_local: `p1-p0 >= 3*(q0-q1)`
  and derivation to E_hard.
- Rejected obligations:
  O2: `(b1>0)->(lambda<=p1/b1)` is FALSE.
  Counterexample in hard regime: n=23, m=8, g6='V???????????_?O?C??_?A??C??C??A???_?{E??^g??'
  lambda=38233/38884, p1/b1=13273/14297, so lambda > p1/b1.

  O3: `exact_slack_B-exact_excess_D = i_m-i_{m-1}` is FALSE.
  Counterexample: n=8, g6='G?`@F_', m=3
  exact_slack_B-exact_excess_D = 0.769580..., while i_m-i_{m-1}=2.

Do not reuse O2 or O3.

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

1) Hard side:
   - Keep only H1_local and H2_local route.
   - Give primitive-step local obligations to prove H1_local and H2_local from DP recurrences.
   - Then give exact algebraic derivation to E_hard.

2) Route-1 replacement R1_new:
   - Must be pure coefficient/bridge inequalities in `{p0,p1,pm,q0,q1,qm,b0,b1,b2,lambda}`.
   - Must NOT mention `mu_P`, `exact_slack_B`, `exact_excess_D`, `P'`.
   - Must NOT use any global geometric bound of form `constant/(1-lambda)` (R1* pattern is rejected).
   - Prove explicitly `R1_new => E_route1`.

3) Minimal obligations (max 4 total):
   - list exact local inequalities needed to prove H1_local, H2_local, and R1_new from one primitive step each.

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
