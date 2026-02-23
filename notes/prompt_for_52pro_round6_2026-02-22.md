You do not have file access. Treat this as complete ground truth.
Stay strictly in packet symbols.

Your previous response produced one useful hard-regime strengthening, but the Route-1 replacement was circular.

Non-circularity rule (strict):
A replacement lemma for E_route1 is INVALID if it is algebraically equivalent to
`mu_P(lambda_m(T)) >= m-2`
or to
`exact_excess_D <= exact_slack_B`.
So do not propose `lambda P'(lambda) >= (m-2)P(lambda)` as replacement; that is equivalent to the target.

Current goals:
- E_route1: exact_excess_D <= exact_slack_B.
- E_hard: rise-neg >= 0 (equiv. p1*db + b1*dq >= 0) under hard regime.

Setup/trusted identities:
1) mu_P - (m-2) = exact_slack_B - exact_excess_D.
2) I(T)=(1+2x)P+(1+x)Q, I(B)=P+Q.
3) i_{m-1}(T)=b1+b0+p0, i_m(T)=b2+b1+p1.
4) lambda=(b1+b0+p0)/(b2+b1+p1).
5) rise-neg = p1*db + b1*dq.
6) if db>0 and b1>0: rise-neg>=0 iff (-dq/db)<=p1/b1.

Useful candidate from your prior output:
- E_hard* assumptions:
  p1>=q1 and db>=-2dq (with dq<0, db>0, nonnegative coefficients)
  imply rise-neg>=0.

Task (exactly 4 sections):
1) E_hard completion:
   - Prove or reduce to minimal lemmas for:
     (hard regime) => [p1>=q1 and db>=-2dq]
   - Then derive E_hard.

2) E_route1 non-circular replacement:
   - Give one replacement lemma R1* that is NOT algebraically equivalent to Route-1 target.
   - R1* must be stated in packet coefficients/bridge quantities (p0,p1,pm,q0,q1,qm,b0,b1,b2,lambda).
   - Prove explicitly that R1* => E_route1.

3) Minimality check:
   - Show why your R1* is strictly stronger than needed but checkable from local coefficient inequalities.
   - Identify exactly which local inequalities must be proved to get R1*.

4) Final closure chain:
   - hard side: (dq<0,db>0) + proved lemmas => E_hard => (-dq/db)<=p1/b1
   - route1 side: R1* => E_route1 => mu_P(lambda_m(T))>=m-2.

Constraints:
- No circular replacements.
- No brute force.
- No alternative route trees.
- No TP2 closure for >=3 factors.
- No mode(P') claim.
- No exact-transfer D<=1/(1+lambda).
