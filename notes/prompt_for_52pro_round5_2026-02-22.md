You do not have file access. Treat this as complete ground truth.
Stay strictly in these symbols.

Important correction before you start:
The sign lemma `(dq<0 and db>0) => (b1-b2<=0)` is FALSE on the certified frontier.
Hard-regime counterexamples from the verified artifact:
1) n=20, m=7, g6=S???????C?G?G?C?@??G??_?@??@?F~_?
   b1=3640, b2=2844, b1-b2=796, db=658, dq=-14
2) n=23, m=8, g6=V???????????_?O?C??_?A??C??C??A???_?{E??^g??
   b1=14297, b2=11314, b1-b2=2983, db=1793, dq=-48
So do NOT use `b1-b2<=0` as a route.

Goal for this round: close exactly two minimal lemmas.

Setup:
- Canonical bridge: B=T-{l,s}, deg(s)=2, u other neighbor of s.
- P=dp_B[u][0], Q=dp_B[u][1], I(B)=P+Q.
- m is leftmost mode of I(T), lambda=i_{m-1}(T)/i_m(T), a=lambda/(1+lambda).
- exact_slack_B = mu_B(lambda) - (m-1-a)
- exact_excess_D = D - (1-a), D=mu_B(lambda)-mu_P(lambda)
- db=b1-b0, dq=q1-q0.

Trusted identities:
1) mu_P - (m-2) = exact_slack_B - exact_excess_D.
2) I(T)=(1+2x)P+(1+x)Q, I(B)=P+Q.
3) i_{m-1}(T)=b1+b0+p0, i_m(T)=b2+b1+p1.
4) lambda=(b1+b0+p0)/(b2+b1+p1).
5) rise-neg = p1*db + b1*dq.
6) if db>0 and b1>0, then rise-neg>=0 iff (-dq/db)<=p1/b1.

Required outputs (only these):
A) E_route1:
   exact_excess_D <= exact_slack_B.
B) E_hard:
   in hard regime (dq<0 and db>0), prove
   rise-neg >= 0,
   equivalently p1*db + b1*dq >= 0.

Do not introduce replacement targets like `b1-b2<=0`.

Output format (exactly 4 sections):
1) Proof attempt for E_route1
   - exact identities
   - inequalities
   - sign/domain assumptions
2) Proof attempt for E_hard
   - exact identities
   - inequalities
   - sign/domain assumptions
3) If blocked: one minimal replacement lemma for E_route1 and one for E_hard
4) Final closure chain:
   - E_route1 => mu_P(lambda_m(T)) >= m-2
   - E_hard + (db>0,b1>0) => (-dq/db)<=p1/b1

Constraints:
- No brute force.
- No alternative routes.
- No TP2 closure for >=3 factors.
- No mode(P') claim.
- No exact-transfer D<=1/(1+lambda).
