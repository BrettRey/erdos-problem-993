You do not have file access. Treat this message as complete ground truth.
Stay strictly in these symbols.

Goal: prove only the two remaining input lemmas needed for full closure.

Setup (already fixed):
- Canonical bridge: B = T-{l,s}, deg(s)=2, u = other neighbor of s.
- P = dp_B[u][0], Q = dp_B[u][1], I(B)=P+Q.
- m = leftmost mode of I(T), lambda=lambda_m(T)=i_{m-1}(T)/i_m(T), a=lambda/(1+lambda).
- p0,p1,pm and q0,q1,qm and b0,b1,b2 are the (m-2,m-1,m) coefficients.
- exact_slack_B = mu_B(lambda) - (m-1-a)
- exact_excess_D = D - (1-a), where D = mu_B(lambda)-mu_P(lambda).
- db=b1-b0, dq=q1-q0.

Trusted identities:
1) mu_P - (m-2) = exact_slack_B - exact_excess_D.
2) I(T)=(1+2x)P+(1+x)Q, I(B)=P+Q.
3) i_{m-1}(T)=b1+b0+p0, i_m(T)=b2+b1+p1.
4) lambda=(b1+b0+p0)/(b2+b1+p1).
5) combined=lc_surplus+mismatch=(rise-neg)+b0*(b1-b2).
6) rise-neg = p1*db + b1*dq.
7) if db>0 and b1>0: rise-neg>=0 iff (-dq/db)<=p1/b1.

Already closed algebraically (do not redo):
- Route 1 is equivalent to E_route1:
  E_route1 := exact_excess_D <= exact_slack_B.
- Hard-ratio is closed once E_sign holds together with combined>=0:
  E_sign := b1 - b2 <= 0 (in the hard regime).

Task for this round:
Prove only these two statements (or reduce each to one even-smaller explicit lemma):

A) E_route1:
   exact_excess_D <= exact_slack_B.

B) E_sign (hard regime):
   (dq<0 and db>0) => (b1 - b2 <= 0).

Output format (exactly 4 sections):
1) E_route1 proof attempt:
   - exact identities used
   - inequalities used
   - sign/domain assumptions used
2) E_sign proof attempt:
   - exact identities used
   - inequalities used
   - sign/domain assumptions used
3) If either fails: one minimal replacement lemma for each (symbolic form only)
4) Final dependency closure statement:
   show A+B imply
   - mu_P(lambda_m(T)) >= m-2
   - and in hard regime (-dq/db) <= p1/b1

Constraints:
- No brute force, no computation suggestions.
- No alternative routes.
- No TP2-3factor closure, no mode(P') claim, no exact-transfer D<=1/(1+lambda).
