You do not have file access. Treat this as complete ground truth.
Stay strictly in packet symbols.

Mode: Constructor-theory framing only.
Interpret statements as tasks that are possible or impossible.

Live bad tasks to rule out:
- T_hard_bad: `dq<0 and db>0 and p1*db + b1*dq < 0`.
- T_route_bad: `exact_slack_B < exact_excess_D`.

Live goals:
- E_hard: in hard regime (`dq<0`,`db>0`,`b1>0`), prove `p1*db + b1*dq >= 0`.
- E_route1: prove `exact_excess_D <= exact_slack_B`.

Trusted:
1) `mu_P-(m-2)=exact_slack_B-exact_excess_D`.
2) `I(T)=(1+2x)P+(1+x)Q`, `I(B)=P+Q`.
3) `i_{m-1}(T)=b1+b0+p0`, `i_m(T)=b2+b1+p1`.
4) `lambda=(b1+b0+p0)/(b2+b1+p1)`.
5) `rise-neg = p1*db + b1*dq`.
6) if `db>0,b1>0`: `rise-neg>=0` iff `(-dq/db)<=p1/b1`.

Non-circularity constraints:
- No postulated convex potentials/flux/cost/transfer functions unless explicitly constructed from packet objects and proved convex from packet identities.
- No replacement lemma equivalent to the target itself.
- No TP2 closure for >=3 factors, no mode(P') route, no exact-transfer `D<=1/(1+lambda)`.

Task (exactly 5 sections):

1) Primitive constructors:
   Define the minimal set of primitive bridge/DP transformations (tasks) you need.
   Each primitive must be written as an exact algebraic rewrite or inequality-preserving map on packet symbols.

2) Monotones / invariants:
   Propose two explicit monotones:
   - M_hard (must imply impossibility of T_hard_bad when nonnegative)
   - M_route (must imply impossibility of T_route_bad when nonnegative)
   They must be explicit formulas in packet symbols.

3) Preservation lemmas:
   For each primitive constructor, state exactly what must be proved to show:
   - M_hard does not decrease below 0
   - M_route does not decrease below 0
   Give each lemma as a falsifiable local inequality at `(m-2,m-1,m)`.

4) Closure theorem statements:
   Give exact theorem statements (not prose):
   - (all primitives preserve M_hard>=0) => impossible(T_hard_bad) => E_hard
   - (all primitives preserve M_route>=0) => impossible(T_route_bad) => E_route1
   Then show final implications:
   - E_hard + (`db>0`,`b1>0`) => `(-dq/db)<=p1/b1`
   - E_route1 => `mu_P(lambda_m(T))>=m-2` via trusted identity.

5) Minimal open obligations:
   Provide the smallest set (max 4) of unproved local inequalities that would complete the constructor route.
   Each obligation must be in packet symbols only and checkable on one primitive step.

Output rules:
- No metaphors.
- No new undefined objects.
- No computation suggestions.
- Must be implementation-ready for Lean-style lemma encoding.
