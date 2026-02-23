You do not have file access. Treat this as complete ground truth.
Stay strictly in packet symbols.

Project context (Erdos 993 closure):
- Hard side route is already reduced via local lemmas (`H1_local`, `H2_local`).
- This round is Route-1 only.

Target:
- E_route1: `exact_excess_D <= exact_slack_B`.

Trusted:
1) `mu_P-(m-2)=exact_slack_B-exact_excess_D`.
2) `I(T)=(1+2x)P+(1+x)Q`, `I(B)=P+Q`.
3) `i_{m-1}(T)=b1+b0+p0`, `i_m(T)=b2+b1+p1`.
4) `lambda=(b1+b0+p0)/(b2+b1+p1)`.

Known status:
- Strict local-only replacement in symbols `{p0,p1,pm,q0,q1,qm,b0,b1,b2,lambda,m}` is underdetermined.
- Rejected global-tail pattern with `1/(1-lambda)` is forbidden.

This round allows ONE explicit tail statistic:
- `TailDef := sum_{k=0}^{m-3} (m-2-k) * p_k * lambda^k`.

You may use the exact identity:
`lambda*P'(lambda) - (m-2)*P(lambda) = sum_k (k-(m-2)) p_k lambda^k`.

Task (exactly 4 sections):

1) Propose one non-circular sufficient lemma `R1_tail` using
   `{p0,p1,pm,q0,q1,qm,b0,b1,b2,lambda,m,TailDef}` only,
   such that `R1_tail => E_route1`.

2) Give full algebraic proof chain:
   `R1_tail => mu_P(lambda)>=m-2 => E_route1`.
   No skipped steps.

3) Primitive obligations (max 4):
   list local one-step DP inequalities needed to bound/propagate `TailDef` and prove `R1_tail`.

4) Binary verdict:
   - SUCCESS with explicit `R1_tail` and closure chain, or
   - BLOCKED with precise reason even with `TailDef`.

Constraints:
- No circular replacement equivalent to target.
- No TP2 closure for >=3 factors.
- No mode(P') route.
- No exact-transfer `D<=1/(1+lambda)`.
- No global geometric-tail bound `const/(1-lambda)`.
