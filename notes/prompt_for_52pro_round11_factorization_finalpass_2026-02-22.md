You do not have file access. Treat this as complete ground truth.
Stay strictly in packet symbols.

Project context (Erdos 993 proof-closure):
- We are proving two remaining lemmas in the canonical degree-2 bridge setup:
  - `E_hard`: in hard regime (`dq<0`,`db>0`,`b1>0`), prove `p1*db + b1*dq >= 0`.
  - `E_route1`: prove `exact_excess_D <= exact_slack_B`.
- This prompt addresses only the `E_hard` route via explicit 2x2 factorization.

This is the final pass for the 2x2-factorization route.
Success criterion is binary:
- either you output explicit formulas for `J_attach(lambda)` and `J_boundary` that `simp` to the required factorization,
- or you output a no-go result that rules out the plausible candidate shapes and we drop this route.

Core normal form (exact identity):
- `p1 = b1 - q1`
- `db = b1 - b0`
- `dq = q1 - q0`
- `rise-neg = p1*db + b1*dq`
- Equivalent determinant form:
  `rise-neg = det([[b1, -p1],[db, dq]])`.
Call
`M_target := [[b1, -p1],[db, dq]]`.

Required factorization target:
`M_target = J_attach(lambda) * J_boundary`.
Then
`rise-neg = det(J_attach(lambda))*det(J_boundary)`.

Non-circularity constraints (strict):
1) `J_attach(lambda)` must be gadget-only: dependence allowed only on `lambda` and fixed attachment constants from `(1+2x),(1+x)`; no direct dependence on `{b0,b1,q0,q1,p0,p1,pm,qm}` beyond lambda.
2) `J_boundary` must be boundary-only from canonical local degree-slice data (`b0,b1,q0,q1,db,dq` and linear transforms like `q-b`).
3) You may NOT define either matrix by solving backwards to encode the target determinant.
4) No postulated convex potentials; no TP2 closure for >=3 factors.

Task (exactly 5 sections):

1) Explicit matrices:
   - Give literal formulas for `J_attach(lambda)` and `J_boundary` entries.
   - No placeholders, no existential wording.

2) Exact multiplication check:
   - Show algebraically that `J_attach(lambda)*J_boundary = M_target` entrywise.
   - Must be simplifiable by direct symbolic algebra.

3) Determinant sign obligations:
   - State exactly two local obligations:
     `det(J_attach(lambda)) >= 0` and `det(J_boundary) >= 0`.
   - Each must reduce to a single local 2x2 inequality in packet symbols.

4) If explicit factorization fails:
   - Provide a no-go theorem for this candidate-shape family:
     define the candidate family clearly (finite list of row/column choices from `{(b1,q1),(b0,q0),(db,dq)}` and `q-b` transform),
     and show why no candidate can satisfy gadget-only `J_attach(lambda)`.
   - Include at least one concrete contradiction pattern (symbolic, not numeric search-only prose).

5) Binary verdict:
   - `SUCCESS` with explicit `(J_attach,J_boundary)` and obligations, or
   - `NO-GO` with the no-go theorem statement and why this route should be dropped.

Output rules:
- No metaphors.
- No alternate routes.
- No references to files or computation.
