# Round 28: Canonical Full-Context Gate (no "INSUFFICIENT" from missing defs)

Context: Erdos #993 closure project, canonical degree-2 bridge setup.

## Why this round exists
Round 27 returned `INSUFFICIENT` due missing formal definitions. This round supplies the exact definitions used by the verifier scripts. You must reason from these definitions only.

## Canonical setup (definition-complete)
For a tree `T` (graph6 `g6`) in the canonical domain:
- Filter: `d_leaf <= 1`.
- Choose canonical leaf triplet `(leaf, support, u)` by:
  - `leaf` minimizes parent degree, tie by smallest leaf index;
  - require `deg(support)=2`;
  - `u` is support's other neighbor.
- Let `B = T - {leaf, support}` and root `B` at `u`.

Rooted DP at `u` (children `c` of `u` in rooted `B`):
- `P(x) = dp0[u] = product_c (dp0[c] + dp1[c])`
- `Q(x) = dp1[u] = x * product_c dp0[c)`
- `I(B;x)=P(x)+Q(x)`
- `f_c(x)=dp0[c]+dp1[c]`, and `P = product_c f_c`.

Mode/fugacity for `T`:
- `I(T;x)=sum_k i_k x^k`
- `m = leftmost mode index of I(T)`
- `lambda = i_{m-1}/i_m`.

## Exact Route-1 tail objects (from verifier)
Let `p_k=[x^k]P(x)`.

- `P(lambda) = sum_k p_k lambda^k`
- `TailDef = Tail_m(P) := sum_{k=0}^{m-3} (m-2-k) p_k lambda^k`
- `R1_tail2 := TailDef <= p_{m-1} lambda^(m-1) + 2 p_m lambda^m`

Equivalent threshold form (exact):
- `S` is telescoped decrement from one-step terms,
- `TailDef = (m-2)P(lambda) - P(lambda) S`, hence
- `Threshold := (m-2) - (p_{m-1} lambda^(m-1) + 2 p_m lambda^m)/P(lambda)`
- `R1_tail2` is equivalent to `S >= Threshold`.

Route-1 closure (accepted):
- `R1_tail2 => mu_P(lambda) >= m-2 => E_route1`.

## Mu4 LP-dual class used by script
Class: `C_{deg,mu,mu2,mu3,mu4}`.
At one step `A_new=A*F` with `A(x)=sum_i a_i x^i`, `F(x)=sum_j f_j x^j`, `a_i,f_j>=0`:
- `w_i := a_i lambda^i / A(lambda)`
- factorial moments `mu_r = sum_i (i)_r w_i` for `r=1..4`.

Step coefficients:
- `alpha_r(F) = (sum_{j>=r+1} f_j lambda^j)/F(lambda)`, `r=0..m-3`
- `beta_i(F,m,lambda) = sum_{r=0}^{m-3} alpha_r(F) * 1_{i<=m-3-r}`
- step decrement `d_m(A,F)=sum_i beta_i w_i`.

LP-dual step lower bound:
- choose `q(i)=c0 + c1 i + c2 i(i-1) + c3 i(i-1)(i-2) + c4 i(i-1)(i-2)(i-3)`
- constraints `q(i) <= beta_i` for all `i=0..deg(A)`
- then `d_m(A,F) >= c0 + c1 mu1 + c2 mu2 + c3 mu3 + c4 mu4`.

`B_max` is the subset-DP optimum over child-factor orders of summed step lower bounds; script checks `B_max >= Threshold`.

## Established computational facts
- Mu4 class frontier pass (canonical `d_leaf<=1`, `n<=23`):
  - checked `931,596`, failures `0`, minimum gap `B_max-Threshold = 0.0002818864`.
- Collisions on key `(deg(P),m,lambda,mu1..mu4)` exist (`n<=20`: 5,893 keys), but observed threshold split count is `0`.
- `n<=21`: observed same-mu4-key / different-mu5 count is `0`.

## Task (binary, strict)
Return exactly one:

A) `SUCCESS`:
Provide one constructive canonical symbolic lemma that is sufficient to force `B_max >= Threshold` (or directly `R1_tail2`) for all canonical instances under the above definitions.

B) `BLOCKED`:
Provide a canonical no-go with one of:
1. explicit canonical pair `(T1,T2)` with identical class inputs used by your lemma but different target side, or
2. a definition-complete proof that canonical constraints still leave a free direction preserving all mu4-class inputs while changing the target side.

If you cannot provide A or B from these exact definitions, output `INSUFFICIENT` and state the single minimal missing canonical invariant needed.

## Hard constraints
- No geometric-series tails (`1/(1-lambda)` pattern).
- No truncation assumptions (`p_k=0 for k>m`).
- No unrestricted perturbation argument unless you explicitly map it into the canonical DP construction above.
- No invented identities for `exact_slack_B` / `exact_excess_D`.

## Required output format
1) `Canonical class inputs used`
2) `Lemma or canonical obstruction`
3) `Derivation` (must reference the exact definitions above)
4) `Binary verdict: SUCCESS or BLOCKED or INSUFFICIENT`
5) `If BLOCKED/INSUFFICIENT: minimal next canonical invariant or class`
