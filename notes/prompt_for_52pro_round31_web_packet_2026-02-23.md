# Round 31 (Web 5.2): Canonical Mu4 Closure From Embedded Packet Only

You are running in web ChatGPT with NO local repo access.
Do NOT claim missing files or inability to import modules.
Work only from the definitions and facts embedded below.

## Embedded definitions (authoritative for this round)

Canonical instance (tree domain):
- `T` is a tree with `d_leaf <= 1`.
- Canonical `(leaf, support, u)` is chosen by:
  - `leaf` minimizes parent degree; tie by smallest leaf index.
  - require `deg(support)=2`.
  - `u` is support's other neighbor.
- `B = T - {leaf, support}`, rooted at `u`.

Rooted DP at `u`:
- `P(x)=dp0[u]=product_c (dp0[c]+dp1[c])`
- `Q(x)=dp1[u]=x * product_c dp0[c)`
- `I(B;x)=P(x)+Q(x)`
- `p_k=[x^k]P(x)`.

Mode/fugacity from `T`:
- `I(T;x)=sum_k i_k x^k`
- `m = leftmost mode index of I(T)`
- `lambda = i_{m-1}/i_m`.

Route-1 retained sufficient inequality:
- `TailDef := sum_{k=0}^{m-3} (m-2-k) p_k lambda^k`
- `R1_tail2: TailDef <= p_{m-1} lambda^(m-1) + 2 p_m lambda^m`.

Equivalent threshold form:
- `TailDef = (m-2)P(lambda) - P(lambda) S`
- `Threshold := (m-2) - (p_{m-1} lambda^(m-1) + 2 p_m lambda^m)/P(lambda)`
- so `R1_tail2` is equivalent to `S >= Threshold`.

Accepted closure chain:
- `R1_tail2 => mu_P(lambda) >= m-2 => E_route1`.

Mu4 local class used in scans:
- `C_{deg,mu,mu2,mu3,mu4}` with one-step LP-dual decrement bound.
- For one step `A_new=A*F`, with `w_i = a_i lambda^i/A(lambda)`,
  moments are factorial moments up to order 4.
- `alpha_r(F) = (sum_{j>=r+1} f_j lambda^j)/F(lambda)`, `r=0..m-3`.
- `beta_i(F,m,lambda) = sum_{r=0}^{m-3} alpha_r(F) * 1_{i<=m-3-r}`.
- step decrement `d_m(A,F)=sum_i beta_i w_i`.
- LP-dual: any polynomial
  `q(i)=c0+c1 i+c2 i(i-1)+c3 i(i-1)(i-2)+c4 i(i-1)(i-2)(i-3)`
  with `q(i)<=beta_i` for all relevant `i` yields
  `d_m(A,F) >= c0 + c1 mu1 + c2 mu2 + c3 mu3 + c4 mu4`.
- `B_max` is the subset-DP optimum over child-factor orders of summed one-step lower bounds.

## Embedded computational facts (take as given)
- Canonical frontier (`d_leaf<=1`, `n<=23`): `931,596` checked, `0` failures of `B_max >= Threshold`.
- Minimum observed gap: `B_max-Threshold = 0.0002818864`.
- Collisions on key `(deg(P),m,lambda,mu1..mu4)` exist (`n<=20`: `5,893` keys), with observed threshold splits `0`.
- Through `n<=21`, observed same-mu4-key but different-mu5 count is `0`.

## Task (single decisive output)
Produce exactly one of:

A) `SUCCESS`:
- Provide one constructive canonical symbolic lemma that is sufficient to force `B_max >= Threshold` (or directly `R1_tail2`) for all canonical instances under the above packet.
- Give explicit dependency chain to `E_route1`.

B) `BLOCKED`:
- Provide a definition-complete canonical no-go argument from this packet (not unrestricted perturbation arguments outside canonical DP).
- Must identify one minimal next canonical invariant/class and a one-step lemma template.

If neither can be justified from this packet, output `INSUFFICIENT` and name the single missing canonical invariant.

## Hard constraints
- No geometric-series tails (`1/(1-lambda)` pattern).
- No truncation assumptions (`p_k=0 for k>m`).
- No claims requiring external file access.
- No fabricated witness strings or runtime outputs.

## Required output format
1) `Canonical class inputs used`
2) `Main lemma or obstruction`
3) `Derivation`
4) `Binary verdict: SUCCESS or BLOCKED or INSUFFICIENT`
5) `If BLOCKED/INSUFFICIENT: minimal next canonical invariant or class`
