# Round 26: Canonical Mu4 Gate (canonical theorem or canonical counterexample)

Context: Erdos #993 closure project, canonical degree-2 bridge setup.

## Established facts
- Target chain:
  - `R1_tail2 => mu_P(lambda)>=m-2 => E_route1`.
- Mu4 LP-dual class:
  - `C_{deg,mu,mu2,mu3,mu4}` with step-local `F` data and LP-dual decrement bounds.
- Empirical frontier:
  - Canonical `d_leaf<=1`, `n<=23`: checked `931,596`, failures `0`, min gap `B_max-Threshold = 0.0002818864`.
- Round-25 no-go:
  - Unrestricted moment-only perturbation gives a valid **non-canonical** no-go.
  - This is not yet a no-go for canonical DP instances.
- New canonical audit:
  - On canonical `n<=20` (`77,141` checked), key collisions on
    `(deg(P),m,lambda,mu1..mu4)` occur (`5,893` keys), but
    **no threshold splitting** observed for same key.

## Required output (binary)
Return exactly one:

A) `SUCCESS`: a canonical symbolic closure step that upgrades mu4 class from empirical to theorem-grade, or
B) `BLOCKED`: a canonical counterexample/no-go (not unrestricted perturbations).

## Strict admissibility
- You must work in the canonical class only:
  - tree DP factorization at canonical bridge root,
  - canonical `deg(s)=2`, `d_leaf<=1` constraints.
- Unrestricted distribution perturbation arguments are admissible only as side remarks; they cannot be the main verdict.
- No geometric-series tail bounds (`1/(1-lambda)`), no `p_k=0 for k>m`, no tautologies.

## What counts as SUCCESS
Provide one explicit canonical implication that is sufficient to force `B_max >= Threshold`, such as:
1. A proof that `Threshold` is a function of canonical mu4-key + fixed local invariants already controlled by the one-step LP,
or
2. A canonical invariant showing `S` is bounded below by a key-determined quantity that always exceeds `Threshold`.

Must end with full symbolic chain to `E_route1`.

## What counts as BLOCKED
Provide a **canonical** obstruction, e.g.:
1. Two canonical instances with the same relevant key/invariants but different `Threshold` or different verdict for `B_max>=Threshold`, or
2. A proof that canonical constraints still leave a free direction that preserves all mu4-class inputs while changing the target side.

If BLOCKED, give one minimal next canonical class (moment or non-moment) with one-step lemma template.

## Output format (strict)
1) `Canonical class statement`
2) `Main lemma or obstruction`
3) `Derivation`
4) `Binary verdict: SUCCESS or BLOCKED`
5) `If BLOCKED: minimal next canonical class`
