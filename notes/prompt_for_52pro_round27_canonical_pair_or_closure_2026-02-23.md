# Round 27: Canonical Pair-Or-Closure Gate (no heuristics)

Context: Erdos #993 closure project, canonical degree-2 bridge setup.

## Current status
- Mu4 LP-dual class (`C_{deg,mu,mu2,mu3,mu4}`) empirical frontier pass:
  - canonical `d_leaf<=1`, `n<=23`: checked `931,596`, failures `0`.
- Round-26 provided a canonical obstruction narrative but no explicit canonical counterexample pair.
- Additional exact audits:
  - canonical `n<=20`: collisions on `(deg(P),m,lambda,mu1..mu4)` exist, but threshold split count is `0`.
  - canonical `n<=21`: no observed cases where same mu4-key yields different `mu5`.

## Required output (binary, strict)
Return exactly one:

A) `SUCCESS`: a constructive canonical symbolic step toward closure, or
B) `BLOCKED`: an explicit canonical counterexample pair (not an unrestricted perturbation argument).

## Branch A (SUCCESS) requirements
Provide one concrete lemma that is fully canonical and non-circular:
1. Statement in canonical symbols only.
2. Proof skeleton showing how it implies `B_max >= Threshold` (or directly `R1_tail2`) for all canonical instances.
3. Complete chain to `E_route1`.
4. No geometric-series tails, no `p_k=0` truncation, no brute-force used as proof.

## Branch B (BLOCKED) requirements
Construct an explicit canonical pair `(T1,T2)` such that:
1. Both are canonical instances (`deg(s)=2`, `d_leaf<=1`).
2. They have identical class inputs (at minimum `(deg(P),m,lambda,mu1..mu4)` and any other invariants you claim mu4-class uses).
3. They differ on target side:
   - either different `Threshold`,
   - or different `B_max-Threshold` sign.
4. Provide exact values (prefer rational where possible), and g6 strings for both trees.

If you cannot produce such a pair, you must not claim canonical no-go.

## Mandatory output format
1) `Canonical class inputs used`
2) `Lemma or explicit pair`
3) `Derivation / verification`
4) `Binary verdict: SUCCESS or BLOCKED`
5) `If BLOCKED without pair: state “INSUFFICIENT”, not BLOCKED`
