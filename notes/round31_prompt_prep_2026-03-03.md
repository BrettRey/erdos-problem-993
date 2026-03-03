# Round 31 Prompt Prep (2026-03-03)

Goal: move from drift-aware monitoring to post-break repair design after n=22.

## Locked baseline at n=22

From completed mod=256 artifacts:

- `results/alpha_bookkeeping_modal_n22_n22_w256.json`
  - `min_alpha_all = 0.1875868239504603`
- `results/lambda_frontier_modal_n22_n22_w256.json`
  - `lambda_max = 0.1968360500404575`
- derived gap:
  - `alpha_front(22) - lambda_front(22) = -0.0092492260899972`

Witness classes:
- alpha-lowering witness: `(a,b)=(2,18)` at step 2
- lambda-raising witness: `(a,b)=(4,16)` at step 2

So the v8 single-envelope k>=2 route (`alpha > lambda`) is now empirically broken at n=22.

## Round 31 objectives

1. Formulate the smallest theorem-level repair compatible with all-diagonal framework.
2. Design targeted local correction around packet-pair templates for stress classes.
3. Upgrade induction to v9 with repaired k>=2 branch while preserving k=1 independence and non-circularity.
4. Convert classifier into an explicit repair planner (actions + priority + stop criteria).

## Prompt design constraints

- Keep odd-only routes excluded.
- Keep broken bookkeeping excluded.
- Use n=22 break as baseline, not as speculation.
- Demand explicit falsification/certification interfaces in every instance.

