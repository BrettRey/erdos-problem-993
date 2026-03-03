# Round 24 Prompt Prep (2026-03-02)

Purpose: validate Round 23 candidate claims before dispatching Round 24 prompts.

Conventions: support-root, boundary-correct indexing, prefix `k < mode(I_new)`, exhaustive trees `n<=19`.

Artifact:
- `results/round24_falsification_checks_n19.json`

## Key outcomes

Full X<0 corpus size: `428,434`.

1. Candidate defect bound from Round 23 output

- Tested: `Defect <= delta * R_shift` with
  - `delta = lambda19 - lambda0 = 0.02942983967920191`
- Result: **fails in 1 case**.
- Worst ratio: `Defect/(delta*R_shift) = 1.0225408414008559`.
- Witness:
  - `n=19`, `g6=R?????????????????C??w?A^_?_~?`, `root=0`, `step=2`, `k=5`, `(a,b)=(3,14)`.

2. Proposed local lemma from Round 23 output

- Tested: `d_t <= e_{t+1} + delta*Lambda_t*K_t`.
- Result: **fails in 5,417 cases**.
- Worst local ratio: `lhs/rhs = 4.470112793366624`.

3. Extra channel proposals at fixed `c=lambda0`

- Axis candidate (`C20+C02`) fails in `12` cases, max ratio `2.014977327129186`.
- Ring candidate (`C20+C02+C21+C12+C22`) fails in `4` cases, max ratio `1.2355506702667347`.

Implication: Round 24 prompts should target repaired bounds and explicit upgraded constants/channels, not repeat the failed local lemma.
