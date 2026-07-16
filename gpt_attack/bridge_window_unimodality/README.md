# Packet: window unimodality from zero geometry (the bridge lemma)
<!-- SUMMARY: Frozen proof-agent packet connecting tree zero geometry to window coefficient monotonicity · status: RESOLVED REFUTED(T3) 2026-07-16, hypotheses proved insufficient, tree-only input required · updated: 2026-07-16 -->

## OUTCOME (2026-07-16): REFUTED(T3), independently verified

Codex (GPT-5.6) resolved the packet in refutation mode within the hour:
`P_s = ((1+x)^6 + 20x)(1+x^2)^{2s}` has positive coefficients, satisfies
every stated zero-geometry hypothesis (unique dominant real zero;
avoidance of the full degree-25 common zero-free set U_24; the 0.9
positive-axis sector; the cardioid-cusp envelope, vacuously at ratio
> 25), and has a valley at 2s-1 < 2s < 2s+1 inside the window, with
alpha = 4s+6 -> infinity and b/alpha -> 1/2. The base Q is exactly the
repo's split-graph control I(K_20 v E_6); the factor (1+x^2)^{2s}
transmits Q's odd/even imbalance (odd coefficient sum 52 > even sum 32)
to the middle of the sequence. All lemmas hand-checked and the valley,
coefficients, and exact root boxes independently replayed by the parent
session (`outcome_2026-07-16/independent_replay.json`).

Consequence: NO pointwise zero-location hypothesis class of this kind
can prove BRIDGE. The missing input is tree structure that excludes
parity-lacunary spectral factors ~-- precisely the tree-DP realizability
invariant of issue #1. The campaign's empirical law that trees never
show rise distance c-b > 1 is the shadow of this exclusion; proving it
is the reformulated bridge.

## What this packet is

The single unclaimed step between the July 2026 spectral program and
Erdos #993. Two literatures that do not cite each other:

- Zero geometry (Shearer; Scott--Sokal; Peters--Regts; Bencs--Buys--
  Peters; Bencs--Csikvari): exact descriptions of where tree/bounded-
  degree independence-polynomial zeros can and cannot be.
- Unimodality combinatorics (Levit--Mandrescu tail; Li; Galvin;
  Kadrawi--Levit; this repo's exhaustive n <= 29 and barrier laws).

Nobody has converted the first into coefficient control in the open
window k < 2 alpha / 3. That conversion is the target; see `target.md`.

## Verification rules (non-negotiable)

1. A claim is a THEOREM only with a complete, hand-checkable proof.
   Numerical evidence, however extensive, is labeled EVIDENCE.
2. Exact certificates (integer/rational arithmetic, replayable scripts)
   may support finite claims; floating-point root-finding is forbidden
   (see DECISIONS.md 2026-07-16, claw-free control failure).
3. Refutation is success: an explicit obstruction showing the stated
   hypotheses cannot imply the conclusion redirects the program and is
   a fully valid outcome (cf. the same-day collar refutation).
4. Cite precisely: theorem numbers from the papers in
   `notes/literature/` (txt conversions available), not from memory.
5. Every asymptotic claim must be effective (explicit alpha_0, explicit
   constants) or flagged as ineffective.

## Contents

- `target.md` — frozen mathematical target, graded sub-targets T1-T4
- `PROMPT.md` — dispatch prompt (Codex / GPT-5.6)
- `data.md` — certified empirical inputs from the 2026-07-15/16 campaign

## Provenance

Distilled from: notes/real_collar_conjecture_2026-07-16.md,
notes/why_trees_resist_2026-07-16.md, DECISIONS.md 2026-07-15/16,
literature intake 2026-07-16. Related packets:
gpt_attack/axiom_fixed_r_certificate/ (Route-2 Lean bridge).
