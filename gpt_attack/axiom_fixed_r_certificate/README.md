# Axiom-Style Fixed-`r` Certificate Packet

Assembled: 2026-05-27

## Objective

Give a proof agent a narrow, certificate-shaped target drawn from the fixed-`r`
Route-2 spider-lane work.  This is not a request to prove Erdős #993.  It asks
for formalizable bridge lemmas that turn existing exact certificate checks into
a theorem statement.

## Files

| File | Purpose |
|------|---------|
| `fixed_r_certificate_target.tex` | Human-readable theorem target and notation. |
| `task.md` | Instructions for a proof agent or formalization run. |
| `problem.lean` | Lean-facing skeleton for the first formal targets. |

## Context Files To Attach

Attach these repo files with this packet:

```text
gpt_attack/conjecture_and_state.md
gpt_attack/SG3_ROUTE2_PACKET.md
notes/fixed_r_proof_note_clean_2026-05-21.md
notes/fixed_r_certificate_lemma_2026-05-21.md
notes/fixed_r_global_f_mode_margin_sublemma_2026-05-21.md
notes/fixed_r_fugacity_shift_sublemma_2026-05-21.md
notes/fixed_r_huboff_reserve_recurrence_certificate_2026-05-21.md
Formal/STP2Closure.lean
```

## Success Criteria

A useful result is one of:

1. A Lean proof of the Gibbs/finite-support mean-shift lemma that replaces the
   packaged Lipschitz assumption in `problem.lean`.
2. A Lean statement and proof of the full fixed-`r` certificate criterion.
3. A sharper theorem statement that still matches the fixed-`r` certificate
   scripts and avoids hidden numerical assumptions.
4. A counterexample to one of the proposed formal bridge lemmas.

The output must distinguish exact certificate facts from floating-point
diagnostics.  Any numerical range must either be proved symbolically, checked by
exact integer/rational arithmetic, or explicitly marked as heuristic.

## Non-Goals

- Do not prove or assume global tree log-concavity.
- Do not assume the mode-mean conjecture.
- Do not assume ECMS.
- Do not restate the fixed-`r` certificate criterion as a proof for all
  `d_leaf <= 1` trees.
- Do not treat Python floating-point scans as proofs.

## Local Check

The skeleton can be syntax-checked independently of the main Lean library:

```bash
lake env lean gpt_attack/axiom_fixed_r_certificate/problem.lean
```

As of 2026-05-27, this file checks with no `sorry`s.
