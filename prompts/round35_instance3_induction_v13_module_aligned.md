# Round 35 (Instance 3): Induction v13 Aligned to ENV/A2_LOCAL/HARDSTEP Modules

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- k=1 micro-lemma branch remains locked and independent.
- v12 hard-step interface is class-indexed and regime-switched.
- Orchestrator now emits three module types: `ENV`, `A2_LOCAL`, `HARDSTEP`.

## Task

Upgrade v12 into v13 so theorem obligations match orchestrator modules exactly:

1. State v13 theorem with one unified k>=2 interface that accepts per-bucket module obligations (`ENV/A2_LOCAL/HARDSTEP`).
2. Define the minimal hard-step class template needed to represent both locked conflict geometries (`n=23` internal `a=2`, `n=24` cross-line `a=2` vs `a=4`).
3. Separate immediate proof obligations (must close now) from deferred bridge obligations (can remain as one lemma).
4. Give a non-circular dependency graph where each module contributes only forward assumptions.
5. Provide deterministic escalation rules when a bucket/module becomes internally inconsistent at higher n.

## Output format

1. `Induction theorem v13 (module-aligned)`
2. `Minimal hard-step class template`
3. `Immediate vs deferred obligations`
4. `Non-circular dependency graph`
5. `Deterministic escalation protocol`
