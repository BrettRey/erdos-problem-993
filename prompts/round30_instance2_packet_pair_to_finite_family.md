# Round 30 (Instance 2): Packet-Pair Bottleneck -> Finite Inequality Family

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Exact lower-face object:
  - `G_k = Lambda_k - D + sum_all`
- Target lower face:
  - `G_k >= alpha_* R_shift`, with `alpha_* = 0.21034113597068071`.
- Proposed mechanism from Round 29:
  - pair diagonals `(2t+1, 2t+2)`,
  - define packet
    `Packet_t = Lambda_old(t)[(10)+(01)+(11)]`,
  - prove local inequality `LHS_t >= alpha_* Packet_t`,
  - summation yields global lower face.

## Task

Turn the packet-pair mechanism into a finite, checkable inequality family:

1. Write the local pair inequality in fully explicit coefficient form.
2. Isolate exactly which terms depend on tree-realizable structure (vs pure algebra).
3. Reduce the local inequality to a finite template family for `a in {2,3,4}` in boundary windows.
4. Give one sharp obstruction mode that would force `alpha_*` lower.
5. Provide a minimal computational certification plan for this local family.

## Output format

1. `Explicit packet-pair inequality (fully expanded)`
2. `Tree-realizable vs purely algebraic dependencies`
3. `Finite template reduction (a=2,3,4)`
4. `Obstruction mode and alpha-lowering mechanism`
5. `Minimal certification plan`

