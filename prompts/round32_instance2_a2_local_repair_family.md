# Round 32 (Instance 2): a=2 Local Repair Family (Dual Witness Mode)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Packet-pair backbone from Round 30–31 is fixed:
  - local inequalities of type `LHS_t >= alpha_* Packet_t`
  - finite boundary templates `a in {2,3,4}` with the same telescoping identity.
- Current drifted alpha calibration in this route:
  - `alpha_* = 0.21034113597068071`.
- n=23 frontier witnesses are both in `a=2` classes:
  - alpha-lowering witness `(2,19)`
  - lambda-raising witness `(2,18)`.
- So the previous “a=2 vs a=4” stress split is obsolete at this cutoff.

## Task

Construct a targeted local repair family for the packet-pair inequalities that is specialized to the new a=2 dual-witness geometry:

1. Keep the packet-pair telescoping backbone unchanged.
2. Identify two distinct failure signatures inside the `a=2` template family (one aligned with `(2,19)`, one with `(2,18)`).
3. Add the smallest explicit correction family (in local `C_ab`/weight language) that can absorb both signatures.
4. Show how the repaired local inequalities still sum to a valid global lower-face statement.
5. Give a witness-first certification plan that calibrates constants on `(2,18)` and `(2,19)` before widening.

## Output format

1. `a=2 dual failure signatures`
2. `Minimal local correction family (a=2 focused)`
3. `Repaired packet-pair inequalities`
4. `Summation to global repaired lower face`
5. `Targeted certification plan`
