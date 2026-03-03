# Round 31 (Instance 2): Targeted Local Repair Family for Packet-Pair Inequalities

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

- Local bottleneck family from Round 30 is fixed:
  - packet-pair inequalities `LHS_t >= alpha_* Packet_t`
  - finite templates for `a in {2,3,4}`
- n=22 break witnesses indicate two stress modes:
  - alpha-lowering witness in class `(2,18)`
  - lambda-raising witness in class `(4,16)`

## Task

Construct a targeted local repair family that modifies only the bottleneck region:

1. Keep the packet-pair template backbone unchanged.
2. Add the smallest explicit correction term family (in `C_ab` language or equivalent local weights) needed to absorb the identified stress modes.
3. Split obligations by template class (`a=2` vs `a=4` stress signatures).
4. Show how the repaired local inequalities still sum to a global lower-face statement.
5. Provide a practical certification plan that focuses only on the two witness classes first.

## Output format

1. `Template-level failure signatures (a=2 and a=4)`
2. `Minimal local correction family`
3. `Repaired packet-pair inequalities`
4. `Summation to global repaired lower face`
5. `Targeted certification plan`

