# Round 26 (Instance 2): Direct Route via sum_err Bound (No Lambda Lower Bound)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

Reliable empirical inequalities through n<=19:
- Single-constant route works:
  `sum_err <= D + 0.08144365672607116 * R_shift`.
- Defect route almost works at `delta=lambda19-lambda0`; one witness miss.
- Tiny residual correction works with explicit second-shift extra channels.

Broken route to avoid:
- Do not use `Lambda_k >= D + R_shift - sum_err` as a base theorem step.

## Task

Build a proof skeleton that treats the sum_err inequality itself as the primary bridge object.

Goal: produce a closure statement that does **not** require any false bookkeeping lower bound.

Requirements:
1. State precisely what must be shown about `sum_err` to imply STP2 closure.
2. Give one primary theorem form (single constant) and one split form (with Extra).
3. Isolate one unresolved lemma per form.
4. Keep constants explicit using locked values.

## Output format

1. `Primary closure theorem (single constant)`
2. `Split-reserve closure theorem`
3. `Implication chain to STP2`
4. `Unresolved lemma(s)`
5. `What remains to verify`
