# Round 23 (Instance 4): Explain Extremal Pattern and Constant Jump (n<=18 -> n=19)

Use only your prior outputs plus this prompt.

## Locked context

- Through `n<=18`: sharp constant `lambda0 = 0.05201381704686925`.
- At `n=19`: required constant jumps to `0.08144365672607116`.
- New hard classes appear: `(2,15)` and `(3,14)`.
- Local neighbor compensation also starts to fail (5 slot failures).

## Task

Give a coefficient-level explanation of the constant jump and new hard classes.

Requirements:
1. Use explicit local quantities (ratios/minors/channels), not only class labels.
2. Explain why `(2,15)` and `(3,14)` are the first new break classes.
3. Relate this to defect growth (`Defect`, `d_t`, `e_t`) and reserve saturation.
4. Produce one predictive criterion for when next jump should occur at larger n.

## Output format

1. `Mechanism for lambda jump`
2. `Mechanism for new hard classes`
3. `Predictive criterion`
4. `What can be proved now vs conjectural`
5. `Minimal verification recipe`
