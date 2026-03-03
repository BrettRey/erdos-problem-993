# Round 27 (Instance 2): Identity-First Rebuild (No Inherited Bookkeeping)

Use only your prior outputs plus this prompt. Do not claim fresh scans.

## Locked context

We now know inherited bookkeeping assumptions were unreliable.

Required reset:
- Start from exact CB expansion only.
- Distinguish identities from inequalities at every step.
- Do not reuse `Lambda >= D + R_shift - sum_err` unless rederived with explicit remainder terms.

## Task

Produce an identity-first decomposition for

`Lambda_k - D`

that is exact and explicitly tracks all remainder pieces.

Requirements:
1. Give one final exact identity in explicit channel notation.
2. Identify the term(s) that were previously collapsed into the wrong sign bound.
3. Propose one viable inequality target that could be true globally.
4. Keep notation compatible with `C_ab`, `R_shift`, `Extra` families.

## Output format

1. `Exact decomposition identity`
2. `Sign-sensitive remainder structure`
3. `One viable inequality target`
4. `One non-viable target (and why)`
5. `Minimal check harness`
