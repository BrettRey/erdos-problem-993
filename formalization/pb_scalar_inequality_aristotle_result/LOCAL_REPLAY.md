# Local replay of the Aristotle return for Proposition 3.1
<!-- SUMMARY: local replay of Aristotle's Lean proof of the PB paper's scalar inequality (Prop 3.1) · status: replay passed; full statement kernel-checked, standard axioms only · updated: 2026-10-03 -->

**Return.**
- Aristotle (Harmonic) project `bd2ed5ea-aff8-4c72-bfcd-4c46aca8957d`, run `16660e17-944d-4715-a182-78bdea524f43`.
- Downloaded 2026-10-03 with `scripts/aristotle_cli.py result`. The archive, kept outside git, has sha256 `1e53600c360337456f181ee8fa7664a81c360f4fe0dfee5d0970e09f0e3e507e`.
- The packet sent was `formalization/pb_scalar_inequality_aristotle/` (committed in `7a3da7e`).

**Self-assessment (Aristotle's claim):** COMPLETE.

## What was replayed here (Claude Code, Opus 5.5, 2026-10-03)

1. **Build.** `lake exe cache get`, then `lake build`, on Lean 4.28.0 (commit
   `7e01a1bf5c70`) with Mathlib v4.28.0 (manifest rev `8f9d9cff6bd7…`):
   "Build completed successfully (8038 jobs)". I rebuilt after moving the
   folder to this path, with the same result.
2. **Axioms.** `lake env lean PBScalar/Axioms.lean` gives
   `[propext, Classical.choice, Quot.sound]` for `scalar_inequality`,
   `scalar_inequality_compact` and `scalar_inequality_large`, and for the two
   cross-check theorems shown there. Neither `Lean.ofReduceBool` nor
   `Lean.trustCompiler` appears.
3. **Escape hatches.** I grepped `PBScalar/` and `PBScalar.lean` for `sorry`,
   `admit`, `axiom`, `implemented_by`, `native_decide`, `ofReduceBool`,
   `trustCompiler`, `unsafe`, `@[extern`, `opaque`, `skipKernelTC`,
   `macro_rules`, `elab` and `set_option`. There were no matches except
   `set_option maxRecDepth 20000`, used twelve times in the cross-check files
   `PaperCells/Cells{A,B,C}.lean`, which is harmless. The three uses of
   `decide +kernel` in `Compact.lean` are reduction by the kernel itself.
4. **Statements.**
   - `a`, `R`, `L`, `w` and `A` in `PBScalar/Defs.lean` are identical, after
     whitespace normalization, to the packet's `PBScalar/Statement.lean`.
   - The three theorem signatures in the returned `PBScalar/Statement.lean` are
     identical to the packet's up to `:=`. Only `by sorry` was replaced, by
     `compact_main …`, `large_main …`, and the assembly proof.
   - `Statement.lean` adds no definitions.
5. **Inputs unchanged.** The certificate JSON,
   `data/verify_pb_large_h_range.py`, `PROMPT.md`, `PROOF_CONTEXT.md` and
   `lean-toolchain` are byte-identical to the packet's.
6. **The new G1 route, checked independently.** In Python with exact
   rationals, written without reference to the Lean, I checked
   $(3+d_{lo})/(4d_{lo}^2)\leq A(d_{hi},K)$:
   - for $K=3$ on 20 equal sub-intervals of $[1/5,1/4]$, all passing, tightest
     ratio $1.0079$;
   - for $K=m$ on $[1/(m+2),1/(m+1)]$, $m=4,\ldots,15$, all passing, with
     ratios from $1.072$ ($m=4$) to $3.517$ ($m=15$).

   I also read `PBScalar/Compact.lean` in full. The argument is: $A$ is
   nondecreasing in nonnegative weights; for fixed $K$ every weight is
   nonincreasing in $\delta$ while $(K+1)\delta\leq1$; the target decreases
   in $\delta$; so one endpoint check per sub-interval suffices. The cell
   $[1/(m+2),1/(m+1)]$ with $K=m$ matches the hypotheses
   $(K+1)\delta<1\leq(K+2)\delta$.

## Verdict

**Replay passed.** Proposition 3.1, the full scalar inequality for
$0<\delta<1/4$ in the formulation of the packet, is kernel-checked in Lean 4
with Mathlib and depends only on the three standard axioms.

Scope and limits:
- The Lean definitions formalize the paper's $A(\delta)$ with $K$ passed
  explicitly and characterized by $(K+1)\delta<1\leq(K+2)\delta$. This matches
  (2.7) in the paper.
- **The compact range is proved by a different route from the paper's** (monotonicity plus 32 rational checks, not the 275 Bernstein coefficients).
  The paper's own certificate is cross-checked separately: all 275 coefficients
  are positive and the identities hold. However, for $m\geq4$ that check is a
  statement about $S_mT_m$ and is not linked to $A$ in Lean.
- **Not formalized:** the deduction of Theorem 1.1 from Proposition 3.1
  (§§2–3: the Hillion–Johnson recurrence, the mass bounds, the
  maximal-mass bound), Example 1.6, and Theorem 1.1 itself.

## Note added 2026-10-03, after the manuscript adopted this route

After this replay, the manuscript's §4.1 switched to the monotonicity route
proved here in `PBScalar/Compact.lean`. The Bernstein certificate became the
supplement's independent second check of 3 < H ≤ 16.

So the comments in Aristotle's files that describe "the paper's" route for the
compact range now describe the earlier draft. That covers `README.md`,
`ARISTOTLE_SUMMARY.md`, `Compact.lean`, `Statement.lean`,
`PaperIdentities.lean` and `PaperCells/Defs.lean`, and their references to the
identity $A-Q=P(H)/(4H^5(H+1)^3)$ and the polynomials $P_m$. These files are
left unchanged, because the provenance claims above depend on that.

The ratios in item 6 are rounded to nearest: 1.0079 and 3.517. The manuscript
prints the same values rounded down, 1.0078 and 3.516. The exact values are
1.00786… and 3.51664….
