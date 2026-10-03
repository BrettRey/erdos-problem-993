# Re-verification after the §4.1 switch and the formalization mentions
<!-- SUMMARY: pass re-checks after §4.1 moved to the monotonicity route and the paper began citing its Lean formalizations · status: independent verifier found no math error; fixes applied · updated: 2026-10-03 -->

Manuscript: `paper/poisson_binomial/variance-local-log-concavity-poisson-binomial.tex`.
After: sha256 `2bfed9dd59fc5032…`. PDF `5f3d6e60e1f0fb35…`, 13 pages.
Supplement zip `0393dd2373a47558…`, 44 members.

## What changed

- **§4.1** now uses the monotonicity route: A_K(δ) is nonincreasing in δ, the
  target is decreasing, and the range reduces to 32 endpoint checks with
  rounded-down margins printed. The P(H) and P_m Bernstein certificate left
  the paper.
- **§4 preamble:** b_r are positive; A is a nonnegative-coefficient quadratic
  in the weights; Q(H); a roadmap.
- **§4.2** starts with the symmetrization (R_r ≥ b_r, A ≥ ST) and the
  Bernstein conversion, both moved from the preamble and the old §4.1.
- **Mentions of the formalizations:** the abstract's last sentence, the
  outline, the supplement description and the acknowledgement (Aristotle:
  both formalizations, plus the monotonicity argument).
- **Files:** a new program, `scripts/verify_pb_compact_monotone.py`. The
  supplement builder now packages the Lean project for Proposition 3.1.
  CERTIFICATE.md, the portal fields and the cover letter are updated.

## Checks

| Pass | Evidence | Result |
|---|---|---|
| numbers-audit | `verify_pb_compact_monotone.py` asserts the 12 printed ratios, the K = 3 minimum 1.0078 on the last piece, and ST < Q at H = 7/2. The verifier recomputed m = 4, 5, 6, 15 and the K = 3 minimum by hand. The Lean kernel checks all 32 inequalities. | Correct. |
| source-reread, external-review-triage | The verifier checked the acknowledgement's two attributions. The ChatGPT Pro reconstruction is of the Bernstein certificates and Example 1.6, per the referee report, l. 5 and 487. "Supplied the monotonicity argument" is supported, because the packet asked only for monotonicity in the weights. No cited claim changed. The ChatGPT Pro review of the old §4.1 now applies to the supplement's second check; the portal-fields audit row says the new §4.1 had no cross-family review. | Supported. |
| negative-claims-audit | New negative: "The rest of Sections 2 and 3 ... is not formalized". It is true, per LOCAL_REPLAY.md and the conditional project. | No issue. |
| level-category-audit, terminological-hygiene | New terms: "monotonicity argument", "Lean proof assistant", "Mathlib", A_K(δ) (defined where used). "Proof assistant" glosses Lean on first mention in the abstract. | No issue. |
| reader-pass, coherence-cohesion, proofread | Read §4 in the rendered PDF. The verifier confirmed every `\ref` resolves and every symbol is defined before §4.2 uses it. Its fine-tier fixes were applied: b_r positive, a componentwise monotonicity statement, a forward pointer to eq:S-T for S and T, $0<\delta_-$, and the K = 3 margins note. | Fixed. |
| build-integrity | 0 undefined references, 0 overfull boxes, 13 pages. Supplement rebuilt; MANIFEST verified; all four programs pass from the extracted zip. | Clean. |
| contribution-alignment, adversarial-cold-read | Only the abstract's method sentence changed ("monotonicity argument and Bernstein expansions … verified in the Lean proof assistant"). The problem statement, theorem and contribution sentences didn't move. | Round-9 verdict stands. No new readers. |
| de-ai-prose, editorial-scar-tissue | §4.1 is written as a fresh argument, with no traces of the replaced certificate. The supplement says plainly that the Bernstein certificate is a second check. | No issue. |

**Independent verifier.** A fresh-context `responsibility-verifier` found no
mathematical error in §4, no dangling reference, and no overclaim about the
formalizations. It confirmed that `Compact.lean` proves what §4.1 describes.
Its report is saved verbatim at
`runs/pb-revision2-20261003/RETURN_verifier_sec41_monotone.md`. It runs on the
parent's model family, so it isn't cross-family evidence. All eight of its
findings were applied:
- finding 1: the supplement wording, "these programs" made specific;
- finding 2: dated provenance notes added to CERTIFICATE.md and LOCAL_REPLAY.md
  instead of editing Aristotle's files;
- finding 3: the rounding note;
- findings 4–8: the fine-tier fixes.
