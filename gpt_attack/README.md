# GPT-5.5 attack bundle — Erdős Problem #993

**Author:** Brett Reynolds
**Assembled:** 2026-04-24
**Target:** A fresh AI-assisted attempt at Erdős Problem #993 (tree independence polynomial unimodality) using GPT-5.5 or any newer frontier model that can handle multi-file context and produce Lean-verifiable proofs.

## What is this bundle?

Six files that together give a frontier model everything it needs to work on the conjecture without repeating prior human/AI effort:

| File | Purpose |
|------|---------|
| `PROMPT.md` | The master prompt to paste into the chat UI (or use as a system prompt via API). Self-contained; the other files are attachments for deeper context. |
| `conjecture_and_state.md` | Precise statement of the conjecture, definitions, and an inventory of what is currently **proved**, **conjectured**, and **empirically verified** in Reynolds (2026). |
| `subgoals.md` | Five specific targets any of which would either close the conjecture or substantially advance it. |
| `literature.md` | Annotated bibliography of the relevant recent work, with the minimal citation chain needed to engage the frontier. |
| `main_v2.pdf` | Reference copy of the current manuscript. Source of truth when the summaries here gloss over details. |
| `README.md` | This file. |

## How to use

1. Paste `PROMPT.md` into GPT-5.5 (or equivalent) as the opening turn, or load it as the system prompt via API.
2. Attach the other five files in the same turn (PDF, markdown all-in-one if the UI allows multiple files).
3. Let the model pick one of the sub-goals in `subgoals.md` and propose an attack.
4. **Verify everything.** The model will confabulate. Every lemma it claims to prove must be independently checked, ideally in Lean against the existing `Formal/` project. Every citation it adds must be verified against the bibliography or arXiv.
5. If the model produces a proof sketch, feed it to Aristotle (harmonic.fun) for Lean verification against the existing `Formal/Basic.lean` and `Formal/P3.lean` infrastructure. Aristotle handled the auxiliary lemmas for the manuscript; the same pipeline can verify any new sub-lemmas.

Before launching a long run, execute:

```bash
python3 gpt_attack/sanity_check.py
```

This catches false global mean/mode targets using exact star calculations and records the current density baseline for the n=28 LC failures / near-misses.

For the next focused proof attempt, start with `SG3_ROUTE2_PACKET.md`. It
targets the degree-2-support leaf bridge and Route-2 compensation inequality,
which are the sharpest currently packaged SG3 proof obligations.

For an AxiomProver-style formalization run, start with
`axiom_fixed_r_certificate/`.  That packet isolates the fixed-`r` Route-2
certificate bridge into small Lean-facing lemmas rather than asking for a broad
proof of the conjecture.

## What success looks like

Any of:
- A proof of Conjecture A in the $d_{\mathrm{leaf}}\leq 1$ regime, together with a formal reduction showing how it closes the relevant PNP mode-control route.
- A proof of ECMS (Edge Contraction Mode Stability).
- A proof of the mode-mean bridge $\mathrm{mode}(I(T))\leq \lceil\mu(T)\rceil$, or of the manuscript's tie-fugacity condition, at least for $d_{\mathrm{leaf}}\leq 1$ trees.
- A quantitative "dip threshold" inequality that rules out large classes of hypothetical counterexamples.
- A fundamentally new structural insight that changes the attack surface.

## What failure looks like (and how to spot it fast)

- **Re-derivation of known results.** The model proves something that is already in the manuscript as a proved theorem or a cited prior result. Red flag: check `conjecture_and_state.md` first.
- **Silent assumption of a conjecture.** The model "proves" Conjecture A by implicitly assuming ECMS, or vice versa. Check that every step reduces to proved results or to genuinely new lemmas.
- **Confabulated citations.** The model cites a non-existent paper or misattributes results. Always check author + year + result against `literature.md` or Google Scholar.
- **Floating-point combinatorics.** The model reasons about "approximately" when exact integer arithmetic is required. The conjecture lives on integer coefficient sequences; there is no approximation tolerance.
- **Hand-wavy asymptotics.** The conjecture is finitary (every tree on every n must be unimodal). An asymptotic argument for large n leaves all small n uncovered. Verified-to-n=29 is the current computational boundary.
- **False global mean/mode bounds.** Do not try to prove $\mu(T)<n/3$ or $\mathrm{mode}(I(T))\leq \lfloor n/3\rfloor+1$ for all trees. Stars $K_{1,d}$ disprove both statements for large d.

## Precedent

The analogous #1014 result (R(k, ℓ+1)/R(k, ℓ) → 1) was announced 2026-04-23 by OpenAI: an internal GPT-5.5 found the proof, Mehtaab Sawhney wrote it up, Boris Alexeev formalized it in Lean. The proof uses dependent random choice, a classical graph-theoretic tool, applied in a novel setting. The lesson: long-open combinatorics problems can yield to clever elementary arguments + frontier-model guidance. The technique itself is density-based and doesn't transfer to trees, but the meta-lesson — that a standard tool may never have been properly pointed at the question — motivates a fresh attempt.

## Contact

- Paper: https://doi.org/10.5281/zenodo.19100781
- Code: https://github.com/BrettRey/erdos-problem-993
- Author: brett.reynolds@humber.ca
