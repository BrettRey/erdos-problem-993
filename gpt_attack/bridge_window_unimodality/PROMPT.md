# Task: prove or refute the bridge lemma (window unimodality from zero geometry)

## Outcome wanted

Work the frozen target in
/Users/brettreynolds/projects/LLM-CLI-projects/papers/queue/erdos-problem-993/gpt_attack/bridge_window_unimodality/target.md
(read it first, then README.md and data.md in the same directory, then
the literature txt files listed in data.md — all read-only; do NOT
modify anything under the repo path).

Success is exactly one of:
(a) a complete, hand-checkable proof of BRIDGE at any sub-target level
    T1/T2/T4, with explicit effective constants;
(b) a rigorous refutation T3: hypothesis-satisfying polynomial
    sequences with window valleys, alpha -> infinity, presented with
    exact verifiable instances;
(c) an honest partial: a proved reduction of BRIDGE to one or more
    sharply-stated lemmas, each with a clear statement of what is
    missing, PLUS a proof of at least one nontrivial step
    (e.g., the saddle localization x_k = O(1) for window k, or the
    oscillation damping estimate under the sector hypothesis).

Grade your own output against this bar at the end. Do not claim (a)
if you delivered (c).

## Constraints

- Verification rules in README.md are binding: theorem vs evidence
  labeling, no float root-finding, precise citations by theorem number
  from the provided txt files (not from memory), effective constants
  or an "ineffective" flag.
- The split-graph control (data.md) shows pointwise spectral bridges
  are FALSE. Your argument must use alpha -> infinity structurally.
- The high-degree hub regime is the known hard case: check every
  claimed estimate against hub-bouquet trees (a center joined to many
  spider gadgets), which realize window saddles near the critical
  activity.
- You may write and run Python for exact sanity checks (big-int only
  for coefficient claims). Put scratch work in your working directory,
  not the repo.

## Deliverable (STRICT)

Write /tmp/claude-agent-output/bridge_lemma_report.md containing:
1. Verdict line: PROVED(T_i) / REFUTED(T3) / PARTIAL(c) / FAILED.
2. The mathematics, in full, self-contained, with numbered lemmas.
3. A "verification map": for each lemma, how a human checks it
   (by hand / by which exact computation / which cited theorem).
4. What you tried that failed, in two or three sentences each.
Work for up to about 60 minutes. An honest PARTIAL with one real
proved step beats a grandiose FAILED-in-disguise.
