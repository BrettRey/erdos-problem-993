# Triage brief: arXiv papers vs Erdős Problem #993

You are triaging arXiv papers for relevance to one research project. You
are reading, not skimming abstracts: every verdict must rest on text you
read in the paper itself.

## The project

Erdős Problem #993 (Alavi, Malde, Schwenk, Erdős 1987): for every tree (or
forest) T, the independence sequence i_0, i_1, ..., i_alpha (i_k = number of
independent sets of size k) is unimodal. Status as of 2 Oct 2026:

- Fang, Lu, Nevo, Yao, Zheng (arXiv:2609.20961) proved it for all forests
  with at least N0 vertices (N0 huge, not computed), via a CLT/local-limit
  theorem for the hard-core model uniform over forests. What is open is
  moderate n: it reduces to log-concavity of i_k on the central window
  [ceil(n/4), q].
- Log-concavity of i_k fails for some trees (Kadrawi–Levit, order 26), so
  only unimodality can hold in general. Tree independence polynomials are
  not real-rooted in general.
- Live or recently live lanes in the project: central-window
  log-concavity; matching-bag / blocked-profile decompositions (depth-3
  window, cross-reserve "Turán-type" margins absorbing negative terms);
  Lorentzian / pre-Lorentzian methods (Bendjeddou–Hardiman edge-replacement
  trees; Liu–Tang 2609.37553; Schweitzer 2608.23262 proving ULC for
  intersections of two matroids and failure at three, which bounds the
  Lorentzian route since a tree's independence system is an intersection of
  n−1 partition matroids); zeros of independence polynomials (certified
  root location, zero-free regions, Lee–Yang / hard-core model); explicit
  families (spiders, caterpillars, brooms, centipedes, Fibonacci trees;
  Galvin, Bautista-Ramos multiple-break families); Mason-type and
  free-count reformulations; LP discharging certificates; Poisson-binomial
  and total-positivity (TP2, Pólya frequency) tools; Turán inequalities.

## Verdict categories

- **relevant**: bears directly on #993. The paper is about independence
  polynomials or independent-set counts of trees, forests, or graph classes
  containing trees, or about the hard-core model on trees/forests, or it
  gives a technique stated for a class of sequences or polynomials that
  provably includes tree independence sequences, or a counterexample or
  boundary result for such a technique. Also relevant: a new general tool
  for unimodality or log-concavity of sums or products of sequences with
  explicit hypotheses that tree independence sequences could plausibly meet.
- **peripheral**: method-ecosystem awareness only. Log-concavity,
  unimodality, real-rootedness, Lorentzian or Hodge-theoretic results whose
  objects are not independence sequences of trees (matroids, Ehrhart
  polynomials, Schur functions, posets, permutation statistics, Kazhdan–
  Lusztig polynomials), with no direct transfer.
- **coincidence**: shares a keyword only ("log-concave" density on R^n,
  "totally positive" in number theory, "unimodal" distributions in
  statistics, chromatic polynomials with no coefficient-shape content).

When unsure between relevant and peripheral, say **relevant?** and explain.
The parent will read those in full. Missing a relevant paper is worse than
flagging a peripheral one.

## What to do for each paper

The text file for each paper is given. For each:

1. Read at least the first three pages: title, abstract, introduction and
   main-results statements. Use the Read tool on the .txt file (the first
   ~250 lines usually cover this).
2. Grep the full text for: `independence polynomial`, `independent set`,
   `hard-core`, `tree`, `forest`, `unimodal`, `log-concav`, `Erd` (Erdős),
   `Alavi`. Report the counts.
3. If any grep hit suggests trees or independence polynomials appear beyond
   a passing mention, read those passages.

## What to return

Return your full report in your final message (do not write report files).
One block per paper, in this exact shape:

```
### <arXiv id>: <short title>
verdict: relevant | relevant? | peripheral | coincidence
objects: <what the paper's main theorem is about, in one line>
reason: <one or two sentences tying the verdict to #993>
evidence: "<a short verbatim quote from the paper, with page or line>"
greps: indep_poly=<n> indep_set=<n> hardcore=<n> tree=<n> forest=<n> unimodal=<n> logconc=<n> erdos=<n> alavi=<n>
read: <which lines/pages you read>
```

If a text file is empty or unreadable, say so for that paper and give no
verdict. Do not guess from the title.
