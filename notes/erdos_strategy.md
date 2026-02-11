# Erdos-style approach (proof-first, extremal)

Goal: avoid heavy computation; prove unimodality by contradiction on a
minimal counterexample and local structural lemmas.

## A. Minimal counterexample framework

Assume there exists a tree T with a non-unimodal independence sequence.
Choose T with minimal |V(T)|. Let k be the first "valley" index with
  i_{k-1}(T) > i_k(T) < i_{k+1}(T).

Key recurrence for a leaf v with neighbor u:
  I(T; x) = I(T - v; x) + x * I(T - N[v]; x).
Thus coefficientwise:
  i_k(T) = i_k(T - v) + i_{k-1}(T - N[v]).

Because T is minimal, both T - v and T - N[v] are unimodal. The task is to
show that their sum cannot create a new valley under suitable structural
conditions. This would contradict minimality.

## B. Structural lemmas to target (Erdos-style “force simplicity”)

1) **Branch-vertex bound.**
   Try to show a minimal counterexample cannot have 3 or more branch
   vertices (degree >= 3). Sketch: remove a leaf from a “remote” branch;
   use the recurrence above plus a monotonicity inequality on the two
   unimodal summands to show the valley cannot appear.

   If successful, the problem reduces to the class C2:
     C2 = {trees with <= 2 branch vertices}.

2) **Compression to a spine.**
   For a C2 tree, the structure is a path between the two hubs with pendant
   paths attached. Try to show that “shortening” a pendant path (or moving
   its length to the main spine) makes the sequence more unimodal (or does
   not create a new valley). In an extremal argument, the worst case should
   then be a broom or double-star.

3) **Leaf-removal monotonicity lemma.**
   Let A_k = i_k(T - v) and B_k = i_{k-1}(T - N[v]). Both are unimodal.
   Find a sufficient condition (e.g., a dominance relation on first
   differences) under which A_k + B_k is unimodal. This is the cleanest
   route to contradiction in the minimal-counterexample setup.

## C. Immediate next targets

T1) Prove a "sum of unimodals is unimodal" lemma under a checkable
    inequality (e.g., A_{k+1} - A_k >= B_k - B_{k+1} across the descent
    window). This mirrors the broom proof but is local and applies to the
    leaf-recurrence decomposition above.

T2) Establish that in a minimal counterexample, every leaf is adjacent
    to a vertex on a longest path (else remove it and use T1 to contradict).

T3) Reduce minimal counterexamples to homeomorphically irreducible trees
    by proving that edge subdivision preserves unimodality (see
    notes/subdivision_lemma.md). This is the most promising “Erdos-style”
    structural reduction found so far.

These are the kinds of "clean inequality" lemmas Erdős would push for:
find the right monotone quantity, then show any deviation can be improved
by a local move, forcing an extremal structure.

## D. Cautionary note from quick experiments

A naive “move a leaf to a diameter endpoint” operation does not appear to
monotonically increase the near-miss ratio in small random tests, so the
compression step likely needs a more subtle inequality than raw near-miss.
