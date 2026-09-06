# Independent mathematical review A

Reviewed artifact: `source/b1_zero_forbidden_core_2026-09-05.md` in this review directory. Line references below refer to that frozen file. I read the entire artifact, both frozen verifier scripts, the frozen forest-poset note, and the relevant Pascal-reserve note and Lean theorem statements/proofs. I did not read another review or prior conversation. No source or result file was modified.

## Understanding of the two results

The first result is a strict depth-three inequality for a specified class of **trees**: if \(33\le n\le38\), \(\alpha\in\{17,18,19\}\), \(2\alpha-n\le5\), and every independent \(\alpha-1\)-set extends to a maximum independent set, then \(s_3^2>s_2s_4\). The proof identifies extendable sets with all independent sets in the allowed induced forest, bounds the forbidden population, and separates a 13-case parameter calculation from one tighter structural calculation on 1,842 rooted forests.

The second result is a sufficient density criterion without the \(b_1=0\) hypothesis. The forest-poset representation makes every blocked set contain a blocked singleton or pair. An incidence count then yields \(6b_3\ge(\alpha-4)b_2\). Together with the corresponding extendable incidence count and the existing Pascal reserve, this implies strict depth-three log-concavity when

\[
\alpha\in\{17,18,19\},\qquad
\frac{b_4}{e_4}\le\frac{\alpha-7}{3(\alpha-3)}.
\]

This criterion itself does not need the window's bounds on \(n\) or \(\delta\). It applies to finite forests with the indicated independence numbers. Its use to restrict the remaining tree window is correct.

**Verdicts: both results are valid as written in their stated tree/forest context.** I found no mathematical repair needed for either conclusion. The scope and notation clarifications below matter for extracting exact standalone statements, especially for Lean. The first result remains a computer-assisted written proof, not a newly kernel-checked theorem.

## Reconstruction and checks of the first proof

### Matching bags and the allowed forest: lines 46–85

König's equality gives exactly \(\alpha\) matching bags. Every independent set occupies at most one vertex of a bag, so every maximum independent set occupies all bags. In particular, unmatched vertices are forced. In a matched bag, a forbidden endpoint has a forced mate, and a forced endpoint has a forbidden mate. This proves \(c=r+\delta\).

The potentially delicate splicing argument at lines 56–67 works. Deleting the matching pair \(xy\) leaves components consisting of whole remaining bags. The two chosen neighbors lie in different components by acyclicity. Each maximum-set restriction occupies every bag of its component, so their union has exactly \(\alpha-1\) vertices. Every other bag is occupied, and both vertices of the sole vacant bag are blocked; hence this union has no maximum extension. This really contradicts \(b_1=0\), rather than merely showing that one selected completion fails.

Once the remaining matching edges are pendant, filling unoccupied pairs at their leaf endpoints is conflict-free. Forced vertices are isolated in the allowed induced forest. Thus every independent subset of \(H\) extends to a maximum set of the original forest. Conversely, a set meeting \(U\) cannot extend. The equality of the extendable polynomial with \(I(H;x)\), used throughout the proof, is justified.

For each forbidden vertex, its flexible neighbors lie in distinct components of \(H\), and each can independently be avoided in a component maximum. Thus the argument producing an independent set of size \(\alpha+1-d_C(v)\) is sound. Excluding \(d_C(v)=0,1,2\) gives three forced neighbors. The acyclic edge bound then yields \(r\le\delta-1\).

### Coefficient bounds and the 13 parameter cases: lines 125–230

For any nonempty \(S\subseteq U\), the bipartite graph consisting of the \(S\)--\(C\) edges has at least \(3|S|\) edges and no cycles. Hence \(|N_C(S)|\ge2|S|+1\). This remains true whether or not \(S\) itself is independent; counting all subsets of \(U\) therefore gives a legitimate upper bound.

The defect shifts in (3) and (4) are correct. A forbidden intersection of size \(t\) contributes a polynomial of degree at most \(\alpha-t-1\). Thus only \(t=1\) contributes at defect two, and only \(t=1,2,3\) contribute at defect four, with exactly the coefficients displayed.

The deletion count (5) has the right direction. A lower-layer independent set leaves \(j+1\) matching pairs unoccupied and has at most \(2(j+1)\) extensions. Dividing by \(\beta_j=2^{-j}\binom mj\) gives a nondecreasing sequence \(t_j\). The cone-ray argument in (7) is valid because a nonnegative nondecreasing sequence is a nonnegative combination of terminal step sequences. Every denominator is positive in the actual parameter range \(m\ge10\).

I independently recomputed all five ray ratios and the endpoint load in all 13 cases. Exactly 12 loads are below \(7/27\), and their largest is

\[
\frac{170838}{691067}
\quad\text{at}\quad (\alpha,\delta,r)=(19,4,3).
\]

The cross-multiplied slack (224843) is correct. The remaining case is exactly \((19,5,4)\). The \(r=0\) case has no blocked sets and follows from the positive reserve.

### Four forbidden vertices: lines 258–409

The equality case in the edge count is used correctly. The \(U\cup C\) graph has 13 vertices and at least 12 forced-incidence edges, so it must be a connected tree with exactly those 12 edges. Therefore each forbidden vertex has exactly three forced neighbors and there are no \(U\)--\(U\) edges. Connectedness of the original tree then gives exactly one attachment for each flexible component.

The shared central forced vertex with two private forced leaves per forbidden vertex attains \(|N_C(S)|=2|S|+1\) for every nonempty \(S\). It therefore dominates every blocked coefficient while preserving \(E=(1+x)^9P\).

I expanded (10) independently. The difference is exactly

\[
ax(A_i-Q_i)(A_j-Q_j),
\]

whose coefficients are nonnegative. Repeated merging therefore increases the blocked polynomial and leaves the extendable polynomial unchanged. The root-change injection at lines 371–377 is also valid: images of sets containing the old base root contain its private leaf, whereas unchanged images contain neither, so there is no collision between the two branches of the map.

The corona reduction covers every flexible forest having a pendant perfect matching. Any chosen leaf endpoint has no interpair edge, so every interpair edge joins base vertices. The remaining base graph is a forest. After root improvement and concentration, rooted forests on ten base vertices are therefore sufficient.

The enumeration is exhaustive at the mathematical level: a rooted tree is a root with an unordered multiset of smaller rooted trees, and a rooted forest is an unordered multiset of rooted trees. The recursive generator implements this decomposition. Adding a new root gives the claimed second enumeration on eleven vertices. The agreement of canonical-code sets is stronger than agreement of the counts alone.

Using a separately written rooted polynomial DP, separate base-subset enumeration, and NetworkX's unrooted-tree enumeration followed by every root choice, I reconstructed all 1,842 certificate rows exactly. The maximum ratio is attained at just one rooted-forest type: the ten-vertex star rooted at its center. I obtained

\[
(e_2,e_3,e_4)=(50841,217392,652212),\qquad
(b_2,b_3,b_4)=(2051,26672,157539),
\]

and \(b_4/e_4=52513/217404\). This independent reconstruction also reproduced the 17 failures of the stronger endpoint-load test. Those failures do not challenge the argument actually used.

### The retained positive term: lines 413–455

For the original tree, deleting its attachment roots leaves independence number ten: there is at most one deleted vertex in each flexible component, and each such vertex is nonforced. The ten original matching pairs continue to bound the extension degree by two for an independent 9-set. Thus \(q_{i,1}\ge5q_{i,0}\).

The singleton-forbidden contribution gives \(b_2=\sum_iq_{i,0}\) and \(b_3\ge\sum_i(6q_{i,0}+q_{i,1})\ge11b_2\). Contributions with two forbidden vertices can only increase \(b_3\). This is correctly proved on the original tree, independently of the coefficientwise upper representative.

The extendable incidence count gives \(e_3\ge e_4/4\), and substituting these two inequalities into the exact margin identity yields (13). The constants \(7/27-v=3851/217404>0\) and \(9/2-v>0\) are correct. Strictness is justified because \(e_2e_4>0\): any maximum independent set has subsets of sizes \(\alpha-2\) and \(\alpha-4\).

## Reconstruction and checks of the second proof

At lines 524–530, the projection property gives precisely the required small obstruction. A partial maximum-assignment word either disagrees with a constant coordinate or assigns conflicting values to a comparable pair of free coordinates. If neither happens, taking an appropriate order closure extends it. The contradiction involves one or two selected graph vertices. This is a statement about subsets of the independent set, so deleting a vertex outside that witness preserves blockedness.

For a blocked set of size \(\alpha-2\), at least \(\alpha-4\) deletions preserve the fixed witness. Each resulting blocked \(\alpha-3\)-set occupies all but three matching bags and has at most six one-vertex extensions. Counting incident pairs yields (14). Singleton bags do not invalidate the bound: they reduce the number of possible extensions.

The extendable incidence count is valid for every finite graph when the relevant layer sizes exist: each extendable \(\alpha-4\)-set has at least four extensions inside a chosen maximum set, and every extendable \(\alpha-3\)-set has exactly \(\alpha-3\) deletions, all extendable.

Consequently,

\[
2e_3b_3\ge
\frac{4(\alpha-4)}{3(\alpha-3)}e_4b_2.
\]

Subtracting \(e_4b_2+b_2b_4\) leaves exactly \((h_\alpha-v)e_4b_2\), so (15) is algebraically correct. The verified differences \(c_\alpha-h_\alpha\) are \(2/63,8/405,1/108\) at independence numbers 17, 18, 19. Thus equality in the density hypothesis is allowed: the first coefficient remains strictly positive and \(e_2e_4>0\). No assumption about \(b_1\), the sign of the correction, or the sign of an individual blocked log-concavity margin is smuggled into this proof.

## Concrete changes and formalization cautions

1. **Define \(D\) locally.** At lines 242–243, 300–304, and especially (16) at line 564, the combined correction is named but never explicitly defined in this artifact. The verifier uses
   \[
   D=2e_3b_3+b_3^2-e_2b_4-e_4b_2-b_2b_4.
   \]
   Adding this definition near (13) or before its first use makes the residual-frontier statement self-contained. With this definition, discarding \(D\ge0\) from the remaining regime is correct by the positive Pascal reserve.

2. **Keep connectedness in the sharp four-forbidden lemma.** Lines 262–265 explicitly invoke it correctly; a formal statement must do so too. The bound \(b_4/e_4\le52513/217404\) is false for general forests. Take the 13-vertex shared-center forbidden core and disjointly add ten isolated edges. This forest has \(n=33,\alpha=19,\delta=5,r=4,b_1=0\), but direct exact calculation gives
   \[
   (e_2,e_3,e_4)=(94464,389376,1125504),\quad
   (b_2,b_3,b_4)=(4096,51200,289792),
   \]
   hence \(b_4/e_4=2264/8793>52513/217404\). This is a counterexample to an accidental generalization, not to the stated tree theorem.

3. **State the graph hypothesis and rank condition on (14) explicitly.** “Without assuming \(b_1=0\)” at lines 517–521 retains the forest context; it does not remove it. A convenient exact statement assumes a finite simple forest and \(\alpha\ge4\). The inequality is false on unrestricted finite simple graphs. For example, take independent vertices \(a,b,c\), three independent pairs \(X_a,X_b,X_c\), all edges between distinct \(X\)-pairs, edges joining each \(a\) to both members of \(X_a\) (and similarly for \(b,c\)), and one isolated vertex \(z\). Maximum independent sets have size five; the only blocked sets are \(\{a,b,c\}\) and \(\{a,b,c,z\}\). Thus \(b_2=1,b_3=0\), contradicting (14). I checked this example by enumerating all 1,024 vertex subsets. Again, the supplied forest proof is unaffected.

4. **Use denominator-free statements and record positivity.** Useful Lean targets are \(217404b_4\le52513e_4\), \(11b_2\le b_3\), and
   \[
   3(\alpha-3)b_4\le(\alpha-7)e_4
   \Longrightarrow s_2s_4<s_3^2
   \quad(\alpha\in\{17,18,19\}).
   \]
   In the window \(e_2,e_4>0\) is immediate from subsets of a maximum independent set, but it should be an explicit supporting lemma. Counts indexed by natural-number subtraction also need the stated lower bounds on \(\alpha\); otherwise small-rank defect indices can silently mean something different.

5. **The certificate still needs a formal completeness bridge.** The Python computation is exact and I independently replayed it, but importing 1,842 verified rows into Lean alone would not prove the graph theorem. The formal proof must connect every rooted ten-vertex forest to an enumerated canonical type, and connect the polynomial computation to independent-set counts. The current note correctly does not claim that this new result has already been Lean-verified.

These are scope, presentation, and formalization requirements, not counterexamples to either stated new conclusion.

## Evidence and limits of this review

I executed the frozen first verifier's `run(14)` in memory, redirecting only its repository lookup to the actual repository. Its counters and complete density certificate matched the frozen JSON exactly: 5,446 small trees, including 594 with \(b_1=0\); 5,438 trees eligible for the general incidence check; 105 disconnected checks; and 109 archived profile replays. No structural or coefficient failure occurred.

I executed the frozen second verifier's `run(1000, 3)` in memory. Its complete certificate rows and counters matched the frozen JSON exactly: 1,842 representatives, 2,692 small decorations, and 1,000 order-33 decorations, with no failures. Calling `run` instead of either script's `main` avoided overwriting its stored result file. Both processes used `python3 -B`.

The separate coefficient reconstruction described above used fresh polynomial multiplication and rooted include/exclude DP code, checked corona polynomials by all subsets of the ten base vertices, reconstructed actual 33-vertex trees, and checked every saved row. I also independently recomputed the 13 rational cone cases, both threshold identities, and the two scope counterexamples above.

The Pascal reserve was checked mathematically against its written Pascal-smoothing proof; the relevant existing Lean statements have the required denominator-cleared constant and rank hypothesis. I did **not** compile Lean or audit its complete axiom closure, so existing kernel-verification status is not newly certified by this review. I did not rerun the full 55-test repository suite or the adaptive search in lines 572–591. The latter is explicitly described as bounded evidence and is unnecessary for either theorem. External bibliographic and novelty claims were not verified; neither verdict rests on one.

**Final verdict:** the \(b_1=0\) tree-window closure is valid as written, with an independently replayed exact finite calculation. The low-\(b_4/e_4\) sufficient criterion is valid as written in the forest context. Neither result closes the residual regime (16), the whole depth-three window, or Erdős #993.
