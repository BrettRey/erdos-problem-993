# Source hook: Schweitzer (2026), ultra log-concavity of matroid intersection
<!-- SUMMARY: arXiv note proving Mason-type ULC for two-matroid intersections and killing it at three partition matroids; maps the sharp boundary of the Lorentzian lane relative to graph independence systems · status: pre-use checks run, cited in main_v2 intro · updated: 2026-08-26 -->

**Source.** Schweitzer, Adam (2026). "A note on the ultra log-concavity of
matroid intersection." arXiv:2608.23262 [math.CO], v1, 24 Aug 2026. KTH
Stockholm. PDF in `~/pdf-inbox/schweitzer_2026_matroid_intersection_ulc.pdf`
(harvested by the arXiv watch, 2026-08-25 digest).

**What it proves.** Corollary 1: for any two matroids $M_1, M_2$ on the same
$n$-element ground set, the numbers $I_k$ of common independent sets of size
$k$ form an ultra log-concave sequence. The proof is Lorentzian (a
differential operator built from an M♮-concave function preserves the
Lorentzian property; Theorem 5 then restricts down). This strengthens
Ardila-Mantilla, Cristancho, Denham, Eur, Huh & Wang, "Tree metrics and
log-concavity for matroids" (arXiv:2601.02547, Theorem 1.4) from M♮-concave to
M♮₂-concave functions. And it stops exactly there: an explicit intersection of
*three* partition matroids fails log-concavity, hence ULC.

**Why this project cares.** Graph independence systems are intersections of
partition matroids, one per edge (block $\{u,v\}$ with capacity 1, singletons
elsewhere), so a tree on $n$ vertices sits at $n-1$ partition matroids. This
paper turns "the Lorentzian/ULC machinery doesn't reach graphs" from a
collection of examples into a theorem boundary: the lane closes at two
matroids, and the known LC failures for trees (Kadrawi-Levit at $n=26$, and
the rest of the failure literature in `gpt_attack/literature.md`) sit far past
it. It also explains structurally why Bendjeddou-Hardiman's pre-Lorentzian
result (cited at `paper/main_v2.tex:96`) needs a special subclass. If the
manuscript ever argues *why* the algebraic route falls short of #993, this is
the citation for the sharp edge.

**The three-matroid witness, read as a graph.** Ground set
$\{a_1,a_2,a_3,b_1,\dots,b_m\}$; matroid $M_i$ has the single nontrivial block
$\{a_i,b_1,\dots,b_m\}$. Common independent sets: all subsets of
$\{a_1,a_2,a_3\}$ plus the singletons $\{b_j\}$. That is the independence
system of the complete split graph (the $b$'s a clique, every $b$ adjacent to
every $a$), sequence $(1, m+3, 3, 1)$, LC failing at $k=2$ for $m \ge 7$.
Note: that sequence is still unimodal, so this is an LC counterexample only,
an Alavi-style shape repackaged; the novelty is the two-matroid theorem, not
the example.

**What must be checked before it is used.**
1. It is a v1 preprint, unrefereed. Before citing, replay the derivation chain
   (Theorem 5, restrict $x=0$, Proposition 1) or wait for a refereed version.
   The acknowledgements say the question came from Georg Loho, with Bérczi,
   Csikvári, Schröter and Brändén thanked, so a fuller paper from that circle
   may supersede this note. Check forward citations.
2. Whether the two-matroid case has any residual use here. Bipartite matchings
   are two-partition-matroid intersections, so Corollary 1 recovers ULC for
   bipartite matching sequences, but that is classical (Heilmann-Lieb
   real-rootedness) and matchings are not vertex independent sets. Probably
   nothing, but worth one thought before dismissing.
3. The companion pointer: arXiv:2601.02547 (Ardila-Mantilla et al., with Huh),
   "Tree metrics and log-concavity for matroids," is not in
   `gpt_attack/literature.md`. Despite the title it is matroid-side, but it is
   recent, Huh-adjacent, and should be skimmed once for anything that touches
   tree structure directly.

**Checks run, 2026-08-26 (same-day session, main context).**

1. *Derivation chain replayed by hand.* Lemma 1 (M♮-concavity is closed under
   complement-duality: direct definition check, sound); Lemma 2 (the RSW23
   criterion applied to $P_{\nu,q}$: $N(x^\kappa P(x^{-1}))$ computes out to
   $Z_{\nu^*,q}/|E|!$, Lorentzian by Lemma 1 + Theorem 4); Theorem 5 (the
   operator applied to $y^{|E|}Z_{\nu_1,q}(x,z)$ reproduces $R_{\nu_1,\nu_2,q}$
   exactly, the $|S_2|!/|E|!$ factors cancelling against
   $\partial_y^{|E\setminus S_2|}y^{|E|}$); Theorem 1 (restricting $x=0$ keeps
   exactly the $S_1=S_2$ terms; BH20 Ex. 2.3 gives ULC). External ingredients
   checked for statement fidelity only: AMCD+26 v2 Theorem 1.5 matches
   Schweitzer's quoted Theorem 4 (read in the full text, now in `literature/`);
   BH20 (published, Annals) and RSW23 taken as cited, not re-proved. The
   $q^\infty = 0$ convention making $I_k$ a plain count is stated in AMCD+26.
   Example 1 witness replayed in exact arithmetic:
   `scripts/verify_schweitzer_witness_20260826.py` (family structure, sequence
   $(1, m{+}3, 3, 1)$, LC failure exactly at $k=2$ iff $m \ge 7$, unimodality,
   the split-graph reading, and the per-edge partition-matroid reduction on 32
   trees $n \le 9$). All pass.
2. *Residual two-matroid use: nothing.* Bipartite matchings are
   two-partition-matroid intersections, so Cor. 1 recovers ULC for matching
   sequences, which is classical (Heilmann-Lieb) and about edges, not vertex
   independent sets. The one genuine consequence runs the other way: ULC
   implies LC for downward-closed families (no internal zeros), so the
   Kadrawi-Levit LC-failing trees have independence systems that are not
   two-matroid intersections on the vertex set. That remark is now in the
   manuscript.
3. *Companion skimmed* (arXiv:2601.02547 v2, 13 Jan 2026; PDF + md promoted to
   `literature/`). "Tree metrics" means ultrametric distance matrices: their
   Theorem 1.8 (a rank-1 upper bound refining Graham-Pollak) is the PSD base
   case for the Lorentzian characterization of M♮-concave functions (Theorem
   1.5, the converse to Eur-Huh). Other results: Mason M3 for M♮-concave
   functions (Theorem 1.4, what Schweitzer strengthens), Dowling 1980 / Zhao
   1985 coefficientwise conjectures (Theorem 1.6). Nothing touches tree
   independence sequences; no import for #993 beyond being the base Schweitzer
   builds on.
4. *Forward-check 2026-08-26 14:35 UTC:* still v1, no replacement, not
   withdrawn (arXiv API `id_list` query). The Loho/Csikvári-circle
   fuller-version risk stands; re-run this check before any submission that
   cites it.

**Cited in the manuscript 2026-08-26:** `paper/main_v2.tex`, new paragraph
after the Bendjeddou-Hardiman paragraph in the introduction, key
`schweitzer2026` in `paper/references.bib`.
