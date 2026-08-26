# The matroid-intersection number of a tree: a new lane, and one route killed
<!-- SUMMARY: Delta(T) is exactly the number of matroids needed to represent a tree's independence system; the first subcubic LC failure (n=34) confirms the three-matroid rung is where tree LC genuinely breaks, refuting this note's own subcubic conjecture and killing the Delta>=4 frontier · status: frame CONFIRMED, subcubic conjecture REFUTED same day, prior art UNCHECKED · updated: 2026-08-26 -->

> **Read §2.5 first.** The subcubic log-concavity conjecture stated in §2 was
> **refuted the same afternoon** by the census this note commissioned, at
> $n = 34$. The refutation *strengthens* the matroid frame rather than damaging
> it. §2 is left standing as written so the sequence of reasoning is legible.

Opened the same day the Schweitzer note (arXiv:2608.23262) was integrated into
`paper/main_v2.tex`. The manuscript sentence uses the crude count (one partition
matroid per edge, $n-1$ for a tree). That count is not tight, and making it
tight turns a passing remark into a lane.

## 1. The representation lemma

For a graph $G$ on vertex set $V$, call the **matroid-intersection number**
$\mathrm{mi}(G)$ the least $k$ such that there are matroids $M_1, \dots, M_k$ on
ground set $V$ with
$$\mathrm{IS}(G) \;=\; \mathrm{IN}(M_1) \cap \dots \cap \mathrm{IN}(M_k).$$

**Lemma.** $\mathrm{mi}(G)$ equals the least number of spanning subgraphs of $G$
covering $E(G)$, each a disjoint union of cliques. For triangle-free $G$ this is
the chromatic index $\chi'(G)$; for a tree $T$ with an edge it is $\Delta(T)$.

*Lower bound.* Every singleton is independent in $G$, so no $M_i$ has a loop. For
$uv \in E(G)$, $\{u,v\} \notin \mathrm{IS}(G)$, so some $M_i$ has $\{u,v\}$
dependent, hence (loopless) a circuit: $u, v$ are parallel in $M_i$. Parallelism
is an equivalence relation, so the dependent pairs of $M_i$ form a disjoint union
of cliques; and each such pair must be a $G$-edge, else it would be independent in
$G$ but not in the intersection. So the $k$ parallel-class partitions cover
$E(G)$ by disjoint-clique subgraphs. This uses nothing about circuits of size
$\ge 3$, so it holds for arbitrary matroids, not just partition matroids.

*Upper bound.* Given such a covering, take the partition matroid of each part
with capacity 1 on every class. A set survives all of them iff it contains no
two vertices of a common class, iff it contains no covered edge, iff it is
independent in $G$.

*Trees.* Triangle-free, so the clique subgraphs are matchings and the count is
$\chi'(T)$, which is $\Delta(T)$ by König's edge-colouring theorem for bipartite
graphs. The sharper degree argument is worth stating directly: at a vertex $v$ of
degree $\Delta$, no two neighbours of $v$ are adjacent, so a single partition can
place $v$ in a class with at most one neighbour, and covers at most one edge at
$v$. Hence $\Delta$ partitions are needed.

**Formally proved** (2026-08-26): the *lower bound*, in the general form "at any
vertex whose neighbourhood is independent, $\deg v \le k$", is machine-checked in
Lean against Mathlib and replayed locally ~-- `lake build` clean at 8,027 jobs, no
proof escape hatches, `#print axioms` reporting only `propext`,
`Classical.choice`, `Quot.sound`. With it the triangle-free and tree corollaries
($k \ge \Delta$), and a stronger form of the ground-set question than was asked:
the bound survives representation by **arbitrary minors** on a larger ground set,
so the contraction case of §6.1 is settled by proof rather than by argument.
Source `formalization/matroid_representation_lower_bound_aristotle/`; packet
`formalization/matroid_representation_lower_bound_aristotle_input_20260826.md`.
The *upper* bound is not formalized and was deliberately out of scope.

**Verified computationally** (2026-08-26, all 200 trees $n = 2..10$, max $\Delta$ seen 9): a
proper $\Delta$-edge-colouring exists, every colour class is a matching, and the
$\Delta$ induced partition matroids have common independent sets exactly
$\mathrm{IS}(T)$. The code is in DECISIONS.md's entry for the day; the crude
per-edge version is `scripts/verify_schweitzer_witness_20260826.py`.

**Prior art: NOT CHECKED.** The construction is, from memory and therefore
UNVERIFIED, what the literature calls the *equivalence covering number*
$\mathrm{eq}(G)$ (Duchet; Alon, "Covering graphs by the minimum number of
equivalence relations", c. 1986). The adversarial fork assigned to search for
this died on a model usage limit before running any query, so **no search has
been run and no novelty may be claimed for anything in this section.** Treat the
lemma as folklore until searched. The likely-classical status does not weaken
the lane: the *consequence* below is what matters, and it needs Schweitzer's
2026 theorem.

## 2. What that buys against #993

Schweitzer proves ULC at two matroids and exhibits LC failure at three. Read
through the lemma, his boundary lands on trees as a statement about maximum
degree:

| $\Delta(T)$ | matroids needed | status |
|---|---|---|
| $\le 2$ (paths) | 2 | ULC by Schweitzer Cor. 1; LC and real-rooted classically |
| $3$ (subcubic) | 3 | first rung past the theorem; Schweitzer's counterexample lives here |
| $\ge 4$ | $\ge 4$ | every known tree LC failure is here |

The last row is not a guess. Over **all 1,230 stored LC failures** (2 at $n=26$
from `results/analysis_n26.json`, 1,228 from the complete $n \le 32$ census
shards), the minimum maximum-degree is **4**, and it is attained (81 failures
have $\Delta = 4$; the distribution is $\Delta = 4$: 81, $5$: 816, $6$: 308,
$7$: 25). Two companion frontiers fall out of the same scan, both tight:

- every failure has **at least 4 branch vertices** (degree $\ge 3$);
- every failure has **at least 8 leaves** (equivalently hub weight
  $\sum_v (\deg v - 2)^+ \ge 6$, since leaves $=$ hub weight $+\,2$).

So three separate structural parameters are pinned, and the max-degree one is
exactly the parameter Schweitzer's theorem speaks to.

**Conjecture (subcubic LC).** Every tree with $\Delta \le 3$ has a log-concave
independence sequence.

This is strictly stronger than what the matroid theorem gives (which stops at
$\Delta \le 2$), and it says trees are better behaved at the three-matroid rung
than general independence systems, where Schweitzer's complete split graph fails.
Acyclicity, not the matroid count, would have to be doing that work.

**Evidence, new today.** A subcubic census (`gentreeg -D3` into
`scripts/lc_census_alpha`, driver `scripts/run_subcubic_census_20260826.py`):

- $n = 33$: **2,076,217,086 subcubic trees, zero LC failures**, total matching
  OEIS A000672 exactly. This is past the general census wall at $n = 32$.
- Gates first, all passed: the two stored $n=26$ failures reproduced as FAIL
  lines with matching polynomials (positive control on the engine); full
  subcubic runs at $n = 26$ and $n = 28$ complete, matching A000672, zero
  failures (as they must, since every stored failure at those orders has
  $\Delta \ge 4$).
- $n = 34$: **one** LC failure (see §2.5, which refutes the conjecture above).
- $n = 35$: **11,077,270,335 subcubic trees, zero LC failures.** So the $n=34$
  witness is isolated so far, and subcubic failures are very sparse rather than
  absent. $n = 36$ was still running at ship time. See
  `results/subcubic_census_20260826/` and the summary JSONs. Restartable,
  self-verifying against A000672.

A000672's offset was validated against locally generated counts at
$n = 10, 15, 20, 24$ before being used as the expected total.

## 2.5. REFUTED, same day, at $n = 34$ — and the frame comes out stronger

The census in §2 kept running. At $n = 34$ it returned **one** log-concavity
failure with $\Delta = 3$:

```
n=34  alpha=18  break at k=17 (dist 5)
par = 0,1,2,3,4,5,4,7,2,9,10,11,1,13,14,15,16,15,18,13,20,21,22,1,24,25,26,27,26,29,24,31,32,33
poly = 1,34,528,4967,31651,144723,490680,1256998,2456187,3669088,
       4172235,3571401,2256341,1018880,311700,58485,5325,72,1
```

Independently replayed in exact arithmetic from the parent array (rebuilt the
tree, recomputed the polynomial by a separate DP, confirmed the census engine's
coefficients byte for byte): it is a tree, $\Delta = 3$ exactly, 7 branch
vertices, 9 leaves, diameter 10, $\alpha = 18$. The break is at $k = 17$, where
$i_{17}^2 = 5184 < i_{16} i_{18} = 5325$, a deficit of 141. **It is still
unimodal** (peak at $k = 10$), so it is not a counterexample to #993.

Three consequences.

1. **The subcubic LC conjecture of §2 is false.** It survived 2.08 billion trees
   at $n = 33$ and died at $n = 34$. Recorded as a caution: an exhaustive check
   one rung short of the answer reads exactly like a theorem.
2. **The $\Delta \ge 4$ frontier is dead.** The minimum maximum-degree of an
   LC-failing tree is **3**, not 4. The other two frontiers from §2 survive this
   witness untouched: it has 7 branch vertices ($\ge 4$) and 9 leaves ($\ge 8$).
   So of the three parameters, the two counting *hubs* hold and the one bounding
   *degree* was an artifact of the $n \le 32$ horizon.
3. **The matroid frame is confirmed, not damaged.** This is the real result. The
   frame predicted that trees pass out of theorem-protected territory exactly at
   $\Delta = 3$, because that is where the representation needs three matroids
   and Schweitzer's guarantee stops. Before today the evidence was equivocal:
   trees at the three-matroid rung looked *better behaved* than general
   three-matroid systems, which would have meant acyclicity, not the matroid
   count, was doing the work. It is now the other way round. The tree LC-failure
   threshold and the matroid-representation threshold **coincide exactly**:

   | $\Delta(T)$ | matroids needed | LC |
   |---|---|---|
   | $\le 2$ (paths) | 2 | always holds (Schweitzer Cor. 1 gives ULC) |
   | $3$ | 3 | **can fail** — first witness $n = 34$ |

   The count is doing predictive work, which is more than a reframing of known
   facts, and it is the strongest reason to think $\mathrm{mi}$ is the right
   invariant to have named.

Also worth noting for the break-depth lane: this break has **dist 5**, whereas
every one of the 1,230 previously known failures has dist 4. It does not help
the dist-3 hunt (it is farther from the threshold, not closer), but it is the
first observed break at that depth and the first evidence that the dist
distribution has a tail.

## 3. The root-sector route, and its death

The obvious algebraic attack on the subcubic conjecture:

**Sector lemma** (elementary). If a polynomial with positive coefficients has all
roots in the sector $|\arg(-z)| \le \pi/3$, its coefficient sequence is
log-concave. (Real-factor $r + x$ and complex-pair factors
$r^2 + 2r\cos\phi\,x + x^2$ with $\phi \le \pi/3$ are each LC without internal
zeros, and LC sequences are closed under convolution.)

If subcubic trees obeyed the $\pi/3$ sector, the conjecture would follow at once.

**They do not.** `scripts/sector_angle_probe_20260826.py`, all roots certified by
Arb via python-flint (never `numpy.roots`, per the project's mandatory rule; an
angle verdict is recorded only when the ball for $\mathrm{dev} - \pi/3$ excludes
zero, and there were no undecided balls):

- *Positive control:* all 21 stored LC-failing trees have certified deviation
  $> \pi/3$, as the lemma's contrapositive requires.
- *But:* exhaustive subcubic sweeps at $n = 14, 16, 18$ and a 1-in-8 sample at
  $n = 20$ (20,512 trees total) contain **269 trees with certified deviation
  $> \pi/3$ and not one LC failure**. Max deviation reaches 1.282 rad against
  $\pi/3 = 1.047$.
- *And it gets worse with size:* complete binary trees run 1.485 ($n=31$), 1.572
  ($n=63$), 1.815 rad ($n=127$), all log-concave.

So the sector condition is violated freely by subcubic trees that are perfectly
log-concave, and the violation grows with $n$. **The route is dead:** no bound of
this shape can prove the subcubic conjecture, because the hypothesis fails on the
very objects the conclusion holds for.

*Retrospective, after §2.5:* this was the first warning and I misread it. Reading
the probe as "the sector route fails" put the fault in the method; the same data
says subcubic trees are not root-geometrically special, which was evidence
against the conjecture itself, hours before the census refuted it. When a
sufficient condition fails on an entire class one is conjecturing good behaviour
for, the class is the thing to doubt.

The probe also maps which subcubic trees are inside the sector, which is the
one reusable thing it produced: binary caterpillars come out at deviation
$0$ (real-rooted) at every size tested to $n = 90$, and binary spiders sit at
$\approx 0.30$–$0.33$ rad, well inside. Path-like subcubic trees are inside the
sector; bushy ones are far outside; both are LC. Any successful route has to be
insensitive to that split, which is a real constraint on the next attempt.

## 4. Where this sits in the queue

Ranked by what today's evidence supports, not by appeal:

1. **The leaf/hub-weight frontier** looks like the strongest of the three, because
   it is a *finite* statement per topological type: a tree with $\le 7$ leaves is
   a subdivision of one of finitely many topological trees, and the campaign
   already owns subdivision machinery (the ECMS lane and the
   $I(T_e) = I(T) + x\,I(T/e)$ identity). The proved two-hub theorem covers the
   $\le 2$-branch-vertex types already. Conjecture: **every tree with at most 7
   leaves is LC** (tight: 8-leaf failures exist).
2. ~~**Subcubic LC**~~ — **refuted at $n = 34$ (§2.5)**. Dead as a conjecture.
   What replaces it is a sharper and more interesting question: the $n=34$
   witness is a *single* failure among 4.79 billion subcubic trees, against 922
   failures among the far smaller set of all trees at $n = 32$. So subcubic trees
   are not immune, merely extremely sparse in failures. **How does the subcubic
   failure density grow?** If it stays near-zero while general density climbs,
   the interesting object is the ratio, not the boolean.
3. **Three-hub LC** ($\le 3$ branch vertices), the next rung above the proved
   two-hub theorem. Not advanced today: the fork assigned to it died on a usage
   limit before producing evidence. Its interim summary line claimed the C2
   single-index laws replicate at three hubs, and **that claim has no artifact
   behind it and must be re-derived from scratch, not cited.**

None of the three is proved, and each of the three bounds ($\Delta \ge 4$,
$\ge 4$ branch vertices, $\ge 8$ leaves) is finite-corpus evidence over
$n \le 32$ plus the subcubic extension, not a theorem.

## 5. Manuscript status

No change made to `paper/main_v2.tex` from this lane. The sentence added earlier
today uses the safe crude count ($n-1$ partition matroids, one per edge), which
is correct and verified. Sharpening it to $\Delta(T)$ should wait until the
prior-art search of §1 has actually been run.

## 6. Independent attack (2026-08-26)

Packet `gpt_attack/matroid_representation_lemma_2026-08-26.md` was run through
`ocx`, against the Z.AI/Zhipu GLM-series preview initially released as "Ox
Alpha". (Z.AI Co. confirmed the provenance on 2026-08-26, reported by Bloomberg;
until an official release name exists, that is the form for durable prose. The
labels "Ox Alpha" and `opencode/x-preview-f-free` are left untouched inside the
packet and the route manifest, which are historical artifacts.) Its verdicts,
with my assessment:

- **T1 `HOLDS-AS-STATED`.** It confirmed the higher-circuit point precisely: in a
  loopless matroid the *only* way $\{u,v\}$ is dependent is by being a circuit,
  so circuits of size $\ge 3$ cannot forbid a pair that parallelism does not
  already forbid. Degenerate cases checked and clean.
- **One real contribution, verified here.** On the ground-set worry I flagged as
  my weakest point, it gave a clean reduction: if the matroids live on
  $E \supseteq V$ and $\mathrm{IS}(G) = \{S \subseteq V : S \in \bigcap
  \mathrm{IN}(M_i)\}$, then $\mathrm{IS}(G) = \bigcap \mathrm{IN}(M_i|_V)$, and
  each restriction is a matroid on $V$. So auxiliary elements plus restriction
  buy nothing and the bound survives. I checked this concretely on partition
  matroids with auxiliary elements before accepting it. It left the
  contraction/projection case open; I expect that to close the same way, since a
  minor of a matroid is a matroid and $((M_i)/C_i)|_V$ is again a matroid on $V$.
  Both parts are now G2 in the Aristotle packet.
- **T2 followed the no-citation rule** and returned concept names with confidence
  labels plus 12 search strings, zero fabricated references. Its best suggestion
  is to search the citing literature of arXiv:2608.23262 directly, on the theory
  that whoever cites Schweitzer will already know this invariant. **The prior-art
  search itself still has not been run.**
- **T3 `PARTIAL`, and overtaken by events.** It independently identified the
  right fork ~-- does a subcubic *triangle-free* system fail LC? ~-- and proposed
  enumerating subcubic bipartite graphs to $n \approx 16$. The census answered it
  while the model was writing: yes, at $n = 34$. Its suggested search range was
  short by 18 vertices, which is a fair measure of how much the answer depended
  on having the engine rather than the idea.

### 6.1 Second round, checked (2026-08-26)

Two further unprompted returns from the same route. Brute-force `eqcov` on small
graphs (exact clique-partition enumeration) settles them:

**Kept, as a theorem with a different proof.** It concluded that the contraction
gap closes, i.e. $\mathrm{mi}^*(G) = \mathrm{mi}(G)$ even when auxiliary ground-set
elements and contractions are allowed. **The conclusion is right and the offered
proof is not.** Its argument reasons about pair-dependences surviving contraction
and then asserts that replacing $M_i/F$ by $M_i|_V$ "reproduces exactly the same
pattern", which does not follow from its own premise that contraction can *add*
dependences. The correct proof is one line and was already in the Aristotle
packet: a minor of a matroid is a matroid, so $((M_i)/F)|_V$ is a matroid on $V$,
and independence of subsets of $V$ is unchanged; hence any representation by
matroids on ground sets containing $V$ collapses to $k$ matroids on $V$. Adopt
the theorem, discard the derivation. This does remove the last definitional
asterisk on §1, which was my stated weakest point.

**Kept, verified.** Three structure facts, all confirmed by brute force:
$\mathrm{mi}(G \sqcup H) = \max(\mathrm{mi}(G), \mathrm{mi}(H))$;
$\mathrm{mi}(G) = 1$ iff $G$ is a disjoint union of cliques;
$\mathrm{mi}(G) \le \chi'(G)$ for every loopless $G$, with equality when
triangle-free. The gap can be large in the presence of triangles:
$\mathrm{mi}(K_4) = 1$ against $\chi'(K_4) = 3$, so $\mathrm{mi}$ genuinely
interpolates rather than tracking $\chi'$.

**Rejected.** Its claim that the triangle-free identity covers odd cycles, "$C_3$
included", giving $\mathrm{mi}(C_{2m+1}) = \Delta + 1$: **false at $m = 1$.**
$C_3$ *is* a triangle, so a single class $\{a,b,c\}$ is a clique covering every
edge and $\mathrm{mi}(C_3) = 1$, against $\chi'(C_3) = 3$. Verified. The claim
holds from $C_5$ on ($\mathrm{mi} = 3 = \chi'$, confirmed for $C_5$ and $C_7$).
It also billed the all-triangle-free identity as strictly stronger than §1, which
misreads §1: that statement was already made for all triangle-free $G$, with the
tree case as the corollary.

**Not usable.** The offered sweep code does not parse and enumerates independent
sets exponentially. Its experimental design also imposes $\delta \ge 2$, which
excludes trees, so it could not have found the $n = 34$ witness. Its `geng`
filter (`geng -c -b -d2 -D3 n`) is however the right invocation for the residual
question it poses, which is now the *cyclic* subcubic bipartite case, the tree
case having been settled in §2.5.

**Posture worth crediting.** It flagged that it could not verify that Erdős #993
or arXiv:2608.23262 exist as described and that everything it said was
conditional on the packet's stipulations. That is the correct stance for a model
with no retrieval, and it is the reason its citation-free output was usable at
all.
