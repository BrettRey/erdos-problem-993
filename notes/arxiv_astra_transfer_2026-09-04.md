# Weekly arXiv and GPT-6 Astra transfer audit
<!-- SUMMARY: The arXiv watcher ran on 1--4 September and queued ten papers. The Astra-inspired obstruction experiment found an exact first-certificate partition, formal prefix closure, and exact induced-forest formulas for unary classes. Direct individual-class closure fails at n=6, and an n=17 counterexample closes the relative-ratio/grouped-summand route. A cross-reserve endpoint bound survives 204,995 relevant trees through n=18. A 600-tree random probe in the actual alpha=17--19, deficiency-at-most-5 rung found no adverse cases, but 500 targeted lifts of the n=17 witness produced 57. All 57 satisfy the endpoint bound. Nine lifts refute blocked-profile log-concavity, so the live sufficient bridge must absorb the mixed and blocked deficits jointly; the resulting stronger conditional bound also survives all 57 adverse lifts. None of these auxiliary results refutes Erdős #993. -->

Date: 2026-09-04. Status: source audit and recommended bounded experiments completed. No proof of Erdős #993, no all-tree claim, and no manuscript change.
## 1. arXiv watcher status
The watcher has been active this week. The loaded LaunchAgent `com.brettreynolds.arxiv-watch` is scheduled daily at 06:15, has last exit code zero, and wrote digest entries on 1, 2, 3, and 4 September. Those four entries queued 1, 3, 2, and 4 papers respectively: ten papers in total.

The 4 September run reported partial coverage because one of 22 queries failed. The digest was still populated, so this is a partial-query warning rather than a failed daily job. The stderr log also contains AppleScript syntax errors from the desktop-notification path; collection and digest writing succeeded, but a successful exit should not be read as proof that every notification displayed. The durable operational sources are `../../../Project-Management/arxiv-watch/digest.md`, `../../../Project-Management/logs/arxiv-watch/launchagent.err`, and `~/Library/LaunchAgents/com.brettreynolds.arxiv-watch.plist`.
## 2. The local target
For the matching-bag decomposition

\[
s_d=e_d+b_d,
\]

the open depth-three margin is

\[
s_3^2-s_2s_4
=(e_3^2-e_2e_4)+(b_3^2-b_2b_4)+C,
\qquad C=2e_3b_3-e_2b_4-e_4b_2.
\]

The extendable term has the proved Pascal reserve. The remaining work is to understand the blocked profile and absorb an adverse mixed term. Blocked partial assignments are precisely assignments that violate a transitive comparison in the forest-poset gauge; at defect at most four, every minimal obstruction is an alternating comparison path. The hard part is overlap: one assignment can contain several obstruction paths.
## 3. Astra solution #74: the useful proof pattern
GPT-6 Astra disproved [Erdős problem #74](https://www.erdosproblems.com/74). The benchmark paper reports a $218, 15-hour default run; the result has a [mathematical exposition](https://www.erdosproblems.com/static/74-proof.pdf) and [Lean repository](https://github.com/tadamcz/erdos74). The present audit read the exposition and source corpus but did not independently replay the Lean build.

The theorem concerns graphs of infinite chromatic number, but three proof moves are relevant here.

1. A failure is represented by a small certificate: a short odd closed walk.

2. The proof recursively selects a witness and splits according to an element that hits it, compressing a large overlapping witness family into a finite controlled core.

3. Sparse exceptional sets are localized at successive scales; a pigeonhole argument finds a clean band across which two local colorings can be glued.


This is not a theorem about matching bags. It suggests a concrete adaptation: freeze an order on missing transitive comparisons (and on their unique cover paths), and charge every blocked assignment to its first violated comparison. That turns a union of overlapping path events into disjoint classes. The analogue of Astra's clean band is the explicit Pascal reserve: the goal is to show that the exceptional classes occupy too little or are shaped too regularly to exhaust it.

The transfer is methodological but unusually well aligned with the existing first-violated-comparison candidate in `theoremgraph_attack_and_bridge_workflow_2026-08-28.md`.
## 4. Astra solution #126: exact translation and failed direct import
GPT-6 Astra proved [Erdős problem #126](https://www.erdosproblems.com/126). The benchmark paper reports a $247, 16-hour default run; see the [mathematical exposition](https://www.erdosproblems.com/static/126-proof.pdf) and [Lean repository](https://github.com/tadamcz/erdos126). Again, the present audit read the proof corpus but did not independently replay its Lean build.

Its core construction is genuinely cross-field useful. For each prime it builds a laminar family of residue supports, forms positive-semidefinite co-membership kernels, gauges them by signs, and separates cross-sign mass from same-sign mass. A conditionally negative kernel controls the cross-sign part. A two-copy Hall injection then compares diagonal mass with off-diagonal mass. The important proof-engineering lesson is to expose positive-semidefinite summands and compare them before addition.

For #993, let (X_v) indicate that vertex (v) is free over a uniformly random (k)-independent set, and put (p_v=E[X_v]). Mason's inequality (b) becomes

\[
2\sum_{\{u,v\}\notin E(T)}\operatorname{Cov}(X_u,X_v)
\le \sum_v p_v^2+2\sum_{\{u,v\}\in E(T)}p_up_v. \tag{1}
\]

Equation (1) makes Astra's diagonal-versus-off-diagonal architecture tempting, but the direct kernel route fails three exact structural tests.

- Positive-covariance supports are not laminar in general; overlapping tree paths need not be disjoint or nested.

- At (k=1) on the five-vertex star, the positive non-edge covariance kernel is not conditionally negative semidefinite. Its leaf-leaf entries are (1/25), and the zero-sum vector ((-4,1,1,1,1)) has quadratic form (12/25>0).

- With four leaves, no single sign assignment can place every positive leaf-leaf pair across the sign cut: some pair necessarily receives the same sign.


Therefore Proposition 1 of the #126 proof cannot simply be applied to the free indicator covariance matrix. What survives is its design pattern: decompose by canonical local certificates, retain positive pieces, and prove comparisons at the summand level before summing.
## 5. Initial \(n\le15\) observation and its \(n=17\) refutation
The existing exact probe was extended to record the cross-multiplied relative normalized log-concavity margin

\[
R=b_3^2e_2e_4-b_2b_4e_3^2. \tag{2}
\]

Over every one of the 13,179 relevant non-isomorphic trees through \(n=15\):

- the mixed term (C) is negative 223 times;

- all 223 adverse cases have (R\ge0);

- (R<0) occurs 49 times, and every one of those 49 cases has (C\ge0).


Thus the initial bounded corpus suggested the implication

\[
C<0\quad\Longrightarrow\quad
b_3^2e_2e_4\ge b_2b_4e_3^2. \tag{3}
\]

The evidence was much thinner than the count 223 suggested: 219 cases had \(b_2=0\), three more had \(b_4=0\), and only one had \(b_2b_4>0\) and hence a genuinely nontrivial right-hand side in (3).

An optimized extension through \(n=17\) refutes (3). Among 81,128 relevant trees it finds 812 cases with \(C<0\), 123 with \(R<0\), and one in their intersection. The first counterexample has graph6 code ``PpG`A?@O??g??@O???I??O??``, \(n=17\), \(\alpha=11\),

\[
(e_2,e_3,e_4)=(584,1408,2241),
\qquad
(b_2,b_3,b_4)=(1,8,44),
\]

and exact margins

\[
C=-5409,
\qquad
R=-3{,}468{,}800.
\]

This does not threaten the actual frontier inequalities: the extendable reserve is 673,720, the cross-reserve margin is 668,311, the blocked Turán margin is 20, and the full depth-three margin is 668,331. Equation (3) was an overly strong proposed bridge and is now closed.

Replay:

```text
python3 scripts/probe_blocked_profile_depth3_20260828.py
python3 scripts/probe_blocked_profile_depth3_20260828.py --max-n 17 \
  --output results/blocked_profile_depth3_probe_n17_20260904.json
```

The bounded counters are in `../results/blocked_profile_depth3_probe_20260828.json`; the refuting certificate is in `../results/blocked_profile_depth3_probe_n17_20260904.json`.
## 6. Weekly corpus triage
The two papers with the best immediate proof-method fit are the following.

- Guo and Kang, [_Palindromic real-rooted polynomials and h-vectors of polytopes_](https://arxiv.org/abs/2609.00781), use nested ballot complexes and construct a divisibility-closed pure multicomplex from summands. If first-violation classes are closed under deleting an erased bag, or admit a ballot completion, their rank profiles may acquire an O-sequence/down-set representation. The applicability test is exact down-closure; absent it, the polytope theorem itself does not transfer.

- Pahari, [_Log-concavity and unimodality of Hodge numbers of Hilbert schemes of points over a surface_](https://arxiv.org/abs/2609.04095), writes coefficients as shifted sums and proves the difficult cross-shift comparison by matching summands and using monotone consecutive ratios before addition. This motivated the class and group tests below, but the \(n=17\) counterexample shows that the direct relative-ratio adaptation is too strong.


One paper is best treated as a discovery tool rather than a source theorem.

- Baldi and Kummer, [_A note on bounded ratios_](https://arxiv.org/abs/2609.03934), identify bounded ratios on a positive semialgebraic set with the dual cone of its tropicalization. Once a fixed obstruction-profile model gives a genuine semialgebraic set of tuples ((e_2,e_3,e_4,b_2,b_3,b_4)), tropicalization could enumerate the ratio inequalities actually valid on that model. This is downstream of the canonical partition; it does not license calling tree independence polynomials Lorentzian.


The remaining seven papers offer secondary tests or analogies, not immediate bridges.

- Xiao, [_The real-rootedness of the toric g-contribution polynomials_](https://arxiv.org/abs/2609.01086): recurrence/interlacing may help only if path classes yield a fixed recurrence.

- Gantert, Löwe, and Schnappinger, [_When does propagation of chaos in the critical Curie--Weiss model stop?_](https://arxiv.org/abs/2609.02464): its conditional-mixture decomposition is a useful warning to localize where unimodality changes, but supplies no tree coefficient inequality.

- Gerling, [_The multiorbital bivariate chromatic polynomial_](https://arxiv.org/abs/2609.01866): orbit decompositions may compress symmetric matching skeletons, not the unrestricted case.

- Liu and Zhang, [_Taylor positivity of Ehrhart polynomials_](https://arxiv.org/abs/2609.03327): shifted bases and Newton inequalities suggest diagnostics around defect depth, but the Ehrhart/real-rooted hypotheses are unavailable here.

- Fouli et al., [_Graded Betti numbers of graded Möbius algebras of uniform matroids_](https://arxiv.org/abs/2609.03827): adjacent Hilbert-series language, but no current map from blocked classes.

- Langharst, Malliaris, and Roysdon, [_On cost-induced Santaló-type inequalities in Polish measure spaces_](https://arxiv.org/abs/2609.01460): even its Hamming-cube application is too analytic for the present exact defect profile without a new transference dictionary.

- Pafis, [_Cramér transform, half-space depth and threshold phenomena for convex bodies_](https://arxiv.org/abs/2608.29972): potentially relevant to a future asymptotic concentration approach, not the bounded depth-three correction.


This ordering is by exact adaptation cost and falsifiability, not by the source field's apparent proximity to graph theory.
## 7. Recommended next experiment
Implement one bounded, deterministic first-violation partition.

1. From each matching-bag instance, construct the forest poset and freeze a total order on missing transitive comparisons and their unique cover paths.

2. Enumerate every feasible defect-2, -3, and -4 partial assignment. Mark it extendable or assign a blocked one to its first violated comparison/path.

3. Require a disjoint partition and exact reconstruction of (b_2,b_3,b_4) on every tree through (n=15). Any mismatch kills the formulation.

4. For each class and useful cumulative union, test closure under deletion of erased bags. If a down-set or ballot completion exists, write its exact rank formula; otherwise record the smallest counterexample and do not invoke O-sequence machinery.

5. Test (3) at class and paired-class level. Following Pahari and the reusable part of Astra #126, seek a summandwise injection or monotone-ratio pairing before adding classes.

6. Only if the valid ratios remain opaque, tropicalize a fixed semialgebraic profile model to enumerate candidate inequalities. Do not tropicalize the unrestricted tree family before such a model exists.


The decisive outputs are therefore not a large census but a small structural table: class key, defect, count, closure result, contribution to (C), and contribution to (R).

### First bounded gate: completed

The disjoint-partition gate was implemented in `../scripts/probe_first_violation_partition_20260904.py` and replayed through \(n\le15\). It passed exact reconstruction for all 13,179 relevant trees and all three defect levels:

- 40,547,882 candidate assignments and 13,254,850 feasible assignments;
- 8,465,993 extendable assignments and 4,788,857 blocked assignments;
- every blocked assignment was assigned to exactly one canonical first obstruction;
- 1,400,970 blocked assignments had multiple available certificates, so the canonical disjointification is doing nontrivial work;
- the assigned first obstruction was unary 4,562,461 times and a pair comparison 226,396 times;
- assigned pair-certificate paths had lengths 2 through 6, with counts 206,132, 19,342, 892, 29, and 1.

The last point is a useful design warning: a later class-level experiment should try a shortest-path-first or minimal-obstruction order before treating long lexicographic paths as structural classes. That refinement, the down-closure test, and the Pahari-style ratio comparison were deliberately not started in this bounded gate.

The exact result is `../results/first_violation_partition_probe_20260904.json`. The full run took 73.89 seconds wall time (73.39 seconds user CPU) on this machine. The LLM portion fit in one focused implementation-and-validation turn; a comparable first gate should be treated as roughly 30--60 minutes of LLM work, not the 2--4 hours estimated below for the broader multi-stage experiment.

### Second structural gate: completed

The follow-up `../scripts/probe_obstruction_class_structure_20260904.py` orders unary certificates first and pair certificates by increasing quotient-skeleton path length, then tests one-bag filling relations and the relative ratio at class, prefix, group, and pair levels. It again exactly reconstructed all profiles through \(n\le15\).

The closure result separates two claims that should not be conflated.

- Individual first-certificate classes are not closed under filling an erased bag. Among 12,059,848 feasible one-bag fillings, 694,084 change the first certificate, affecting 615,159 source assignments. The smallest example is the six-vertex tree `EsOG`.
- Cumulative prefixes have zero closure failures. This is not merely empirical: filling a bag preserves every previously violated certificate and may only add earlier violations, so the first-certificate rank can only decrease. The resulting closure is in the feasible partial-assignment poset; an ordinary multicomplex or ballot-complex representation still requires a separate dictionary.

The ratio tests give a more promising Pahari-style target.

- All 46,557 nonempty tree-local class profiles are log-concave. All 41,427 profiles present at defects 2 and 4 satisfy
  \[
  c_3^2e_2e_4\ge c_2c_4e_3^2.
  \]
  There are zero class-level failures, including all 335 classes on the 223 trees with \(C<0\).
- Cumulative prefixes remain log-concave, but their relative ratio fails 125 times. All 125 failures have \(C\ge0\), so prefix closure alone does not explain the inequality.
- Of 71,445 cross-class pairs, 9,991 have a negative pair contribution. Exactly one lies on a tree with \(C<0\): an \(n=15\) tree where the first unary class contributes 3,587,661 units of relative reserve and its cross term with the second class is -2,645,000, leaving 942,661.
- Grouping by unary versus pair obstruction, and refining pair obstructions by path length, removes that exceptional interaction on the bounded adverse corpus. The 223 adverse trees split into 220 unary-only and three pair-only cases. Two pair-only cases use path length 2; one uses path lengths 2 and 3. All 224 resulting atomic group profiles have nonnegative relative margin, and the sole cross-group pair term is also nonnegative.

Consequently, on the exact \(n\le15\) corpus, the grouped diagonal terms plus grouped pair terms provide a summand-before-addition certificate of \(R\ge0\). This was a useful killable adaptation of Pahari and Astra #126, but it is superseded by the \(n=17\) counterexample in Section 5. The direct individual-class down-set route is independently closed by the \(n=6\) example.

### Proportional proof extraction: exact unary formulas

Let \(U\) be the vertices of \(T\) that belong to no maximum independent set, and put

\[
J_d=[x^{\alpha-d}]I(T-U;x).
\]

The coarse obstruction groups have the exact coefficientwise filtration

\[
e_d\ \le\ J_d\ \le\ s_d,
\qquad
u_d=s_d-J_d,
\qquad
p_d=J_d-e_d, \tag{4}
\]

where \(u_d\) is the unary group and \(p_d\) the pair group. Indeed, an independent set is assigned to the unary group exactly when it contains a vertex in \(U\); if it avoids \(U\) but is blocked, its first obstruction is a pair comparison.

There is also an exact disjoint formula for individual unary classes. Order \(U=(v_1,\ldots,v_r)\) by the unary certificate order and put

\[
F_i=T-N[v_i]-\{v_j:j<i\}.
\]

Then the class whose first forbidden vertex is \(v_i\) has

\[
c_{i,d}=[x^{\alpha-d-1}]I(F_i;x). \tag{5}
\]

This follows by deleting \(N[v_i]\) after requiring \(v_i\), and deleting the earlier forbidden vertices to enforce first occurrence. Thus the dominant unary part of the obstruction problem is now expressed by ordinary independence polynomials of smaller forests, rather than opaque certificate counts. The formula-augmented replay checked (5) for 41,685 unary certificates and both identities in (4) on all 13,179 trees, with zero mismatches.

The individual unary-class inequality suggested by the \(n\le15\) data becomes

\[
\bigl([x^{\alpha-4}]I(F_i;x)\bigr)^2e_2e_4
\ge
[x^{\alpha-3}]I(F_i;x)\,[x^{\alpha-5}]I(F_i;x)e_3^2. \tag{6}
\]

No general forest inequality presently proves (6), and unrestricted forest independence polynomials are not log-concave. More importantly, proving (6) class by class would not establish the aggregate target: the \(n=17\) counterexample has two unary classes with profiles \((1,8,28)\) and \((0,0,16)\). Their individual relative margins are nonnegative, but their cross-class contribution is \(-31{,}719{,}424\), overwhelming the first class's reserve 28,250,624 and producing \(R=-3{,}468{,}800\).

The exact filtration (4) and unary formula (5) remain valid and potentially useful for a direct cross-reserve argument. The relative-ratio consequence, individual-to-group summation, and natural grouped-summand route do not.

The proved extendable reserve plus blocked log-concavity also cannot rescue this route algebraically. At \(\alpha=5\), the positive numerical profiles

\[
(e_2,e_3,e_4)=(10,20,20),
\qquad
(b_2,b_3,b_4)=(1,10,100)
\]

satisfy the exact extendable reserve bound and blocked log-concavity, yet have \(C=-620\) and \(R=-20{,}000\). These numbers are not asserted to arise from a tree; they show that any successful proof must use additional tree structure and must target the direct cross-reserve inequality rather than (3).

### Direct cross-reserve candidate

The proved Pascal estimate gives

\[
e_3^2-e_2e_4
\ge
\frac{5\alpha+17}{27(\alpha-3)}e_2e_4. \tag{7}
\]

Since \(2e_3b_3\ge0\), a stronger sufficient condition for the live cross-reserve target is

\[
e_2b_4+e_4b_2
\le
\frac{5\alpha+17}{27(\alpha-3)}e_2e_4, \tag{8}
\]

or equivalently

\[
\frac{b_2}{e_2}+\frac{b_4}{e_4}
\le
\frac{5\alpha+17}{27(\alpha-3)}. \tag{9}
\]

This discards the favorable middle term rather than requiring a relative log-concavity comparison between the whole profiles. The optimized probe tests (8) only where it is needed, namely \(C<0\).

Through \(n=18\), all 1,539 negative-cross cases satisfy (8). The census covers 204,995 relevant trees and has zero failures of (8), the weaker Pascal-plus-full-\(C\) condition, the actual cross-reserve inequality, blocked-profile depth-three log-concavity, or the full depth-three inequality. The worst adverse \(-C\) consumes \(9/35\) of the guaranteed lower reserve; the smallest scaled margin in the Pascal-plus-full-\(C\) test is 585.

Equation (8) is a live cross-reserve component. It asks for a weighted bound on the blocked densities at defects 2 and 4 and fits the original obstruction-path charging architecture. The first-violation partition and exact unary filtration remain possible tools for proving it. The targeted experiment below shows that it is not sufficient by itself, because blocked-profile log-concavity is false in the intended rung. This is bounded evidence, not an all-tree claim, and neither the \(n=17\) auxiliary counterexample nor the \(n=18\) checks refute or prove Erdős #993.

The exhaustive \(n\le18\) census does not enter the intended high-independence, low-deficiency rung: the live finite reduction has \(\alpha\in\{17,18,19\}\), deficiency \(\delta=2\alpha-n\le5\), and \(33\le n\le38\). A bounded follow-up therefore generated random labeled quotient trees, singleton locations, and endpoint ports, conditioned on the proposed bag matching being maximum. It accepted 50 trees in each of the 12 admissible \((\alpha,\delta)\) cells, for 600 accepted trees from 1,560 trials. None had \(C<0\), so none exercised (8). There were also no blocked-profile depth-three, cross-reserve, full depth-three, or unimodality failures, but those zeroes are only biased random evidence. The run took 240.37 seconds wall time. Its value is diagnostic: untargeted bag-tree sampling is too unlikely to reach the adverse-cross regime, so scaling the same sampler is not warranted. Any future computational test of (8) should deliberately optimize for \(-C\) while retaining the actual-rung constraints.

The exact random-probe record is `../results/cross_reserve_target_scope_probe_20260904.json`; its generator is `../scripts/probe_cross_reserve_target_scope_20260904.py`. Before the final run, a 12-tree smoke test took 4.78 seconds and a 120-tree pilot took 47.27 seconds; both likewise contained no \(C<0\) sample. Total local compute for developing and running this target-scope probe was 292.42 seconds.

### Targeted witness lifts and the corrected joint target

The first \(n=17\) negative-cross witness has \((n,\alpha,\delta)=(17,11,5)\). Appending a pendant \(P_2\) raises the order by two and the maximum-matching size by one, hence raises \(\alpha\) by one and preserves \(\delta\). Eight such matched-pair lifts therefore land exactly in the live \((n,\alpha,\delta)=(33,19,5)\) cell. The targeted probe preserved the witness's chosen maximum matching, varied the attachment vertices and ports through fan, chain, distributed, and adaptive mutations, deduplicated by exact unrooted-tree isomorphism, and rechecked maximality on every lift.

Among 500 nonisomorphic lifts, 57 have \(C<0\). Equation (8) has zero failures; the worst endpoint load is only 0.0250234 of the guaranteed Pascal reserve. There are zero cross-reserve and full depth-three failures and zero non-unimodal independence polynomials.

However, nine lifts refute blocked-profile depth-three log-concavity. The first exact witness has

\[
(e_2,e_3,e_4)=(325632,1245184,3337984),
\qquad
(b_2,b_3,b_4)=(1,16,4216),
\]

so

\[
b_3^2-b_2b_4=-3960.
\]

Its mixed term is \(C=-1{,}336{,}356{,}608\), but its cross-reserve margin is \(462{,}192{,}431{,}360\) and its full depth-three margin is \(462{,}192{,}427{,}400\). A fresh maximum matching reconstructed the same \(e\), \(b\), and \(s\) profiles directly from the saved graph6 certificate. This refutes only the auxiliary blocked-profile lemma.

The live sufficient bridge must therefore treat the blocked and mixed deficits together. Put

\[
D=(b_3^2-b_2b_4)+C.
\]

If \(D\ge0\), the full margin is already nonnegative. When \(D<0\), a deliberately stronger sufficient condition is

\[
e_2b_4+e_4b_2+b_2b_4
\le
\frac{5\alpha+17}{27(\alpha-3)}e_2e_4. \tag{10}
\]

Indeed, (10) discards both favorable terms \(2e_3b_3\) and \(b_3^2\) before charging the remaining endpoint load to the Pascal reserve. The 500-lift corpus has 57 cases with \(D<0\), exactly the 57 negative-cross cases, and zero failures of (10) there. The worst combined load is 0.0250394 of the guaranteed reserve. Equation (10) fails 380 times outside the \(D<0\) regime, so it is a conditional bridge only and must not be promoted to a universal inequality.

The final corrected run took 56.66 seconds wall time. Exact generator and certificate: `../scripts/probe_cross_reserve_witness_lifts_20260904.py` and `../results/cross_reserve_witness_lifts_20260904.json`. The live proof-frontier graph was updated to mark `blocked_profile_depth3_lc` refuted and replace the old two-lemma assembly by the joint Pascal-margin target.

The formula-augmented final pass took 110.81 seconds wall time. Four full structural development passes took 440.37 seconds (about 7.3 minutes) of laptop compute in total.

## 8. Time and compute estimate on this machine
Machine: 2026 MacBook Pro, Apple M4 (4 performance plus 6 efficiency cores), 32 GB RAM, Python 3.14.2.

Measured 4 September:

- the optimized existing (n\le15) probe: 4.73 seconds wall time;

- the optimized extension through \(n\le17\): 46.98 seconds wall time;

- the optimized direct-reserve kill test through \(n\le18\): 159.83 seconds wall time;

- the actual-rung random bag-tree probe: 600 accepted trees from 1,560 trials, 240.37 seconds wall time;

- the targeted actual-rung witness-lift probe: 500 nonisomorphic lifts, 56.66 seconds wall time;

- a deliberately naive explicit enumerator over all defect-2/3/4 assignments: 40,547,882 candidates, 13,254,850 feasible assignments, exact reconstruction against the polynomial DP, 39.25 seconds wall time on one process.


Likely cost for the recommended bounded experiment:

- **LLM working time:** about 2--4 hours for implementation, exact validation, counterexample inspection, and a concise result note. A first decisive pass/fail signal should appear in 45--90 minutes. If the partition survives and the goal expands to a proof-quality Pahari-style summand comparison, budget another 4--12 hours of mathematical work; success is not implied.

- **Local compute per full (n\le15) run:** approximately 1--3 minutes for an unoptimized Python implementation with path certification and closure tests; probably below one minute after caching and bitset cleanup. Memory should stay comfortably below 1 GB.

- **Cumulative laptop compute during development:** roughly 10--45 minutes for repeated full runs. CPU time is not the bottleneck at this scope.


The optimized profile scan through \(n=18\) is complete. Extending the full assignment-class census to that order would be a different, much more expensive validation tier and is not warranted by the current target. Enumerating the full original \(n=33\)--38 tree window locally is not part of this experiment and should not be inferred from these timings.
## 9. Warrant boundary
- The 13,179-instance statements and timings above were replayed locally in exact integer arithmetic on 4 September 2026.

- The Astra solutions are source inputs and proof-method exemplars; their Lean repositories were not independently rebuilt in this audit.

- Equation (3) is false; its first exact counterexample in this scan occurs at \(n=17\).

- Equation (8) has no failure through \(n=18\), but remains an unproved conditional bound.

- The 600-tree actual-rung random probe found no \(C<0\) input and therefore supplied no direct evidence for or against Equation (8).

- The targeted lift probe found 57 \(C<0\) inputs in the actual rung and no failure of Equation (8), but it refuted blocked-profile log-concavity in nine cases.

- Equation (10), the corrected conditional joint bound, has no failure among the 57 lifts with adverse combined correction; it remains unproved.

- None of the weekly papers directly proves either open matching-bag frontier node.

- No source here changes the submission or manuscript claims.

- The obstruction-class experiment supplies exact formulas and a formal prefix-closure observation. Its relative-ratio proof route is closed by the \(n=17\) certificate.

---
comments:
  c1:
    body: this sounds like a go, but I don't want to burn a big chunk of my weekly
      chatGPT tokens budget on this
    by: user
    at: 2026-09-04T13:40:49.769Z
  c2:
    body: I ran only the bounded disjoint-partition gate and stopped before down-closure, class-ratio, proof-search, or tropicalization work. It passed in 73.89 seconds of local compute.
    by: Codex
    at: 2026-09-04T13:49:34.000Z
    re: c1
  c3:
    body: On request I continued through the closure and summandwise-ratio gate. Individual closure fails, cumulative-prefix closure is automatic, and the natural grouped ratio decomposition has no adverse-case failure through n<=15. I stopped before proof development or tropicalization.
    by: Codex
    at: 2026-09-04T14:07:04.000Z
    re: c1
  c4:
    body: The proportional proof follow-up derived exact unary-class and coarse-group formulas in terms of smaller forest independence polynomials, verified them with zero mismatches, and stopped at the first unsupported relative coefficient inequality.
    by: Codex
    at: 2026-09-04T14:23:11.000Z
    re: c3
  c5:
    body: The proportional n<=17 validation found an exact counterexample to C<0 => R>=0. The grouped relative-ratio route is closed; the original cross-reserve remains positive. Exact unary filtration formulas are retained for any future direct attack.
    by: Codex
    at: 2026-09-04T14:32:15.000Z
    re: c4
  c6:
    body: The direct cross-reserve follow-up found a stronger sufficient blocked-density bound with zero failures in 1,539 adverse cases through n<=18. This is an auxiliary proof target only; no counterexample to Erdős #993 was found.
    by: Codex
    at: 2026-09-04T14:42:58.000Z
    re: c5
  c7:
    body: The actual-rung random probe accepted 600 trees across all 12 alpha=17--19, deficiency-at-most-5 cells, but found no C<0 samples. It therefore did not test the conditional endpoint bound. Untargeted scaling stops here; no counterexample to Erdős #993 was found.
    by: Codex
    at: 2026-09-04T15:02:00.000Z
    re: c6
  c8:
    body: Targeted matched-pair lifts produced 57 adverse actual-rung cases. The endpoint bound survives, but nine lifts refute blocked-profile log-concavity. The corrected live bridge charges the mixed and blocked deficits jointly and survives all 57 adverse lifts. No counterexample to Erdős #993 was found.
    by: Codex
    at: 2026-09-04T15:10:13.000Z
    re: c7
