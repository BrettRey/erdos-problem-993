# Independent mathematical review B

Reviewed artifact: `source/b1_zero_forbidden_core_2026-09-05.md` in this review directory, SHA-256 `4d77890b084b1a7aec056217b84d2606f17c78b4ec5772aff2401426d300a8bb`.

All line references below refer to that frozen artifact unless another file is named. I read the complete artifact, the frozen forest-poset note, both frozen audit scripts, and the relevant local Pascal-smoothing argument and Lean theorem statements. I did not read other reviews or prior conversation, modify source files, or dispatch another reviewer.

## Claims and verdicts

The first result is a strict upper-layer log-concavity theorem for a specified family of finite trees: if `33 <= n <= 38`, `alpha` is 17, 18, or 19, `2*alpha-n <= 5`, and every independent `(alpha-1)`-set extends to a maximum independent set, then `s3^2 > s2*s4`. Its new structural content is that the allowed induced forest is well-covered and every forbidden vertex has at least three forced neighbours. This confines the forbidden set to at most four vertices. Twelve parameter cases are absorbed by a coefficient-density estimate; the last case is reduced to 1,842 rooted forests and finished using a positive lower bound on `b3`.

The second result is the blocked-set incidence inequality `6*b3 >= (alpha-4)*b2`, without the assumption `b1=0`. Together with the extendable incidence bound and Pascal reserve, it yields a sufficient density threshold for strict depth-three log-concavity. For `alpha=17,18,19`, those thresholds are respectively `5/21`, `11/45`, and `1/4`. This is a sufficient condition, not a proof that every tree in the remaining branch meets it.

| Result | Verdict |
| --- | --- |
| Section 1 structural lemma for finite forests with `b1=0` | **Valid as written.** |
| Sections 1–7 strict depth-three theorem in the stated tree window | **Valid as written as a computer-assisted mathematical proof.** The finite certificate was replayed and independently recomputed. |
| Section 7 bounds `b4/e4 <= 52513/217404` and `b3 >= 11*b2` | **Valid as written for the connected tree family specified there.** |
| Section 9 incidence inequality and sufficient density threshold | **Valid as written in the stated window.** A standalone general statement should explicitly say finite forest and `alpha >= 4`, and retain the stated rank restrictions for its strict-threshold consequence. |

I found no mathematical error requiring a change to either intended result. The changes recommended below concern explicit theorem scope, definitions, and reproducibility. They are not evidence against the proof.

## Reconstruction and examination of the proof

### 1. The structural lemma is sound

At lines 46–54, a maximum matching and its unmatched vertices partition the graph into exactly `alpha` bags. Every independent set occupies at most one vertex per bag, and every maximum independent set occupies every bag. Thus an unmatched singleton is forced; a matched bag either has two allowed endpoints or has one forbidden endpoint and one forced endpoint. This justifies `|C|=r+delta`. It does not assume that all maximum sets make the same choice in a flexible pair.

The splicing argument at lines 56–67 is valid. Removing the endpoints of one matching edge leaves components made from whole remaining bags: another matching edge cannot cross two of those components. A path connecting the chosen neighbours outside the removed edge would give a cycle. Restrictions of maximum sets on the separate components can therefore be combined, each contributing the number of bags in its component. The resulting independent `(alpha-1)`-set blocks both endpoints of the omitted bag. That is precisely a contradiction to `b1=0`.

The extension construction at lines 69–72 also works. Each unoccupied flexible matching pair has a degree-one endpoint in the allowed induced graph. Adding those endpoints cannot introduce an edge between different chosen pairs, and the forced vertices are isolated there. Consequently every independent set avoiding `U` is extendable. This is the needed identification of the entire extendable polynomial with `I(H;x)`.

At lines 74–81, acyclicity ensures that a forbidden vertex has at most one neighbour in each component of `H`. A nonforced vertex in a component can be avoided by a maximum set of that component; otherwise it would be globally forced. Componentwise choices therefore construct the claimed set of size `alpha-d_C(v)`. Adding `v` produces the three separate contradictions for `d_C(v)=0,1,2`. Finally, the forest edge bound on `U union C` gives `3r <= 2r+delta-1`. All four structural conclusions follow for disconnected forests as well as trees.

### 2. The coefficient-density estimates are justified

For nonempty `S subset U`, the bipartite subgraph on `S union N_C(S)` has at least `3|S|` edges and is a forest, giving `|N_C(S)| >= 2|S|+1`. Ignoring the restrictions on flexible vertices and independence within `S` can only increase the coefficient counts. This gives the polynomial upper bound at lines 132–138, and extracting its upper coefficients gives equations (3) and (4).

The factor two in equation (5) is correct: an independent set one layer lower leaves `j+1` of the fixed two-vertex matching bags empty, with at most two choices per empty bag. Thus `p_j / (2^-j * binom(m,j))` is nondecreasing. Equation (7) then follows by writing that monotone sequence as a nonnegative combination of suffix indicators. Every denominator is positive in the nonempty-forbidden window, because `m >= 10` and the term at `j=4` is positive.

I independently recomputed all 13 parameter triples and their five rational suffix ratios. The 12 cases with `r<=3` all satisfy the stronger common comparison with `7/27`. The largest load is exactly `170838/691067` at `(alpha,delta,r)=(19,4,3)`, and its cross-multiplied slack is exactly `224843`. The sole relaxation failure is `(19,5,4)`, with load `1012328/3244617`. There is no omitted admissible positive-`r` case at `alpha=17`: its window constraints force `delta<=1`.

The final algebra at lines 232–240 is correct. Discarding the nonnegative `2*e3*b3+b3^2` term is legitimate when the joint endpoint load is already below the reserve.

### 3. The four-forbidden reduction covers the tree family

The equality case of the edge count at lines 257–265 forces a connected 13-vertex core with precisely 12 forbidden-to-forced edges. Every forbidden vertex has exactly three forced neighbours, and there are no forbidden-to-forbidden edges. Because the original graph is a tree, every component of `K` has exactly one edge to this core. A forced attachment endpoint is impossible; two attachments would give a cycle.

The central forced vertex plus two private forced leaves per forbidden vertex attains `|N_C(S)|=2|S|+1` simultaneously for every nonempty `S`. Its replacement therefore supplies the stated coefficientwise domination while preserving `E`.

I expanded equation (10) directly. Its difference is exactly `a*x*(A_i-Q_i)*(A_j-Q_j)`. Both polynomial differences count sets containing at least one deleted attachment root, so every coefficient is nonnegative. Multiplication by the unaffected factors and subtraction of the unchanged `a^4*P` preserve the domination. Concentrating attachments is therefore justified coefficient by coefficient.

The corona description at lines 366–369 follows from the pendant perfect matching: a selected leaf endpoint has no edge to another matching pair, so all such edges join the selected base vertices. The root-switch injection at lines 371–377 is also valid. A switched set contains the private leaf; an unchanged set avoids it, so the two types of image cannot collide. Different components have separate attachment roots, allowing the injections to be applied componentwise.

Finally, rooted forests of total order ten are in bijection with rooted trees of order eleven, by adjoining one new root. The verifier compares full canonical-code sets, and my separate calculation also recovered exactly those 1,842 codes. The reduction preserves the extendable polynomial attached to each representative; it does not require a single polynomial to dominate representatives with different base forests.

The sharp ratio was recomputed as `157539/652212 = 52513/217404`, attained by the ten-vertex star rooted at its centre. Every frozen coefficient row agreed with the independent calculation described below. The 17 failures of the stronger discarded joint bound were also reproduced. Those failures do not undermine the argument actually used.

### 4. The positive defect-three terms are used correctly

For the original tree, deleting the flexible neighbours of one forbidden vertex leaves independence number ten: each deleted neighbour lies in a different flexible component and is individually avoidable by a maximum set. A nine-set occupies nine original matching bags and has at most two one-vertex extensions, so the deletion count gives `q_i,1 >= 5*q_i,0`.

A defect-two blocked set has exactly one forbidden vertex, and the six other forced vertices are available. Hence `b2=sum(q_i,0)` and `b3 >= sum(6*q_i,0+q_i,1) >= 11*b2`. This is proved for the original tree, not inferred from its dominating representative.

The extendable incidence inequality `16*e3 >= 4*e4` has the correct direction: every extendable 15-set has at least four extensions, while each extendable 16-set has precisely sixteen deletions. Substitution gives the two coefficients in (13), `7/27-v` and `9/2-v`, with `7/27-v=3851/217404>0`. Strictness follows because `e2*e4>0`, independently of whether any blocked coefficient is positive.

### 5. The general incidence argument and density threshold are sound

The frozen forest-poset note, Section 2, supplies the projection fact used at lines 524–530. A partial assignment fails to extend only if it specifies a forbidden value on a constant coordinate or specifies two free values violating a comparison in the full induced poset. Taking the down-closure of a consistent partial ideal supplies an extension. Thus a blocked independent set has a blocked singleton or pair contained in it.

Fixing such a witness inside a blocked `(alpha-2)`-set leaves at least `alpha-4` deletions that preserve the witness. A blocked `(alpha-3)`-set leaves exactly three matching bags empty, including possible singleton bags, so it has at most six one-vertex extensions. Counting the same deletion incidences on both sides proves (14).

The subsequent algebra is correct:

`2*e3*b3 >= [4*(alpha-4)/(3*(alpha-3))]*e4*b2`.

After subtracting `e4*b2+b2*b4`, the coefficient is exactly `h_alpha-v`, where `h_alpha=(alpha-7)/(3*(alpha-3))`. The reserve leaves `(c_alpha-v)*e2*e4`. Also `c_alpha-h_alpha=4*(20-alpha)/(27*(alpha-3))`. For the three stated ranks this is strictly positive, so the weak condition `v<=h_alpha` indeed proves a strict margin.

## Changes recommended for an exact theorem and reusable proof

1. **State the graph and rank hypotheses explicitly in Section 9.** Lines 517–558 inherit the forest setting and current window, so the intended claim is sound. A standalone theorem should nevertheless say “finite simple forest with `alpha>=4`” before (14)–(15), and separately restrict the threshold conclusion to `alpha in {17,18,19}`. This avoids negative defect sizes, natural-number subtraction conventions, and the zero denominator at `alpha=3`. For formalization, the threshold can be written without division as `3*(alpha-3)*b4 <= (alpha-7)*e4` in the stated ranks.

   The forest assumption cannot simply be dropped. I checked the finite graph `C9 disjoint union K1`: it has `alpha=5`, `(b0,b1,b2,b3,b4)=(0,3,3,0,0)`, and three minimal blocked triples. Equation (14) would say `0>=3`. This is a scope counterexample, not a counterexample to the forest theorem.

2. **Retain connectedness in the Section 7 finite theorem.** Lines 262–265 correctly use an original tree. Section 7.2’s “one distinguished base vertex in each component” is not valid for an arbitrary disconnected forest, whose flexible components may have no attachment. The numerical bound itself fails in that larger class. The disjoint union of the 13-vertex shared-four core and ten isolated edges has `n=33`, `alpha=19`, `delta=5`, `r=4`, `b1=0`, but `b4/e4=2264/8793 > 52513/217404`. I computed this example. It does not refute the main theorem, which is explicitly about trees. Add “Let T be a tree” to the Section 7 hypothesis paragraph to protect the statement when extracted for Lean.

3. **Define `D` in this document.** The symbol first appears at lines 242–243 and becomes part of the remaining regime at lines 560–565, but no displayed definition is supplied. Insert

   `D = 2*e3*b3 + b3^2 - e2*b4 - e4*b2 - b2*b4`,

   so that `s3^2-s2*s4 = (e3^2-e2*e4)+D`. With that definition, excluding `D>=0` from the unresolved regime is justified by the positive reserve. The displayed negative and positive values of `D` that I arithmetically recomputed agree.

4. **Make elementary side conditions explicit once.** Before the ratio calculations, record `e_d >= binom(alpha,d)>0` for `0<=d<=4` when `alpha>=4`: the subsets of one maximum set already provide these sets. At the start of Section 3, say that its density calculation assumes `r>0`; the separate `r=0` argument at line 230 already handles the other case. These observations justify all divisions, strict conclusions, and the inequality `c>=2r+1` used there.

5. **Identify the dependency theorem and complete the replay package.** Lines 202–210 correctly identify the reserve but refer indirectly to a status graph. The exact denominator-cleared imported statement is

   `32*(alpha-2)*e2*e4 <= 27*(alpha-3)*e3^2`, for `alpha>=4`.

   I inspected the local written Pascal proof and the matching theorem statements `erasureProfile_depth_three` in `PascalBridge.lean` and `erasure_depth_three_reserve` in `TreeCodeBridge.lean`. The latter supplies the tree-code version. I did not build or kernel-replay Lean, so this review is not a new formal-verification claim.

   The frozen Python files depend on repository modules outside `source/`, and their `REPO = Path(__file__).resolve().parents[1]` calculation does not locate the repository when run from the frozen directory. Their command-line `main()` functions also write certificates. I replayed their `run()` functions in memory after setting `REPO` to the actual repository and loading the frozen first audit under its import name. A durable standalone package should freeze the imported helpers and archive inputs, identify the Python/NetworkX environment, and provide a read-only replay entry point. This is reproducibility work; the current replay and independent reconstruction succeeded.

## Checks actually performed

- Executed the frozen first audit as `run(14)` without calling its writing `main()`. All result fields matched the frozen JSON, except runtime. This reran 5,446 trees, 105 disjoint forests, 109 archived graph6 profiles, and all 13 density records. It included independent-set enumeration on the small graphs and 5,438 eligible small-tree incidence checks. Runtime was approximately 13.923 seconds.
- Executed the frozen concentration audit as `run(1000,3)`, loading the frozen first audit as its dependency. Every non-runtime result field matched its frozen JSON. This reran all 1,842 representatives, 2,692 small decorations, and 1,000 deterministic order-33 decorations. Runtime was approximately 3.489 seconds.
- Independently generated all rooted eleven-vertex tree types using NetworkX’s unrooted tree enumeration and every possible root, then removed the root to obtain each rooted forest. Independently enumerated subsets of the ten base vertices to obtain the corona polynomials `P,Q`. For every defect `d=2,3,4`, I extracted blocked coefficients from the direct forbidden-subset expansion

  `sum_{t=1}^4 x^t*(1+x)^(8-2*t) * (binom(3,t)*P + binom(3,t-1)*Q)`.

  This agrees with equation (11) but does not call the supplied concentration-polynomial implementation. All 1,842 canonical codes and all seven numerical entries in every frozen coefficient row matched. The maximum ratio, its slack, and all 17 stronger-bound failures matched.
- Independently recomputed the 13 rational parameter cases and all cone-ray ratios using `Fraction` arithmetic.
- Independently enumerated every independent set in all 986 nonisomorphic trees of orders 2–12, formed the blocked family directly from maximum-set containment, and checked every inclusion-minimal blocked set. Every such witness had size at most two. The eligible instances satisfied (14).
- Checked the two scope counterexamples above. Recomputed the Section 8 displayed correction `-178212783` and full margin `94423068753`, and the Section 7 extremizer correction `2637714832`, from their displayed profiles.

I did not rerun the full repository test suite, the Section 8 leaf-reattachment search, the Section 9 adaptive search, or a Lean build. Their historical execution claims are not independently certified by this report. The adaptive search supplies no universal premise to either accepted proof. No external literature source was used as verified evidence; the proof dependencies discussed here were local sources that I actually read.

## Final assessment

The two intended mathematical results survive this independent review. The first is a valid computer-assisted theorem for the stated tree window; the second is a valid forest incidence lemma with the stated sufficient thresholds. There is no demonstrated closure of the unrestricted `b1>0` branch. The most useful next preparation for formalization is to extract explicit theorem statements with their rank and connectedness hypotheses, spell out positivity and `D`, and make the finite enumeration certificate reproducible independently of the live repository layout.
