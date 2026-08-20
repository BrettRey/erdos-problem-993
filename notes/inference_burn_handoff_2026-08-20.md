# Erdős 993 inference-burn handoff

Date: 2026-08-20. This file records the durable endpoint of the intensive
attack session. It does not alter `paper/main_v2.tex` or any submitted
surface.

## Proved or proof-drafted advances

1. **LAW-V with one trivial hub.** The path-extension likelihood-ratio
   order, the spider Turán gap, and two shifted LR summands prove the law
   for arbitrary arms. See `law_v_one_hub_proof_2026-08-20.md`.
2. **Adjacent two-hub base.** A disjoint clan-block decomposition repairs
   the collision in the naive map injection and gives a complete proof
   draft that the adjacent polynomial \(B\) is log-concave. See
   `two_hub_B_logconcavity_proof_2026-08-20.md`; the failed injection is
   retained in `two_hub_clan_cancellation_attack_2026-08-20.md`.
3. **All connector lengths.** Odd connectors inherit the same clan
   argument. Even connectors reduce to one centrally-unimodal Laurent
   block proved by a coefficient-range split and Vandermonde bounds. This
   yields a proof draft for log-concavity of every tree with at most two
   vertices of degree at least three. See
   `c2_connector_clan_reduction_2026-08-20.md`.
4. **Maximum-matching bag model.** Contracting a maximum matching turns a
   tree into an \(\alpha\)-bag binary CSP with
   \(\delta=2\alpha-n\) singleton bags. The current rung has 17--19 bags
   and at most five singletons. Minimal blocked residuals through depth
   four are alternating-port paths. See
   `matching_bag_csp_attack_2026-08-20.md`.
5. **Forest-poset gauge.** Maximum assignments form an order-ideal code of
   a poset whose cover graph is a forest, times constant coordinates. Its
   projection polynomial is the independence polynomial of the canonical
   bipartite poset graph and equals the full-whiskering transform of the
   antichain polynomial. See
   `matching_bag_poset_reduction_2026-08-20.md`.
6. **Universal Pascal-smoothing lemma.** For every downward-closed family,
   normalized face densities decrease. Two applications of Pascal's
   identity then give the abstract ratio
   \(S_m^2\ge(8/9)S_{m-1}S_{m+1}\). Hence a full-whiskering transform has
   all upper-tail inequalities through defect depth eight at every rank
   and is wholly log-concave through
   rank 33. The extendable profile in the remaining \(\alpha\le19\) rung
   is completely proved, with
   \[
   e_3^2\ge\frac{32(\alpha-2)}{27(\alpha-3)}e_2e_4.
   \]
   See `pascal_smoothing_shadow_lemma_2026-08-20.md`.

## Killed or superseded routes

- The naive clan-level injection is noninjective; its exact collision is
  preserved rather than hidden.
- Ultra-log-concavity of arbitrary binary erasure shadows is false. The
  exact seven-bit counterexample has profile
  \((1,14,80,187,220,145,52,8)\), failing by 56 in the normalized
  inequality at index 6.
- The proposed forest-poset log-concavity lemma is no longer an open step
  in the bounded-rank problem. The universal Pascal lemma proves more than
  the current application needs.
- Full-whiskering log-concavity at arbitrary rank remains a separate open
  problem; the bounded-rank lemma must not be advertised as resolving it.

## Exact remaining bottleneck

Write

\[
s_d=e_d+b_d,
\]

where \(s_d=i_{\alpha-d}(T)\), \(e_d\) counts partial assignments that
extend to a maximum assignment, and \(b_d\) counts blocked assignments.
The term \(e_d\) is finished. The next proof task is solely to express
and bound \(b_2,b_3,b_4\), whose minimal cores are the already classified
alternating paths, and charge their overlap-corrected contribution against
the explicit reserve above using \(\alpha\le19\) and \(\delta\le5\).

Before treating the C2 theorem as final, independently audit the clan-state
partition/normalization and the even-connector coefficient-range split.

## Short replay commands

```text
python3 verify_law_v_one_hub_20260820.py
python3 verify_two_hub_clan_attack_20260820.py
python3 verify_local_li_block_20260820.py
python3 verify_two_hub_B_logconcavity_20260820.py
python3 verify_matching_bag_csp_20260820.py
python3 verify_erasure_shadow_poset_20260820.py
python3 verify_pascal_smoothing_20260820.py
```

The connector core is included in
`verify_two_hub_B_logconcavity_20260820.py` (14,400 cases). No long census
was launched in this session.
