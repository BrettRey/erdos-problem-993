# Will Blair rooted-bush comparison

Date: 2026-07-16
Status: exact comparison of all publicly materialized witnesses; grammar-level comparison of the full V2 search space; no duplicate enumeration launched

## Decision

Do **not** launch another broad rooted-bush or spider-bouquet enumeration.

Blair's V2 rooted-bush grammar is an exact special case of the existing local spider-bouquet grammar in `scripts/valley_search.py`. The 25 witnesses actually serialized in the public result were reconstructed and checked exactly. They add mixed-gadget parameter points and five improve the available comparator at exactly the same order, but none improves the size-unrestricted maximum or the order-at-most pressure Pareto frontier. They therefore add useful bounded data without adding a new recursive mechanism.

If full archival deduplication later matters, the bounded next step is to materialize Blair's missing 4,420 polynomial records once and feed them through the importer. That is a corpus-recovery task, not a reason to repeat the overlapping search.

## Provenance and public scope

The source is Will Blair's [`verified-combinatorics/erdos-993`](https://github.com/willblair0708/verified-combinatorics/tree/main/erdos-993) directory.

- Repository commit checked: `3f84ebc061123ec29a5147204b269ab83162af34`.
- Last commit touching `erdos-993`: `59911f92a9b47aee1aa69a055edf18d1830d36da`, dated 2026-06-16.
- `results.json` SHA-256: `f13c19edf3759304ea6c0696e9a8ddf3d0170c4de5e38fbd89d9d97c6073d911`.
- The comparator reran the upstream bounded verifier and captured its command, output, and successful return code in the JSON. That verifier checks its kernel, two literature seeds, a bounded forest resweep, and its scanner; it explicitly does not rerun V2's 112,916-tree scan.

The banked V2 report says that 112,916 single trees were scanned and 4,445 distinct non-log-concave seed polynomials were found. It also reports a larger product/power/path-product surface. The committed `results.json`, however, does **not** serialize all 4,445 witnesses. It contains summary counts, two explicit order-26 seed trees, and only the 25 most severe non-log-concave rooted bushes.

This audit did not run `search_993_v2.py` or `generate_results.py`. Statements about all 4,445 objects are therefore limited to what follows from the common grammar and the upstream summary; exact graph-by-graph collision counts below apply to the 25 published records.

The reported forest total *is* independently reproducible without regenerating the library. V2's loops give

```text
C(81,2) + C(82,3) + 4,445*19
  + 80*(C(16,1) + C(17,2) + C(18,3))
= 253,695.
```

This verifies the summary's product/power/path-product arithmetic, not the underlying 4,445-library cardinality.

## Grammar translation

A Blair label has the form

```text
bush(token_1,...,token_r)|n=N.
```

The tokens translate to local spider gadgets as follows.

| Blair token | Local gadget attached to the common root |
|---|---|
| `u c p L` | `c` equal center-to-leaf legs of length `L+1` |
| `s c p L` | `c-1` legs of length `L+1` and one leg of length `L+2` |

The local `bouquet_adj` and `bouquet_poly` constructors already allow an arbitrary list of leg lengths for every spider gadget. They also allow additional root paths and leaves that V2 does not use. Hence every V2 bush is literally in the local grammar; this is containment by construction, not an empirical resemblance.

Two familiar slices are visible immediately.

- Repeated `u c p1` gadgets are Galvin trees.
- Three uniform `p1` gadgets with child counts `3,m,n` are the unstarred Li/Kadrawi--Levit family.

There is also an upstream scope mismatch. V2's `s` token lengthens one leg by **one** edge. The actual `T*` builder in upstream `search_993.py`, and the local `make_li_tree(..., starred=True)`, lengthen the distinguished leg by **two** edges. The V2 bank therefore does not literally contain the advertised starred family, although the separate `seed_Tstar_3_3_4.json` uses the correct builder.

The later `search_993_v3_wide.py` has wider ranges but no matching committed result artifact. The reported 112,916/4,445 figures belong to V2 and should not be attributed to V3.

## Exact comparison performed

The rerunnable tool is `scripts/compare_blair_rooted_bush.py`. It performed the following checks.

1. Parsed all 25 public labels and reconstructed each adjacency list.
2. Recomputed every independence polynomial twice: by the generic tree DP and by the local bouquet formula.
3. Verified the stored order, exact coefficient vector, unimodality/non-log-concavity classification, and rounded severity.
4. Replayed 5,823 committed family-grid records:
   - 759 Galvin;
   - 64 Bautista--Ramos;
   - 2,500 Li;
   - 2,500 Li-star.
5. Ingested 747 unique retained graph6 trees, representing 992 source references, from the eleven Ramos--Sun summaries and the exact `n=26,27,28` artifacts.
6. Replayed all 41 explicit-edge cases in `results/literature_root_stress_20260716.json`, independently recomputing and checking their stored exact polynomials. These include the Jerrum--Patel and Bautista root-stress lanes and add four order-matched candidates to the unrooted-tree comparison.
7. Replayed eight saved evolutionary trees from every committed `best_evolutionary_tree*.json` artifact. The two polynomial-only `lc_breaker_evo_n26*.json` records were resolved exactly to the retained order-26 Li witness already in the index.
8. Built an index of 6,619 local tree records, 4,128 exact coefficient tuples, and 3,800 exact metric signatures.
9. Compared tree identity by nauty `labelg`, polynomial identity by exact integer tuple, and post-descent signatures by reduced rational arithmetic.
10. Checked all pairs and triples of indexed components, direct powers `2,...,20` of every indexed component, and every product of one to three path polynomials `P_1,...,P_16` with an indexed residual.

The two standalone seed files were also replayed from their edge lists:

- `seed_T_3_4_4.json` is exactly local unstarred Li `(m,n)=(4,4)`;
- `seed_Tstar_3_3_4.json` is exactly local Li-star `(m,n)=(3,4)`.

## Results for the 25 public rooted bushes

| Test | Result |
|---|---:|
| Labels/orders/polynomials exactly verified | 25 / 25 |
| Exact unrooted-tree collisions with the local index | 1 / 25 |
| Exact independence-polynomial collisions | 1 / 25 |
| Exact metric-signature collisions | 1 / 25 |
| Known-family members detected syntactically and isomorphically | 1 / 25 |
| Known indexed product, power, or path-product closures | 0 / 25 |
| Polynomials having a single path factor but an unindexed residual | 15 / 25 |
| Structurally new grammars | 0 / 25 |
| Witnesses beating the available same-order pressure comparator | 5 / 25 |
| Witnesses extending the order-at-most pressure Pareto frontier | 0 / 25 |
| Witnesses beating the size-unrestricted pressure maximum | 0 / 25 |

The collision is

```text
bush(u6p1,u6p1,u6p1)|n=40
```

which is exactly Galvin `(m,t)=(3,6)`, both as an unrooted tree and as an independence polynomial.

The other 24 coefficient tuples are absent from the bounded local index. That is exact-data novelty relative to the committed/reconstructible index, not new-family novelty: all 24 remain ordinary mixed spider bouquets.

Fifteen of the 25 polynomials have an exact factor `I(P_1;x)=1+x`, `I(P_2;x)=1+2x`, or `I(P_3;x)=1+3x+x^2`. In every case the residual exact coefficient tuple is absent from the 6,619-record index. Thus these are algebraic path factorizations, but not collisions with a known local path-product closure. No indexed product or power collision was found either.

Every public Blair witness has exactly one log-concavity defect. The largest exact post-descent adjacent ratio among them is

```text
454170198467 / 536024379389
= 0.8472939215650911...
```

at `bush(u4p2,u4p3,u6p1,u6p1)|n=57`. The maximum in the reconstructed local single-tree family grids is

```text
0.9808181929611438...
```

at Galvin `(m,t)=(21,11)`. Existing local forest/product audits are tighter still: the committed product-power audit reaches `0.9988814019` without a non-unimodal row. Blair's large normalized log-concavity severities are deep-tail effects; they are not competitive on the near-mode valley-pressure metric driving the current proof program.

After the 41 literature root-stress trees are included, the full comparison index reaches `0.9813215689727014...` at `jp-exact-k0-h8`. Thus the broader committed hard-family corpus strengthens, rather than changes, the frontier comparison.

Order must nevertheless be controlled: adjacent ratios tend to increase with size, so the size-unrestricted maximum alone is not a fair local comparison. Five Blair witnesses set new records relative to the available local trees of exactly the same order:

| Order | Blair ratio | Previous same-order ratio |
|---:|---:|---:|
| 52 | 0.831790 | 0.817336 |
| 54 | 0.838807 | 0.834338 |
| 55 | 0.840281 | 0.813880 |
| 58 | 0.844314 | 0.834586 |
| 59 | 0.846028 | 0.834134 |

None is a Pareto improvement when compared with every indexed tree of order at most its own: the order-27 retained near miss already reaches `0.8571425...`, and the order-57 saved evolutionary witness reaches `0.9241694...`. The right conclusion is therefore “useful same-order frontier data inside a known grammar,” not “no numerical improvement” and not “new mechanism.”

## What is and is not new

The correct classification is:

- **new public exact examples:** yes, for 24 of the 25 published bushes relative to the bounded local coefficient/tree index;
- **new construction grammar:** no;
- **new known-family coverage:** no; the public bank includes a Galvin duplicate, and both standalone seeds are existing Li-family members;
- **new same-order records:** yes, five within the bounded index;
- **new order-versus-pressure Pareto frontier:** no;
- **reason to launch duplicate enumeration:** no.

The full 4,445-way exact collision count is unresolved because 4,420 coefficient records are not committed. The raw 50,000-row Ramos--Sun Prüfer corpus is likewise no longer local, so an exhaustive pairwise identity audit would require recovering both sides. Neither missing artifact affects the grammar-containment or frontier decision.

## Enumeration gate

The gate is closed. A new computational campaign is justified only if it changes at least one of the following:

1. the rooted-tree DP grammar, rather than widening parameters inside `bouquet_adj`;
2. the near-mode post-descent pressure frontier;
3. the location/coupling of log-concavity defects relative to that pressure;
4. a proof-relevant root or conditional-law state not already represented.

Blair V2 supplies better same-order parameter points but does not change the grammar or the order-versus-pressure Pareto frontier. Do not run another broad bouquet scan. Preserve the five same-order records as comparator data. If the unpublished 4,420 records become available, ingest them with the existing script and report collision counts; do not regenerate the same local search space under a new name.

## Reproduction

With the pinned upstream checkout available locally:

```bash
python3 scripts/compare_blair_rooted_bush.py \
  --blair-repo /private/tmp/willblair-verified-combinatorics-20260716 \
  --output results/blair_rooted_bush_comparison_20260716.json
```

The generated JSON contains every public coefficient vector, canonical graph6 certificate, exact rational metric, collision source, validation record, and limitation.
