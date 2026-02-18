# Injection Exploration for Erdős #993

## Date: 2026-02-14

## Summary

We explored three untouched corners of Stanley's taxonomy of unimodality proof
techniques: direct injection, chain decomposition, and sign-reversing involution.
The injection exploration yielded a striking new empirical result.

## Key Finding: Augmented Hall Property

### The Obstruction (Containment-Only)

The simple "add one vertex" injection (containment bipartite graph between
level k and level k+1) **fails** for a large fraction of trees:

| n  | Trees with matching failure | Cause |
|----|---------------------------|-------|
| 5  | 1/3 (33%)    | K_{1,4}: {center} is maximal IS of size 1 |
| 8  | 11/23 (48%)  | Stars, near-stars: small dominating IS |
| 10 | 21/106 (20%) | Trees with high-degree vertices |
| 14 | 1718/3159 (54%) | Grows with n |

The failures occur because **maximal independent sets of size < mode** have
no containment neighbors in the next level. They can't be extended by adding
a vertex.

All 6 canonical injection rules (smallest, largest, near_root, far_root,
min_degree, max_degree) give identical success rates (~27%), confirming the
failure is structural, not a consequence of the choice function.

### The Resolution (Containment + Swap)

Augmenting the bipartite graph with **swap edges** resolves the obstruction.
A swap edge connects S ∈ level_k to T ∈ level_{k+1} when
T = (S \ {u}) ∪ {v, w} for some u ∈ S, v, w ∉ S, with T independent.

| n  | Trees | Levels checked | Augmented matching failures |
|----|-------|---------------|---------------------------|
| 3-4 | 3   | 3             | 0 |
| 5-9 | 87  | 230           | 0 |
| 10 | 106  | 338           | 0 |
| 11 | 235  | 854           | 0 |
| 12 | 551  | 2,238         | 0 |
| 13 | 1,301 | 5,448        | 0 |
| 14 | 3,159 | 14,416       | 0 |
| 15 | 7,741 | ~37,000      | 0 |
| 16 | 19,320 | ~100,000    | 0 |
| **Total** | **32,506** | **~160,000** | **0** |

**Zero failures in ~60,000 level checks across 13,000+ trees.**

Additionally: **zero levels with min-degree 0 in the augmented graph**.
Every independent set below the mode has at least one augmented neighbor
in the next level.

### Swap Edges Are Abundant

The swap mechanism provides massive connectivity:

| n  | Containment edges | Swap edges | Ratio |
|----|------------------|------------|-------|
| 10 | 34,435    | 177,804    | 5.2x  |
| 12 | 678,271   | 5,191,836  | 7.7x  |
| 14 | 11,910,081 | 126,437,523 | 10.6x |

Swap edges grow much faster than containment edges with n.

## Conjecture

**Augmented Hall Property:** For any tree T and any k < mode(I(T)),
the augmented bipartite graph between level_k and level_{k+1}
(containment + swap edges) satisfies Hall's condition.

**Corollary:** The independence polynomial of every tree is unimodal.

Proof sketch (assuming AHP):
- Ascending (k < mode): AHP gives injection level_k → level_{k+1}, so i_k ≤ i_{k+1}.
- Descending (k > mode): Reverse containment (remove vertex) has left degree k+1,
  making Hall's condition easy.

## Why This Approach Is New

In Stanley's taxonomy of unimodality proof techniques, this is a
**non-containment injection**: the map is not S ↦ S ∪ {v} but sometimes
S ↦ (S \ {u}) ∪ {v, w}. This "replacement injection" mechanism appears
not to have been explored for independence polynomials of trees.

The key insight: maximal IS of small size (the obstruction to containment
injection) have excellent swap properties because removing any vertex u
"frees" all of u's tree neighbors, providing abundant candidates for the
two replacement vertices.

## Toward a Proof: The Private Neighbor Bound

### Theorem (Private Neighbor Bound)

For any tree T on n vertices and any dominating independent set S of
size k, the total number of private neighbors is P ≥ n - 2k + 1.

**Proof.** S is independent, so all edges from S go to V\S. S is dominating,
so each v ∈ V\S has ≥1 S-neighbor. A vertex v ∈ V\S is *private for u*
if u is v's only S-neighbor. Let P = number of private vertices,
Q = n - k - P = number of shared vertices (≥2 S-neighbors each).
Edge count: edges from S to V\S ≥ P·1 + Q·2 = P + 2(n-k-P) = 2(n-k) - P.
But edges from S to V\S ≤ n - 1 (total edges in tree).
So 2(n-k) - P ≤ n - 1, giving P ≥ n - 2k + 1. ∎

The bound is tight: achieved at K_{1,n-1} with S = {center}, k=1, P=n-1.

### Lemma (Tree Swap)

In a tree, if u ∈ S has two private neighbors v, w, then v and w are
non-adjacent (since a triangle u-v-w is impossible in a tree). Therefore
(S \ {u}) ∪ {v, w} is an independent set of size k+1: a valid swap.

### Pigeonhole Route (works for most cases)

If P ≥ n - 2k + 1 > k (i.e., n > 3k - 1, i.e., n ≥ 3k), then by
pigeonhole some u ∈ S has priv(u) ≥ ⌈P/k⌉ ≥ 2.

**Verified for 99.6% of cases** (n ≤ 18). The remaining 0.4% have n < 3k.

### The n < 3k Violations: Double-Star Structure

All 20 violation cases (n = 14, 16) share the same structure:

- Tree is a **double star** S(a,b): two hubs connected by a path, each
  with many leaves
- The IS consists of **one hub vertex** (degree 7-10) plus **several leaves
  from the other hub** (degree 1 each)
- The hub IS vertex alone provides 7-9 private neighbors (all its
  leaf-neighbors)
- **Actual P = 7-9, far exceeding the bound P ≥ 5**
- max_priv ≥ 7 in all cases (swap trivially available)

The bound P ≥ n-2k+1 is loose because it treats all edges from S equally;
in reality, a single high-degree IS vertex concentrates most private
neighbors.

**Key observation**: In all 13,042 cases (including all 20 violations),
max_priv ≥ 2. This is strictly stronger than the pigeonhole argument.

### The Leaf-Neighbor Property (CORRECTED: fails at n=18)

**Original claim (n ≤ 16, 13,042 cases, 0 failures):** For any tree T
and any maximal IS S of size k < mode(I(T)), some vertex u ∈ S has ≥ 2
leaf-neighbors in the tree.

**Status: FALSE.** Fails at n=18 for a spider tree (center degree 5,
four arms of length 4, one pendant leaf). The tree has arms:
17-4-13-9-0, 17-5-14-10-1, 17-6-15-11-2, 17-7-16-12-3, 17-8.
IS S = {0, 1, 2, 3, 17} (center + four distant feet), k=5, mode=6.
Center vertex 17 has only 1 leaf-neighbor (vertex 8, the pendant).
But priv(17) = 5: all 5 neighbors are private because the "next" nodes
down each arm (13, 14, 15, 16) are not in S.

Discovered by Gemini (2026-02-14). Also see `notes/leaf_neighbor_proof.md`.

### The Private Neighbor Property (correct statement)

**Conjecture (PNP, n ≤ 18, 99,147 cases, 0 failures):** For any tree T
and any maximal IS S of size k < mode(I(T)), some u ∈ S has priv(u) ≥ 2.

Why this suffices for swap existence:
- Two private neighbors v, w of u ∈ S are non-adjacent in a tree
  (acyclicity prevents triangle u-v-w).
- Therefore (S \ {u}) ∪ {v, w} is independent of size k+1: valid swap. ∎

For the n=18 counterexample: max_priv = 5 (vertex 17 has 5 private
neighbors). The swap is trivially available via private (non-leaf)
neighbors.

PNP is proved for 99.6% of cases via pigeonhole (when n ≥ 3k).
The remaining 0.4% (all n < 3k): if PNP fails, then priv(u) ≤ 1 for
all u, giving P ≤ k, hence k ≥ (n+1)/3. For trees with mode ≤ ⌊n/3⌋+1
this already gives k ≥ mode. For high-mode trees, no such IS exists
(verified exhaustively through n = 18; 99,147 IS, 0 failures).

**Note on Bollobás-Cockayne:** The claim priv(u) ≥ 1 (external) for all
u in a maximal IS is FALSE in general (counterexample: K_{1,4}, S = all
leaves, priv = 0 for every leaf). The proof does NOT require this: it
only needs P ≤ k when all priv ≤ 1 (PNP failure). Verified: 34,055
maximal IS below mode at n=18 have some vertex with priv = 0, yet
max_priv ≥ 2 in all cases.

### The Remaining Gaps

Two theoretical questions remain open:

1. **PNP for n < 3k**: Prove that every maximal IS of size k < mode(I(T))
   has some u with priv(u) ≥ 2, even when n < 3k.
   - 405 empirical cases, 0 failures
   - All are double-star-like trees with one high-degree hub in S
   - Likely provable via degree concentration arguments

2. **From individual degree to Hall's condition**: Even proving every IS
   has augmented degree ≥ 1 is not enough for a full matching. Hall's
   condition requires |N(A)| ≥ |A| for ALL subsets A ⊆ level_k.
   - Need expansion properties of the augmented bipartite graph
   - Augmented min-degree ≥ 3 empirically (from swap degree analysis)
   - Swap edges grow as ~10x containment edges, suggesting rich expansion

## Proof Strategy Summary

To prove the Augmented Hall Property, we need to show: for any subset
A ⊆ level_k, |N_aug(A)| ≥ |A|.

Possible approaches:
1. **Vertex expansion in trees**: use tree sparsity (n-1 edges) to bound
   expansion of swap neighborhoods.
2. **Separate cases**: non-maximal IS (containment suffices) vs maximal IS
   (swap suffices), then combine.
3. **Double counting**: count edges in the augmented bipartite graph from
   both sides; show average degree is high enough.
4. **Private Neighbor Bound + pigeonhole** (new): covers 99.85% of cases.
   The 0.15% double-star exceptions need a separate structural argument.

## Sign-Reversing Involution (SRI)

### Motivation

The augmented matching (Hopcroft-Karp) proves a matching EXISTS but doesn't
construct one canonically. A sign-reversing involution gives a CANONICAL
pairing phi: IS_k u IS_{k+1} -> IS_k u IS_{k+1} where phi^2 = id,
phi swaps sizes (sign-reversing), and fixed points are all at size k+1.
This directly proves i_{k+1} >= i_k.

### Strategies Tested (n=5..12, 982 trees, 3,679 levels)

| Strategy | Perfect Levels | % |
|----------|---------------|---|
| toggle/label | 0 | 0.0% |
| toggle/rev_label | 0 | 0.0% |
| toggle/deg_asc | 0 | 0.0% |
| toggle/deg_desc | 0 | 0.0% |
| toggle/bfs | 0 | 0.0% |
| toggle/dfs_post | 0 | 0.0% |
| leaf_swap alone | 0 | 0.0% |
| contain_first/label | 3,679 (weak test) | see correction below |
| contain_first/deg_asc | 3,679 (weak test) | see correction below |
| hybrid/label | 982 | 26.7% |

### Containment-First Involution (CORRECTED)

**Definition.** For k < mode(I(T)), define phi on IS_k u IS_{k+1}:

**Forward (size k -> size k+1):**
- If S is **non-maximal**: phi(S) = S u {v}, where v is the first vertex
  (label order) not in S with no neighbor in S.
- If S is **maximal**: phi(S) = (S \ {u}) u {v, w}, where (u, v, w) is the
  canonical swap triple (leaf-swap or private-neighbor swap).

**The collision problem.** The original test (`explore_sri.py`) checked each
forward-reverse pair (S -> T -> S') individually using the known forward rule
type, but did NOT check whether multiple left elements map to the same right
element. A stronger test (`diagnose_sri_collision.py`) reveals **massive
collisions**: containment and swap forward maps frequently target the same T.

| n | Levels | Collision levels | % |
|---|--------|-----------------|---|
| 5 | 6 | 3 | 50.0% |
| 8 | 64 | 41 | 64.1% |
| 10 | 338 | 232 | 68.6% |
| 12 | 2,238 | 1,687 | 75.4% |
| 13 | 5,448 | 4,147 | 76.1% |

**Example (K_{1,4}, k=1):** T = {leaf1, leaf2} is simultaneously:
- Containment image of {leaf1} (adding leaf2 as first free vertex)
- Containment image of {leaf2} (adding leaf1 as first free vertex)
- Swap image of {center} (removing center, adding leaf1 + leaf2)

Three left elements target the same right element, so the forward map is
not injective and phi is not a valid involution.

**Conclusion:** The containment-first SRI as designed does NOT produce a
valid sign-reversing involution. The original "100% success through n=16"
claim was an artifact of an insufficient correctness test. Designing a
collision-free SRI for tree independence polynomials remains open.

### Why Collisions Are Fundamental

The collision arises because:
1. Multiple non-maximal IS can share the same "first free vertex," so
   their containment images coincide.
2. A swap image T = (S \ {u}) u {v, w} can also be the containment image
   of a non-maximal IS S' = T \ {x} for some x in T.

Resolving collisions would require a modified forward map that redirects
some IS to alternative targets when the canonical target is already claimed.
This requires global coordination (knowing all forward targets) rather than
a local per-element rule, which is antithetical to the SRI framework.

### What DOES Work

The **augmented bipartite matching** (Hopcroft-Karp algorithm) correctly
finds an injection IS_k -> IS_{k+1} for every k < mode, verified through
n=18 across 204,909 trees. This non-constructive existence proof is valid,
even though we lack a canonical/constructive injection.

### Open Question

Can a valid SRI be designed for tree independence polynomials? This would
require either:
1. A collision-free forward rule (seems difficult given the structural
   reasons for collisions)
2. A fundamentally different involution design (e.g., operating on a
   different combinatorial object)
3. An SRI on a related polynomial (e.g., matching polynomial) with a
   transfer argument

## Implications for the Taxonomy Mapping

| Stanley category | Status in this project | New insight |
|-----------------|----------------------|-------------|
| Generating function | Heavily explored (approaches 2,4,6,8) | -- |
| Real-rootedness | Blocked | -- |
| Log-concavity | Blocked (Galvin 2025) | -- |
| **Direct injection** | **Now explored** | **Containment fails; augmented works** |
| sl_2 representations | Not explored | -- |
| **Sign-reversing involution** | **Explored; collision problem** | **Naive SRI fails; valid SRI open** |
| Chain decomposition | Not explored | -- |

## Files

- `explore_injection.py`: Containment-only matching + canonical rules
- `explore_augmented_injection.py`: Containment + swap matching
- `explore_why_augmented_works.py`: Structural analysis of swap mechanism
- `explore_sri.py`: Sign-reversing involution (4 strategies, 6 orderings)
- `verify_mode_bound.py`: Checks mode > (n-1)/3 and P bound tightness
- `verify_n3k_claim.py`: Detailed n >= 3k analysis
- `analyze_violations.py`: Tree structure of n < 3k cases
- `explore_sri_priv.py`: SRI with private-neighbor swaps (strong test)
- `diagnose_sri_collision.py`: Collision diagnosis for containment-first SRI
- `verify_pnp_extended.py`: PNP verification through n=18
- `results/injection_exploration.json`: Containment results data
- `notes/pnp_n_lt_3k_proof.md`: 1-Private Conjecture (Gemini analysis)
