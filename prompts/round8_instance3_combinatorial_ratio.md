# Round 8, Instance 3: Why is J_k/E_k nonincreasing? Combinatorial and structural approaches

## The Problem

We're proving that the independence polynomial of every tree is unimodal (Erdős Problem #993). The full problem reduces to showing ratio dominance E ≽ J at every support vertex.

**E ≽ J means:** E_{k+1}·J_k ≥ E_k·J_{k+1} for ALL k ≥ 0.

Equivalently: the ratio J_k/E_k is nonincreasing in k (for k where E_k > 0).

## Combinatorial meaning

At a support vertex r of a tree T:
- E_k = #{independent sets of size k that exclude r}
- J_k = #{independent sets of size k+1 that include r}

So J_k/E_k = #{IS of size k+1 containing r} / #{IS of size k excluding r}.

**Why should this ratio decrease in k?** As k grows, it becomes "harder" to include r (because including r forces all neighbors excluded, using up one vertex and blocking deg(r) others). The ratio measures the "ease" of including r at size k+1 versus not including r at size k.

Verified: 907M+ checks across 9.1M trees n ≤ 22, 0 failures.

## Product structure at support vertices

Root at support vertex r with ℓ leaf neighbors and non-leaf subtrees T_1,...,T_s.

```
E = (1+x)^ℓ · ∏ I(T_j)     (each subtree IS poly contributes independently)
J = ∏ E(T_j)                (including r forces all subtree roots excluded)
```

Incrementally:
```
Stage 0: E^{(0)} = (1+x)^ℓ,  J^{(0)} = [1]
Stage t: E^{(t)} = E^{(t-1)} · I_t,  J^{(t)} = J^{(t-1)} · E_t
```

## What's proved and what's not

| Property | Status |
|----------|--------|
| E is PF2 (nonneg LC) | **PROVED** (products of PF2) |
| J ≤ E coefficientwise | **PROVED** |
| Karlin main part ≥ 0 | **PROVED** (E^{(t-1)}·E_t ≽ J^{(t-1)}·E_t) |
| P3: e_k ≥ j_{k-1} | **PROVED** (leaf-swap injection) |
| E ≽ J | VERIFIED (0 fails), need proof |
| SCC: (1+x)I ≽ E | VERIFIED (0 fails), need proof |

## Approach 1: Diagonal convolution matrix

The incremental step acts on (E, J) by:
```
[E_new]   [I_t   0 ] [E_old]
[J_new] = [0    E_t] [J_old]
```
in convolution sense. This is a DIAGONAL matrix M_t = diag(I_t, E_t).

**Question 1:** Is there a total-positivity condition on M_t that preserves E ≽ J?

If we could show that M_t preserves ≽ (ratio dominance), the proof would close by induction. But the standard condition would be I_t ≽ E_t, which FAILS at ~30% of factors.

However, E^{(t)} ≽ J^{(t)} ALWAYS holds despite factor-level failures. The product of multiple factors rescues the property. Why?

**Observation:** At s=1 (one non-leaf child), E ≽ J follows from SCC of the subtree via Karlin + transitivity. This handles 63% of support vertices. At s ≥ 2, the product of two or more "imperfect" factors somehow preserves ratio dominance.

## Approach 2: Combinatorial injection

E ≽ J at index k says: #{IS of size k+1 excluding r} · #{IS of size k+1 including r at index k} ≥ #{IS of size k excluding r} · #{IS of size k+2 including r}.

This is a "cross ratio" comparison. Can we construct a combinatorial injection or coupling?

**Idea:** The leaf-swap injection (for P3) maps S ↦ (S\{r}) ∪ {leaf}. This maps "include r, size k+1" injectively to "exclude r, size k+1", proving e_{k+1} ≥ j_k. But E ≽ J requires a MULTIPLICATIVE comparison, not additive.

**Question 2:** Is there a measure-preserving coupling between the product sets that witnesses E ≽ J?

Concretely, define:
- A_k = {IS of size k excluding r}, |A_k| = e_k
- B_k = {IS of size k+1 including r}, |B_k| = j_k

E ≽ J says: |A_{k+1}| · |B_k| ≥ |A_k| · |B_{k+1}|.

Equivalently: there is an injection φ: A_k × B_{k+1} → A_{k+1} × B_k.

Can you construct φ? Each element of A_k × B_{k+1} is a pair (S, T) where S is an IS of size k not containing r, and T is an IS of size k+2 containing r.

## Approach 3: The s=1 → s=2 lifting

**s=1 is PROVED.** When r has one non-leaf child c and ℓ leaves:
```
E = (1+x)^ℓ · I_c,  J = E_c
E ≽ J follows from: SCC at c → Karlin → transitivity
```

**s=2 is the key battleground (30% of vertices).** When r has two non-leaf children c₁, c₂ and ℓ leaves:
```
E = (1+x)^ℓ · I₁ · I₂,  J = E₁ · E₂
```

After stage 1: E^{(1)} = (1+x)^ℓ · I₁ ≽ E₁ = J^{(1)} (by s=1 proof).

At stage 2: need (1+x)^ℓ · I₁ · I₂ ≽ E₁ · E₂.

By Karlin: E^{(1)} ≽ J^{(1)} and E₂ PF2 → E^{(1)}·E₂ ≽ J^{(1)}·E₂ = E₁·E₂ = J^{(2)}.

But E^{(2)} = E^{(1)}·I₂ = E^{(1)}·E₂ + x·E^{(1)}·J₂.

So Δ_k(E^{(2)}, J^{(2)}) = Δ_k(E^{(1)}·E₂, J^{(2)}) + Δ_k(x·E^{(1)}·J₂, J^{(2)}).

First term ≥ 0 (Karlin). Second term is the correction from x·J₂.

**Question 3:** At s=2, can you use SCC at c₂ (i.e., (1+x)I₂ ≽ E₂) together with the s=1 ratio dominance E^{(1)} ≽ J^{(1)} to control the correction?

The SCC at c₂ says I₂ has a specific structural relationship to E₂. Since I₂ = E₂ + x·J₂, and SCC gives (1+x)(E₂+xJ₂) ≽ E₂, this constrains J₂ relative to E₂. Maybe this constraint makes the correction small enough?

## Approach 4: Reciprocal sequence / Stieltjes continued fraction

If E is PF2 and J ≤ E, then the ratio r_k = J_k/E_k is a sequence in [0,1] with r_0 = 1 (since J_0 = E_0 = 1 at support vertices). E ≽ J says r_k is nonincreasing.

**Question 4:** Can you characterize when r_k is nonincreasing for product structures?

For (1+x)^ℓ:  if J = [1], then r_k = 1/C(ℓ,k), which is strictly decreasing. ✓

When we multiply by factor (I_c, E_c):
```
r_k^{new} = (E_c * J_old)_k / (I_c * E_old)_k
```

This is a RATIO of convolutions. The ratio r_k^{new} is a weighted average of products of old ratios and factor contributions. Under what conditions on the factor is this ratio still nonincreasing?

## Key constraints (available for any approach)

1. E_t is PF2 (nonneg, LC) at every factor — PROVED
2. J_t ≤ E_t coefficientwise — PROVED
3. (1+x)I_t ≽ E_t (SCC of subtree) — VERIFIED, 0 failures
4. E^{(t-1)} is PF2 — PROVED
5. J^{(t-1)} ≤ E^{(t-1)} — PROVED
6. E^{(t-1)} ≽ J^{(t-1)} — inductive hypothesis
7. Products of PF2 are PF2 (Karlin) — standard
8. Karlin: PF2 kernel preserves ≽ order — standard

## Dead ends

- Factor-level E_t ≽ J_t: FAILS ~14% of factors. Cannot argue factorwise.
- Factor-level I_t ≽ E_t: FAILS ~30%. Transitivity chain fails.
- Generic product closure (non-tree pairs): FALSE. Tree-realizability essential.
- HWZZ partial synchronicity: FALSE at n≥12.

## Notation

| Symbol | Definition |
|--------|-----------|
| E ≽ J | E_{k+1}·J_k ≥ E_k·J_{k+1} for all k |
| PF2 | nonneg + log-concave coefficients |
| SCC | (1+x)I ≽ E (Strong Condition C) |
| I_t = E_t + x·J_t | IS polynomial of t-th non-leaf subtree |
| Karlin | PF2 kernel preserves ≽ (likelihood ratio order) |

## Verification data

Star K_{1,4}, root at center: E = (1+x)^4 = [1,4,6,4,1], J = [1].
- r_k = J_k/E_k: r_0 = 1/1 = 1, r_k = 0 for k ≥ 1. Nonincreasing. ✓

Path P_5 = 0-1-2-3-4, root at vertex 1: E = [1,4,4,1], J = [1,2].
- r_0 = 1, r_1 = 2/4 = 0.5, r_2 = 0/4 = 0. Nonincreasing. ✓

Pendant-star n=7, root at support vertex: E = [1,6,11,10,5,1], J = [1,4,6,4,1].
- r_0=1, r_1=4/6=0.667, r_2=6/11=0.545, r_3=4/10=0.4, r_4=1/5=0.2, r_5=0/1=0. ✓
