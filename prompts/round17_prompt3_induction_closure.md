# Round 17, Prompt 3: Full Induction Structure and STP2(I,E) Derivation

**Target model:** GPT 5.2 Pro (theory-heavy)

---

## Context

We are proving the independence set sequence of every tree is unimodal via a support-vertex DP induction. The core mechanism is the Cauchy-Binet expansion of the **ladder minor** Λ_k that encodes STP2.

Round 16 established:
- The CB expansion Λ_k^{new} = Σ_{i,j} Δ_{i,j}(A,B)·P(k-i)Q(k-j)
- Ladder-minor monotonicity (5.5M checks, 0 decreases)
- Condition C is DEAD (θ diverges)
- The induction step 6 (deriving STP2 from P2) is broken as stated

The proof now needs a clean induction that doesn't rely on Condition C or on deriving STP2 from E ≽ J.

## The Induction Target

**Conjecture (Chain STP2)**: For every tree T and every rooting at vertex v, the triple (I_v, E_v, J_v) satisfies:
1. STP2(E_v, J_v): ladder minor Λ_k(E,J) ≥ 0 for all k
2. STP2(I_v, E_v): ladder minor Λ_k(I,E) ≥ 0 for all k

**Verified**: 0 failures across all rootings of all trees n ≤ 22.

**Consequence at support vertices**: Chain STP2 + the W-form identity gives w_k ≥ 0, which gives E ≽ J prefix, which gives unimodality of I = E + xJ.

## Task 1: Structure the induction precisely

The DP processes children c_1, ..., c_s of root r (non-leaf children only; leaves contribute (1+x)^ℓ to E_acc):

**Initial state** (after ℓ leaf children):
```
E_acc = (1+x)^ℓ,  J_acc = 1
```
Both trivially satisfy STP2(E_acc, J_acc) (Λ_k = C(ℓ,k)·0 - C(ℓ,k-1)·0 ≥ 0 for k ≥ 2; need to check k=1).

Actually: E_acc(1) = ℓ, E_acc(0) = 1, J_acc(0) = 1, J_acc(1) = 0.
Λ_1 = E_acc(1)·J_acc(1) - E_acc(0)·J_acc(2) = 0 - 0 = 0 ≥ 0. ✓
(For J_acc = [1], all Λ_k = 0 for k ≥ 1.)

**Step t**: process non-leaf child c_t with subtree triple (I_t, E_t, J_t).
```
E_new = E_acc · I_t,  J_new = J_acc · E_t
```

**Inductive hypothesis**: (I_t, E_t, J_t) satisfies chain STP2 (by induction on subtree size).

**To show**: STP2(E_new, J_new) holds, i.e., Λ_k^{new} ≥ 0 for all k.

**After all children**: the final (E, J) pair satisfies STP2(E, J).

**Then show**: STP2(I, E) where I = E + xJ.

**State the complete induction claim** with all hypotheses and what needs to be verified at each step. Be precise about:
- Base case (single vertex, path of length 2, star)
- Which properties are hypothesized vs which are consequences
- Whether we need STP2(I,E) at EVERY rooting or just support vertices
- The role of LC(E), LC(J), J ≤ E (are these hypothesized or proved independently?)

## Task 2: STP2(I,E) from STP2(E,J)

The STP2(I,E) diagonal condition (under LC(I), which holds since I is a product of LC factors):
```
E(k+1)·I(k-1) ≤ E(k)·I(k)
```

Using I(k) = E(k) + J(k-1):
```
E(k+1)[E(k-1) + J(k-2)] ≤ E(k)[E(k) + J(k-1)]
```
Rearranging:
```
c_k(E) + [E(k)J(k-1) - E(k+1)J(k-2)] ≥ 0
```
where c_k(E) = E(k)² - E(k-1)E(k+1) ≥ 0 is the LC gap.

**The correction term** E(k)J(k-1) - E(k+1)J(k-2) **can be negative**. Example: star S_n, J = [1], k = 2: gives -E(3) = -C(n,3) < 0.

But c_k(E) rescues. At the star: c_2((1+x)^n) = C(n,2)² - nC(n,3), which is ≥ C(n,3) for n ≥ 3. ✓

**Key question**: can STP2(I,E) be derived from STP2(E,J) + LC(E) + J ≤ E?

### Sub-question 2a: What does STP2(E,J) give about E(k)J(k-1) - E(k+1)J(k-2)?

STP2(E,J) diagonal: E(m)J(m) ≥ E(m-1)J(m+1).

At m = k-1: E(k-1)J(k-1) ≥ E(k-2)J(k). Does this help bound E(k)J(k-1) - E(k+1)J(k-2)?

At m = k: E(k)J(k) ≥ E(k-1)J(k+1). This gives E(k)/E(k-1) ≥ J(k+1)/J(k) (when J(k) > 0).

We need E(k)J(k-1) ≥ E(k+1)J(k-2), i.e., E(k)/E(k+1) ≥ J(k-2)/J(k-1).

By LC(E): E(k)/E(k+1) ≥ E(k+1)/E(k+2) (nondecreasing reciprocal ratios).
By LC(J): J(k-2)/J(k-1) ≤ J(k-3)/J(k-2) (nonincreasing J ratios ⟹ nondecreasing J^{-1} ratios).

Can STP2(E,J) bridge the gap? STP2(E,J) at m = k-1, n = k-2: J(k)/J(k-1) ≤ E(k-1)/E(k-2).

**Work out**: does c_k(E) + E(k)J(k-1) - E(k+1)J(k-2) factor nicely using STP2(E,J)?

### Sub-question 2b: Alternative approach via direct CB

Instead of deriving STP2(I,E) from STP2(E,J), can we prove STP2(I,E) directly via its own CB expansion?

At a step where E_acc → E_acc · I_t and we want STP2(I_new, E_new):
```
I_new = E_new + x·J_new = E_acc·I_t + x·J_acc·E_t
```

The ladder minor Λ_k(I_new, E_new) has its own CB structure. Since I_new involves BOTH E_acc and J_acc (through I_t = E_t + xJ_t), this might be more complex.

**Alternatively**: if STP2(I,E) only needs to hold at the FINAL state (after all children), and if it follows from STP2(E,J) + LC(E) + some condition that only holds at the end, it might not need step-by-step preservation.

### Sub-question 2c: Does STP2(I,E) follow from STP2(E,J) at support vertices?

At support vertices, we additionally have E ≽ J (prefix ratio dominance, from w_k ≥ 0).

E ≽ J prefix: J(k)/E(k) ≤ J(k-1)/E(k-1) for k ≤ mode. This gives E(k)J(k-1) ≥ E(k-1)J(k).

**Check**: does E ≽ J imply E(k)J(k-1) ≥ E(k+1)J(k-2)?

E ≽ J gives E(k)J(k-1) ≥ E(k-1)J(k). By LC(E): E(k-1)/E(k) ≥ E(k)/E(k+1), so E(k-1) ≥ E(k)²/E(k+1). Thus:
```
E(k)J(k-1) ≥ E(k-1)J(k) ≥ [E(k)²/E(k+1)]J(k)
```
We need E(k)J(k-1) ≥ E(k+1)J(k-2), i.e., [E(k)²/E(k+1)]J(k) ≥ E(k+1)J(k-2)?
That gives E(k)²J(k) ≥ E(k+1)²J(k-2), i.e., [E(k)/E(k+1)]² ≥ J(k-2)/J(k).

By LC(J): J(k-2)/J(k) ≤ [J(k-1)/J(k)]². And by E ≽ J: J(k)/E(k) ≤ J(k-1)/E(k-1), so J(k-1)/J(k) ≥ E(k-1)/E(k). Thus [J(k-1)/J(k)]² ≥ [E(k-1)/E(k)]².

So J(k-2)/J(k) ≤ [E(k-1)/E(k)]². And we need [E(k)/E(k+1)]² ≥ J(k-2)/J(k) ≤ [E(k-1)/E(k)]².

By LC(E): E(k)/E(k+1) ≥ E(k-1)/E(k). So [E(k)/E(k+1)]² ≥ [E(k-1)/E(k)]² ≥ J(k-2)/J(k). ✓

**If this chain is correct, STP2(I,E) follows from E ≽ J + LC(E) + LC(J) at support vertices!**

**Verify this carefully.** The argument chains three inequalities and each step needs to be checked for edge cases (J(k) = 0, indices near 0, prefix vs full range).

## Task 3: Non-support vertex STP2(I,E)

At non-support vertices, E ≽ J may fail. We still need STP2(I,E) for the parent's CB expansion to use.

**Option A**: Prove STP2(I,E) directly via CB, analogous to STP2(E,J).

**Option B**: Show STP2(I,E) follows from STP2(E,J) + LC(E) + c_k(E) bounds, WITHOUT E ≽ J.

For option B, the condition is c_k(E) + E(k)J(k-1) - E(k+1)J(k-2) ≥ 0.

This can be rewritten as E(k)I(k) - E(k+1)I(k-1) = E(k)I(k) - E(k+1)I(k-1).

Using I = E + xJ:
= E(k)[E(k)+J(k-1)] - E(k+1)[E(k-1)+J(k-2)]
= c_k(E) + [E(k)J(k-1) - E(k+1)J(k-2)]

**The correction** E(k)J(k-1) - E(k+1)J(k-2) resembles a "shifted" version of the STP2(E,J) diagonal. Specifically:

STP2(E,J) at m=k-1: E(k-1)J(k-1) ≥ E(k-2)J(k).

**Can we relate** E(k)J(k-1) - E(k+1)J(k-2) to STP2(E,J) minors via an identity?

Try: E(k)J(k-1) = E(k)/E(k-1) · E(k-1)J(k-1) ≥ E(k)/E(k-1) · E(k-2)J(k) [by STP2 at m=k-1].

Then: E(k)J(k-1) - E(k+1)J(k-2) ≥ [E(k)/E(k-1)]·E(k-2)J(k) - E(k+1)J(k-2).

Does this lead somewhere useful?

## Task 4: The minimal induction package

Based on the analysis above, state the **minimal set of properties** that need to be proved at each step:

1. Λ_k(E_new, J_new) ≥ 0 for all k [STP2(E,J) preservation — via CB]
2. ??? [STP2(I,E) — derived or proved in parallel?]
3. LC(E_new), LC(J_new) [automatic from products]
4. J_new ≤ E_new [proved independently]

Is property 2 free from property 1 + LC + J≤E? Or does it require its own CB argument?

If STP2(I,E) follows from STP2(E,J) + LC (as sketched in Task 2c for support vertices), then the induction is:
- **Induct on**: STP2(E,J) at every rooting
- **Free consequences**: LC, J ≤ E, STP2(I,E)
- **The single thing to prove**: Λ_k^{new} ≥ 0 via CB (i.e., X_k ≥ 0)

**State this cleanly.** If some step doesn't work at non-support vertices, identify what additional property is needed there.

## Task 5: The base cases

Verify the induction base cases:

1. **Single vertex** (n=1): E = 1, J = 1. I = 1+x. STP2(E,J): trivial (length 1). STP2(I,E): I(1)E(0) ≥ I(0)E(1) ⟹ 1·1 ≥ 1·0 ✓.

2. **Edge** (n=2): Root at one vertex, one child (a leaf).
   E = 1+x, J = 1. I = 1+2x. STP2(E,J): Λ_1 = E(1)J(1) - E(0)J(2) = 0 ≥ 0. ✓
   STP2(I,E): Λ_1 = I(1)E(0) - I(0)E(1) = 2·1 - 1·1 = 1 ≥ 0. ✓

3. **Star** S_n (root = center): E = (1+x)^n, J = [1]. All Λ_k(E,J) = 0 for k ≥ 1 (J is constant). STP2(I,E): c_k(E) - E(k+1)·0 + E(k)·0 + ... = c_k(E) ≥ 0. ✓ (for k ≥ 2; need separate check for k=1).

4. **Path** P_n (root at support vertex, one leaf child): STP2 verified computationally.

**Verify these algebraically** and confirm no edge-case issues.

## Deliverables

1. Complete induction statement with all hypotheses, base cases, and inductive step
2. STP2(I,E) derivation from STP2(E,J) — does the chain argument in Task 2c work?
3. Non-support vertex STP2(I,E): does it follow from STP2(E,J) + LC, or need separate argument?
4. The minimal induction package: what's proved at each step, what's free?
5. Base case verification
6. Any logical gaps or circular dependencies in the proposed induction
7. Overall assessment: is the induction now self-contained (modulo X_k ≥ 0)?
