# Round 8, Instance 2: Curvature rescue — prove T3 ≥ |T1| always

## The Problem

We're proving that the independence polynomial of every tree is unimodal (Erdős Problem #993). The full problem reduces to showing Strong Condition C (SCC) at every support vertex.

## The 3-term identity (PROVED)

At each incremental product stage, after multiplying the t-th factor into the accumulated pair (I_acc, E_acc):

**b_{k-1} · Δ_k = T1 + T2 + T3** where:
- T1 = b_{k-1} · d_k (current LR minor of (I_acc, E_acc), can be negative)
- T2 = b_k · d_{k-1} (memory from adjacent index)
- T3 = a_{k-1} · c_k (LC curvature amplified by a/b ratio, always ≥ 0)

Here a_k = I_acc_k, b_k = E_acc_k, d_k = a_{k+1}b_k - a_kb_{k+1}, c_k = b_k² - b_{k-1}b_{k+1}.

## The key computational finding

Through n = 22 (7.77M trees, 739.7M (stage,k) checks, 19.6M d_k < 0 events):

| Property | Value |
|----------|-------|
| SCC failures | **0** |
| d_k < 0 at k = 1 | **0** (never) |
| d_k < 0 at k ≥ 2 | 19,559,141 events |
| **T3 ≥ \|T1\| (curvature alone)** | **100%** of d_k < 0 events |
| T2 needed for rescue | **0** (never) |
| Min T3/\|T1\| at n=21 | 7.09 |
| Min T3/\|T1\| at n=22 | **4.23** |

**Curvature alone ALWAYS suffices.** T2 (memory) is never needed. The min ratio is decreasing but still comfortable (4.23x at n=22).

## What T3 ≥ |T1| means algebraically

When d_k < 0 (so T1 = b_{k-1} · d_k < 0):

T3 ≥ |T1|  ⟺  a_{k-1} · c_k ≥ b_{k-1} · |d_k|  ⟺  (a_{k-1}/b_{k-1}) · c_k ≥ |d_k|

Since a_{k-1} = I_k ≥ E_k = b_{k-1} (from J ≤ E coefficientwise, so I = E + xJ ≥ E), the ratio a_{k-1}/b_{k-1} ≥ 1. So it suffices to show:

**c_k ≥ |d_k| whenever d_k < 0**

i.e., the LC gap of E at index k is at least as large as the magnitude of the negative LR minor.

## Product structure

At stage t:
```
E_acc = (1+x)^ℓ · I_1 · ... · I_t       (accumulated E)
J_acc = E_1 · ... · E_t                   (accumulated J)
I_acc = E_acc + x·J_acc
```

The LC gap c_k of E_acc is given by Cauchy-Binet:
```
c_k = Σ_{over all factor pairs} (positive contributions from factor curvatures)
```
Products of PF2 polynomials are PF2, so c_k ≥ 0 always. The question is: is c_k large ENOUGH?

The LR minor d_k = a_{k+1}·b_k - a_k·b_{k+1} where a = I_acc = E_acc + x·J_acc:
```
d_k = (b_{k+1} + j_k)·b_k - (b_k + j_{k-1})·b_{k+1}
    = j_k·b_k - j_{k-1}·b_{k+1}
```

So d_k < 0 means j_{k-1}·b_{k+1} > j_k·b_k, i.e., J_acc is "growing faster" than E_acc at index k.

## The tightest cases

All tightest cases at n=22 have:
- k = 10, stage 1 (single factor), s = 1
- T1 = -128,384, T2 = +3,954,608, T3 = +542,640
- T3/|T1| = 4.23

The tightest case is always s=1 (one non-leaf child). At s=1 with ℓ leaves:
- E_acc = (1+x)^ℓ · I_c, J_acc = E_c
- c_k comes from (1+x)^ℓ · I_c (rich product)
- d_k = j_k · b_k - j_{k-1} · b_{k+1} where j = E_c, b = (1+x)^ℓ · I_c

## s=1 reduction (PROVED, 63% of cases)

For s=1, E ≽ J follows from SCC of the subtree via Karlin + transitivity. So T3 ≥ |T1| at s=1 is already proved indirectly. But a DIRECT proof of T3 ≥ |T1| would be stronger — it would give SCC directly via the 3-term identity.

## Your task

### Priority 1: Prove c_k ≥ |d_k| when d_k < 0

Since c_k = b_k² - b_{k-1}·b_{k+1} (LC gap) and d_k = j_k·b_k - j_{k-1}·b_{k+1}, when d_k < 0:

|d_k| = j_{k-1}·b_{k+1} - j_k·b_k

Need: b_k² - b_{k-1}·b_{k+1} ≥ j_{k-1}·b_{k+1} - j_k·b_k

Rearranging: b_k(b_k + j_k) ≥ b_{k+1}(b_{k-1} + j_{k-1})

i.e., b_k · a_k ≥ b_{k+1} · a_{k-1}

This is exactly d_{k-1} ≥ 0!

**Wait — is this circular?** We need d_{k-1} ≥ 0 to prove T3 ≥ |T1| at index k. But d_{k-1} ≥ 0 is exactly E ≽ J at index k-1. So T3 ≥ |T1| at k follows from E ≽ J at k-1.

**Check:** the computational data says d_k < 0 never at k=1 (only k ≥ 2). At k ≥ 2, we need d_{k-1} ≥ 0. If d_{k-1} is also negative, we'd need T3 ≥ |T1| at k-1 too, which needs d_{k-2} ≥ 0...

**Key question:** Is there a forward induction in k?
- d_0 ≥ 0 always (trivially: a_1·b_0 ≥ a_0·b_1 follows from... what?)
- d_1 ≥ 0 always (computationally verified: d_k < 0 never at k=1)
- For k ≥ 2: if d_{k-1} ≥ 0, then T3 ≥ |T1| at k, and with T2 ≥ 0 (from d_{k-1} ≥ 0), we get Δ_k ≥ 0.

**Does this chain work?** Need to verify:
1. d_0 ≥ 0 (base case)
2. d_1 ≥ 0 (seems to always hold)
3. d_{k-1} ≥ 0 ⟹ T3 ≥ |T1| at k (the rearrangement above)
4. d_{k-1} ≥ 0 ⟹ T2 ≥ 0 at k (immediate from definition T2 = b_k · d_{k-1})
5. Steps 3+4 give T1+T2+T3 ≥ 0, i.e., Δ_k ≥ 0 (SCC at k)

**But this proves SCC, not E ≽ J.** The chain proves Δ_k ≥ 0 for all k (SCC) assuming d_0, d_1 ≥ 0. It does NOT prove d_k ≥ 0 for all k (which is E ≽ J).

**Can we extract d_k ≥ 0 from Δ_k ≥ 0?** Yes, if we also know T2 ≥ 0 and T3 ≥ 0:
Δ_k ≥ 0 means b_{k-1}·d_k + b_k·d_{k-1} + a_{k-1}·c_k ≥ 0.
If d_{k-1} ≥ 0 and c_k ≥ 0 (always), then even if d_k < 0, we can't conclude d_k ≥ 0 from Δ_k ≥ 0 alone.

So the CIRCULARITY question is critical. Please investigate:
1. Is c_k ≥ |d_k| equivalent to d_{k-1} ≥ 0?
2. If yes, can we set up a forward induction that proves BOTH d_k ≥ 0 and Δ_k ≥ 0 simultaneously?
3. What is the base case? Why is d_0 ≥ 0 and d_1 ≥ 0?

### Priority 2: The amplification factor

Even if c_k < |d_k|, the full T3 = a_{k-1}·c_k might still exceed |T1| = b_{k-1}·|d_k| because a_{k-1}/b_{k-1} ≥ 1. The amplification factor is:

a_{k-1}/b_{k-1} = (b_{k-1} + j_{k-2})/b_{k-1} = 1 + j_{k-2}/b_{k-1}

This is ≥ 1 and can be significantly larger (especially at small k where J contributes heavily). Can this amplification compensate for c_k < |d_k|?

### Priority 3: Forward induction on d_k

Does d_k satisfy a recursion? From d_k = j_k·b_k - j_{k-1}·b_{k+1}:

d_{k+1} = j_{k+1}·b_{k+1} - j_k·b_{k+2}

Can we relate d_{k+1} to d_k via the coefficients? Using LC of E (c_{k+1} ≥ 0) and J ≤ E:

b_{k+2} ≤ b_{k+1}²/b_k (from LC), so
d_{k+1} ≥ j_{k+1}·b_{k+1} - j_k·b_{k+1}²/b_k = b_{k+1}(j_{k+1} - j_k·b_{k+1}/b_k)

This needs j_{k+1}/j_k ≥ b_{k+1}/b_k, which is NOT guaranteed (it's the opposite of E ≽ J!).

## Notation

| Symbol | Definition |
|--------|-----------|
| a_k = I_k | coefficient of I = E + xJ |
| b_k = E_k | coefficient of E |
| j_k = J_k | coefficient of J |
| d_k = a_{k+1}b_k - a_kb_{k+1} = j_kb_k - j_{k-1}b_{k+1} | LR minor |
| c_k = b_k² - b_{k-1}b_{k+1} | LC gap of E (always ≥ 0) |
| Δ_k = e_{k+1}b_k - e_kb_{k+1} where e = (1+x)I | SCC quantity |
| T1 = b_{k-1}·d_k, T2 = b_k·d_{k-1}, T3 = a_{k-1}·c_k | 3-term SCC decomposition |

## Verification

Path P_5 rooted at vertex 1: E = [1,4,4,1], J = [1,2], I = [1,5,6,1].
- d_0 = j_0·b_0 - j_{-1}·b_1 = 1·1 - 0·4 = 1 ≥ 0 ✓
- d_1 = j_1·b_1 - j_0·b_2 = 2·4 - 1·4 = 4 ≥ 0 ✓
- d_2 = j_2·b_2 - j_1·b_3 = 0·4 - 2·1 = -2 < 0
- c_2 = 4² - 4·1 = 12, |d_2| = 2, c_2/|d_2| = 6. ✓ T3 rescues.
- Check: d_1 = 4 ≥ 0, and c_2 ≥ |d_2| ⟺ 12 ≥ 2 ⟺ d_1 ≥ 0. ✓ consistent.
