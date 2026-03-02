# Round 11, Prompt 1: Diagonal/Cross Decomposition Verification + Aggregated Condition Test

**Target model:** Codex 5.3 (exhaustive computation)

---

## Context

We are proving E ≽ J (ratio dominance) at every support vertex of every tree. At each incremental step when processing non-leaf child t, the P2 minor is:

```
w_k = d_k(E_new, J_new) = E_new_{k+1} · J_new_k - E_new_k · J_new_{k+1}
```

where E_new = E_acc · f, J_new = J_acc · g, with f = I_t = E_t + xJ_t and g = E_t.

## The Diagonal/Cross Decomposition (NEW, to be verified at larger n)

Expanding w_k = Σ_{i,j} E_acc_i · J_acc_j · φ(i,j) where:
```
φ(a,b) = f_{k+1-a} · g_{k-b} - f_{k-a} · g_{k+1-b}
```

Split into diagonal (i=j) and cross (i≠j):

**Diagonal:** D_k = Σ_i E_acc_i · J_acc_i · φ(i,i) = Σ_i E_acc_i · J_acc_i · d_{k-i}(f, g)

Note: d_p(f,g) = f_{p+1}·g_p - f_p·g_{p+1} is the factor-level LR minor. This CAN be negative (factor E≽J fails ~14% of support vertices).

**Cross:** X_k = w_k - D_k = Σ_{i≠j} E_acc_i · J_acc_j · φ(i,j)

### Prior result (n ≤ 18, 410K support vertices with s ≥ 2):
- **X_k ≥ 0 in ALL cases (0 failures)**
- **D_k ≥ 0 in ALL cases (0 failures)**
- Both parts are independently non-negative, so w_k = D_k + X_k ≥ 0 follows trivially.

This decomposition was also verified at the T_{3,4} broom (n=28, s=1):
```
k= 0: w_k=2, D_k=1, X_k=1
k= 7: w_k=2657771116, D_k=58794144, X_k=2598976972
k=14: w_k=0, D_k=0, X_k=0
All non-negative.
```

## Your Tasks

### Task 1: Extend diagonal/cross verification to n = 20

For ALL trees n = 3..20, ALL support vertices with s ≥ 2, ALL incremental steps, ALL k from 0 to mode-1:

1. Compute w_k = E_new_{k+1}·J_new_k - E_new_k·J_new_{k+1}
2. Compute D_k = Σ_i E_acc_i · J_acc_i · d_{k-i}(f, g)
3. Compute X_k = w_k - D_k
4. Report: any D_k < 0? any X_k < 0?

The key finding to confirm or refute: **D_k ≥ 0 and X_k ≥ 0 both hold universally through n=20.**

Use geng at `/opt/homebrew/bin/geng` for tree enumeration. Core DP infrastructure:
```python
from indpoly import _polymul, _polyadd  # available in cwd
```

### Task 2: Test the Aggregated Condition (Agg)

GPT 5.2 Pro identified this as the key condition for a transport proof. At each incremental step, for each k and each valid i:

```
(Agg): r_{k-i} · J_{k-1} ≤ r_{k-1-i} · J_k   for all relevant (i, k)
```

where r = J_t (the current child's J), J = J_new = J_acc · g (the accumulated J after this step).

Equivalently: r_{k-i}/r_{k-1-i} ≤ J_k/J_{k-1}.

This says the accumulated product J has "heavier tails" (in the ratio sense) than the individual factor r.

For each incremental step, for each k, for each i where both r_{k-i} > 0 and r_{k-1-i} > 0:
- Check whether (Agg) holds.
- Report: total checks, total failures, failure rate.
- For failures: report the tree (g6), n, root, step t, k, i, and the values.

### Task 3: Test the Shifted TP2 Condition (STP2)

At each incremental step, for q = E_t and r = J_t:

```
(STP2): r_{m+1} · q_n ≤ r_m · q_{n+1}   for all m > n ≥ 0
```

This is a stronger condition than (Agg) — it controls the mixed bracket pointwise without needing the accumulated J.

Check (STP2) for all tree-realizable (q, r) pairs at support vertices, n = 3..20.
Report: total checks, total failures, failure rate, smallest failing tree.

### Task 4: Profile D_k and X_k relative magnitudes

At each s ≥ 2 support vertex, record:
- min(D_k) over k
- min(X_k) over k
- ratio D_k/X_k at the k where w_k is minimized (the "tightest" index)

Report by s: the distribution of D/X ratios and whether D or X is typically dominant.

## Output Format

Print results clearly with section headers. Include total counts and failure counts. If any failures are found, print full details (g6 string, n, root, step, k, i, values).

For n=19 and n=20, partial results are fine if time is limited — just report how many trees were scanned.
