# Squeeze Scan Analysis (n <= 23)

**Data:** `results/squeeze_n12_networkx.json`, `results/squeeze_n23_geng.json`.

## Key findings

1) The strict squeeze condition

   `first_descent >= tail_start = ceil((2α-1)/3)`

   fails starting at n = 9.

2) Worst margin through n <= 23 is **-3**, i.e.

   `first_descent = tail_start - 3`.

3) The margin -3 first appears at n=21 and persists through n=23.

   - n=21: 15 trees total (7 double-stars, 8 other).
   - n=22: 267 trees total.
   - n=23: 3,513 trees total.

So any deterministic squeeze proof must allow a margin of at least -2.
The data now force margin at least -3.

## Implication (empirical only)

For the scanned range, the relaxed condition

  `first_descent >= tail_start - 3`

is tight at n <= 23.

But this is **not** a global theorem for all trees.

Counterexample (star):
  \(S_{28}=K_{1,28}\) has \(\alpha=28\),
  \(t=\lceil(2\alpha-1)/3\rceil=\lceil 55/3\rceil=19\),
  and first descent \(d=15\), so
  \(d=t-4<t-3\).

Therefore the global target \(d(I)\ge t-3\) is false.

## Boundary check at t-2 and t-1 (empirical)

Let t = ceil((2α-1)/3).

1) Margin = -2 (n <= 20): when d(I)=t-2, the boundary step always satisfies
   i_{t-1} >= i_t. No counterexample found.

2) Margin = -3 (n <= 23): when d(I)=t-3, the boundary steps satisfy
   i_{t-2} >= i_{t-1} >= i_t in all cases checked.
   - n=21: 15 cases, all pass.
   - n=22: 267 cases, all pass.
   - n=23: 3,513 cases, all pass.

## Reframed target: MBI + leaf-heavy reduction

Since \(d(I)\ge t-3\) fails globally, use a minimal-bad-index setup:

1) Assume there exists a tree with \(d \le t-4\), and choose one with minimal
   order (MBI witness). Let \(b\) be the first index with \(\Delta i_b<0\),
   so \(b \le t-4\).
2) Apply leaf recurrence
     \(I(T)=I(T-v)+xI(T-N[v])\)
   at index \(b\) to force strong local negativity in one (or both) summands.
3) Prove this forces a leaf-heavy structure (large leaf mass near one hub,
   star/double-star-like extremal geometry).
4) Handle that reduced leaf-heavy class directly by binomial-style inequalities
   at boundary indices.

So squeeze remains a useful finite-n signal, but the proof target should be a
structural reduction rather than a universal \(d\)-vs-\(t\) bound.
