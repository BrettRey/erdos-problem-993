# Lamport-style transition for unimodality (proof obligation)

Goal: recast unimodality as an inductive safety property under the rooted
composition step used in the tree DP. This avoids any global “B_T is LC”
assumption (which is false).

## Rooted composition step

Let the current parent state be (P,Q), and a child subtree state be (U,V):

- P, Q: parent polynomials (exclude/include root)
- U = P_child (child excluded)
- V = Q_child (child included, with the x factor)

Then after attaching the child,

  P' = P (U + V),    Q' = Q U,
  I' = P' + Q' = I U + P V,   where I = P + Q.

## Safety property (unimodality)

Let ΔX_k := X_{k+1} - X_k. The sequence is unimodal iff:

  if ΔX_k < 0 at some k, then ΔX_j <= 0 for all j >= k.

Equivalently, once the sequence starts decreasing, it never increases again.

## Sufficient condition for safety preservation

Assume IU is unimodal, and let d = first descent index of IU (Δ(IU)_d < 0).
If for all k >= d,

  Δ(PV)_k <= -Δ(IU)_k,      (difference dominance)

then Δ(IU + PV)_k <= 0 for all k >= d, and I' is unimodal.

A stronger sufficient condition is:

  Δ(PV)_k <= 0 for all k >= d,

since Δ(IU)_k <= 0 for k >= d.

## Empirical status

Early small-n checks suggested the Lamport tail inequalities might hold
broadly, but this is false in general.

For the concrete broom-root / leaf-child model
\[
M(p,s):\ \Delta(xP)_k \le -\Delta I_k\ \text{for all }k\ge d(I),
\]
the large table in `results/knuth_leaf_table_p120_s700.json` gives explicit
counterexamples (e.g., `p=3, s=3` fails at `k=d(I)=2`).

So: the global difference-dominance condition remains a sufficient condition,
but not a universally true invariant.

## Next possible lemmata

1) **Mode ordering**: show first descent of PV occurs no earlier than the
   first descent of IU in rooted composition.
2) **Ratio interlacing** for IU and PV to ensure the difference dominance.
3) **Local injection**: map independent sets counted by PV in the tail to
   sets in IU to compare first differences.

If any of these is proved, the inductive safety argument would yield
unimodality for all trees by structural induction on the DP.

## Restricted-case lemma (broom root + leaf child)

Let $T_s=\broom(p,s)$ with $p \ge 2$ and $s \ge p$, rooted at the hub.
Write
  $P_s(x)=(1+x)^s A(x)$, $Q_s(x)=xB(x)$, $I_s=P_s+Q_s$,
where $A(x)=I(P_{p-1};x)$ and $B(x)=I(P_{p-2};x)$.
For Lamport composition with a leaf child, $(U,V)=(1,x)$ (since $P_{\text{leaf}}=1$,
$Q_{\text{leaf}}=x$), so $IU=I_s$ and $PV=xP_s$.
Let $d(F)=\min\{k:\Delta F_k<0\}$.

**Lemma (mode ordering for broom-root leaf step).**
Fix integers $p\ge 2$ and $s\ge p$. Define
\[
P_s(x)=(1+x)^s I(P_{p-1};x),\qquad
I_s(x)=P_s(x)+x\,I(P_{p-2};x).
\]
For the Lamport child update with a leaf $(U,V)=(1,x)$ (hence $IU=I_s$ and
$PV=xP_s$), one has
\[
d(PV)\ge d(IU),
\]
where $d(F)=\min\{k:\Delta F_k<0\}$.

**Sketch:** For $k \ge t+1$ where $t=\deg B=\lfloor (p-2)/2\rfloor$, we have
$\Delta I_{s,k}=\Delta P_{s,k}$ since $b_k=b_{k-1}=0$. Also
$\Delta P_{s,k}\ge 0$ for $k \le \lfloor s/2\rfloor-1$ by binomial monotonicity,
hence $d(P_s)\ge \lfloor s/2\rfloor\ge t+1$ when $s\ge p$.
Thus $d(I_s)\le d(P_s)$ and $d(PV)=d(xP_s)=d(P_s)+1\ge d(I_s)$.

This gives **mode ordering** for this specific Lamport step, but only in the
broom-root / leaf-child regime.

## Conjectural extension: broom root + path child

Let $T_s=\broom(p,s)$ rooted at the hub with $s \ge p$, and let the child
subtree be a path of length $\ell$ (i.e., $\ell$ edges, $\ell+1$ vertices)
rooted at its endpoint. Then the Lamport step uses
  $U = P_{\text{child}} = I(P_\ell;x)$ and $V = Q_{\text{child}} = x I(P_{\ell-1};x)$.

**Conjecture (mode ordering):** $d(PV) \ge d(IU)$.

**Stronger conjecture (difference dominance):**
for all $k \ge d(IU)$, $\Delta(PV)_k \le -\Delta(IU)_k$.

**Evidence:** exhaustive parameter scan for
  $2 \le p \le 6$, $p \le s \le 10$, $0 \le \ell \le 6$
found no violations of either inequality.

The proof attempt likely mirrors the leaf case:
  - express $P_s = (1+x)^s A$ and $I_s = P_s + x B$,
  - use binomial monotonicity for $P_s$ up to $\lfloor s/2\rfloor$,
  - show the $xB$ term is supported strictly below the first descent once
    $s$ is large relative to $\deg A$ and $\deg U$.

## Lemma (eventual mode ordering for broom + path child)

Fix $p \ge 2$ and path length $\ell \ge 0$. There exists
  $s_0(p,\ell) \ge p+\ell+2$
such that for all $s \ge s_0(p,\ell)$, the Lamport step for the broom root
with path child satisfies $d(PV) \ge d(IU)$.

**Proof (formal at the level of asymptotic control).**
Write $P_s=(1+x)^s A$, $Q_s=xB$ with $A=I(P_{p-1})$, $B=I(P_{p-2})$.
Let $C=I(P_\ell)$ and $D=I(P_{\ell-1})$, so $U=C$, $V=xD$.

1) The correction term $xBC$ has degree at most
   $1+\deg B + \deg C = \lfloor p/2 \rfloor + \lfloor \ell/2 \rfloor
   \le \lfloor (p+\ell)/2 \rfloor$.
   For $s \ge p+\ell+2$, we have $\lfloor s/2 \rfloor > \deg(xBC)$, so
   in the central window the coefficients of $IU=(P_s+xB)C$ agree with those
   of $P_s C$. Hence $d(IU)=d(P_s C)$ for all such $s$.

2) For any fixed polynomial $H$ with nonnegative coefficients and no internal
   zeros, there exists $s_0(H)$ such that for all $s \ge s_0(H)$,
   the first descent of $(1+x)^s H$ occurs at
     $k = \lfloor s/2\rfloor + m_H - 1$,
   where $m_H=\lceil \mu_H + \tfrac12\rceil$ and $\mu_H=H'(1)/H(1)$.
   This follows from the ratio expansion
     $r_k = 1 - (4y+2-4\mu_H)/s + O_H(1/s^2)$
   with $k=\lfloor s/2\rfloor+y$, together with the strict monotonicity of
   $r_k$ in $y$ for large $s$ (the $O_H(1/s^2)$ term can be dominated by the
   $4/s$ slope by taking $s$ sufficiently large).
   Applying this to $H=A C$ and $H=A D$ gives
   \[
     d(P_s C) = \lfloor s/2\rfloor + m_{AC} - 1,\quad
     d(P_s D) = \lfloor s/2\rfloor + m_{AD} - 1
   \]
   for all sufficiently large $s$.

3) Since means add under convolution with nonnegative coefficients,
   $\mu_{AC}=\mu_A+\mu_C$ and $\mu_{AD}=\mu_A+\mu_D$.
   For paths, $\mu_n$ is strictly increasing and $\mu_n-\mu_{n-1}<1$:
   writing $F_n=I(P_n;1)$ and $G_n=I'(P_n;1)$, we have
     $F_n=F_{n-1}+F_{n-2}$ and
     $G_n=G_{n-1}+G_{n-2}+F_{n-2}$,
   hence
     $\mu_n=\frac{F_{n-1}\mu_{n-1}+F_{n-2}(1+\mu_{n-2})}{F_n}$.
   Since $\mu_{n-1}\ge\mu_{n-2}$, this gives
     $\mu_n-\mu_{n-1}\le F_{n-2}/F_n<1$.
   Therefore $m_{AC} \le m_{AD}+1$.

4) Therefore $d(P_s C) \le d(P_s D)+1$, and since $PV=xP_s D$ has
   $d(PV)=d(P_s D)+1$, we obtain $d(PV) \ge d(IU)$ for all large $s$.

This establishes mode ordering in the broom-root + path-child regime
for sufficiently large $s$.

## Knuth-style leaf table scan (p <= 120, s <= 700)

Computed truth table
\[
M(p,s):\ \Delta(xP)_k \le -\Delta I_k\ \text{for all }k\ge d(I)
\]
for the broom-root parent state with leaf child.

Artifact:
- `results/knuth_leaf_table_p120_s700.json`

Raw parity classes over scanned box (`s = 0..700`):
- `both`: 65 values of `p`
- `only_even`: 41 values of `p`
- `only_odd`: 13 values of `p`
- `neither`: 0 values of `p`

In the relevant broom regime (`s >= p`), the classes are:
- `both`: 7 values of `p` (`{2,4,33,42,80,89,118}`)
- `only_even`: 57 values of `p`
- `only_odd`: 55 values of `p`
- `neither`: 0 values of `p`

Low-mod periodic fits (majority bucket model over full table) are not exact:
- `mod2_p_mod2_s`: 40,965 mismatches / 83,419
- `mod4_p_mod2_s`: 39,610 mismatches / 83,419
- `mod2_p_mod4_s`: 40,958 mismatches / 83,419
- `mod4_p_mod4_s`: 39,610 mismatches / 83,419

So there is no exact small-period (`mod 2` / `mod 4`) rule in `(p,s)` on this box.
Also, all observed violations occur at the boundary index `k=d(I)` (none for
`k>d(I)` in the scanned box).

### Concise conjecture (eventual parity form)

For each fixed `p >= 2`, there exist a threshold `S_p` and a nonempty parity set
`E_p ⊆ {0,1}` such that for all `s >= S_p`,
\[
M(p,s)\iff s \bmod 2 \in E_p.
\]
Moreover, empirical table evidence gives `E_p = {0,1}` for `p in {2,4}`,
and `|E_p| = 1` for every `p >= 5` (checked for `p <= 120`).

### Short proof sketch attempt (incomplete)

Write `I = P + xB` with fixed core polynomial `B` (depends on `p`) and
`P = (1+x)^s A`.

1. Empirically, every violation occurs at the boundary index `k = d(I)` (none
   were found at larger `k` in the full table).
   Hence the global condition appears equivalent to a single inequality:
   \[
   \Delta(xP)_{d(I)} \le -\Delta I_{d(I)}.
   \]
   Substituting `I = P + xB` gives
   \[
   2\Delta P_{d(I)} + \Delta(xB)_{d(I)} \le 0.
   \]

2. For fixed `p`, `d(I)` is in a bounded window around `s/2` for large `s`,
   while `B` has fixed degree. So `\Delta(xB)_{d(I)}` is an `O_p(1)` correction.

3. The main term `\Delta P_k` near `k = s/2 + y` is governed by the binomial
   ratio expansion and changes sign across adjacent `k` values in a parity-sensitive
   way when the first descent sits on a near-tie.

4. Therefore the sign of
   `2\Delta P_{d(I)} + \Delta(xB)_{d(I)}`
   should stabilize by parity for large `s`, yielding the eventual parity rule above.

Missing step: a rigorous asymptotic reduction from the full tail condition to
the boundary inequality, plus a precise sign analysis of the boundary expression.

## Shannon boundary principle (empirical, reproducible)

For a Lamport step
  \(I' = IU + PV\), let \(d=d(IU)\) be the first descent index of \(IU\).

Define the two properties:

- `Boundary-only DD failure`: any violation of
  \(\Delta(PV)_k \le -\Delta(IU)_k\) for \(k\ge d\) occurs only at \(k=d\).
- `No interior positive tail`: \(\Delta(I')_k \le 0\) for all \(k>d\).

**Empirical result (exhaustive):**
`results/shannon_transition_scan_n17_plus_samples_18_20.json` reports, for all
rooted composition steps from all trees with \(1\le n\le 17\):

- `lamport_steps`: 20,380,920
- `steps_with_descent`: 14,237,515
- `dd_failures`: 1,414,591
- `dd_failures_boundary_only`: 1,414,591
- `dd_failures_interior`: 0
- `positive_tail_after_boundary`: 0
- `pv_positive_after_d`: 13,169
- `pv_positive_after_d_plus_1`: 0
- `pv_positive_offset_hist`: `{1: 13,169}`
- `dd_boundary_fail_and_pv_positive_after_d`: 13,169
- `dd_boundary_hold_and_pv_positive_after_d`: 0

So in this full range, every DD failure is boundary-only, and no step creates
a positive first difference after the boundary.
Also, whenever \(PV\) has a positive first difference past \(d\), it occurs
only at \(k=d+1\), and only when boundary DD already fails at \(k=d\).

**Larger-\(n\) stratified samples (geng partitions):**

- \(n=18\): 1,749,708 steps, `dd_failures_interior=0`, `positive_tail_after_boundary=0`
- \(n=19\): 2,188,800 steps, `dd_failures_interior=0`, `positive_tail_after_boundary=0`
- \(n=20\): 2,432,000 steps, `dd_failures_interior=0`, `positive_tail_after_boundary=0`

The same one-step \(PV\)-overshoot pattern appears in these samples:
`pv_positive_offset_hist = {1: ...}` and `pv_positive_after_d_plus_1 = 0`
for each of \(n=18,19,20\).

Additional stratified samples for \(n=21,22,23,24\) in
`results/shannon_transition_samples_n21_24.json` show the same pattern:
`dd_failures_interior=0`, `positive_tail_after_boundary=0`,
and `pv_positive_offset_hist={1:...}` with `pv_positive_after_d_plus_1=0`
in every sampled \(n\).

### Working lemma candidate (target)

For rooted composition states arising from trees, with \(d=d(IU)\):
\[
\Delta(I')_k \le 0 \quad \text{for all } k>d.
\]
Equivalently, any possible Lamport obstruction is concentrated at the single
boundary index \(k=d\).

If proved, this collapses the Lamport proof obligation from a tail family of
inequalities to one boundary inequality.

### Sharpened structural target

A stronger empirical statement is:
\[
\Delta(PV)_k \le 0 \quad \text{for all } k \ge d+2,
\]
and if \(\Delta(PV)_{d+1}>0\), then boundary DD already fails at \(k=d\).

This isolates all nontrivial behavior to the two-index window \(\{d,d+1\}\).

### Two-index reduction lemma (deterministic)

Let \(A,B\) be coefficient sequences and define \(\Delta X_k=X_{k+1}-X_k\).
Fix an index \(d\) and assume:

1. \(\Delta A_k \le 0\) for all \(k\ge d+1\),
2. \(\Delta B_k \le 0\) for all \(k\ge d+2\).

Then:
\[
\Delta(A+B)_k \le 0 \;\; \forall k\ge d+1
\quad\Longleftrightarrow\quad
\Delta B_{d+1} \le -\Delta A_{d+1}.
\]

Proof: for \(k\ge d+2\), both differences are nonpositive by (1)-(2), so only
\(k=d+1\) can fail. At \(k=d+1\), the condition is exactly
\(\Delta B_{d+1} \le -\Delta A_{d+1}\). \(\square\)

Applied to Lamport (\(A=IU,\;B=PV\)), once we prove
\(\Delta(PV)_k\le 0\) for \(k\ge d+2\), the entire interior tail control reduces
to one check at \(k=d+1\). The scanned data supports this split:

- `pv_positive_after_d_plus_1 = 0` (so \(\Delta(PV)_k\le 0\) for \(k\ge d+2\)),
- `dd_failures_interior = 0` (so the \(k=d+1\) inequality always held in range).

### Quantitative \(k=d+1\) margin (empirical)

Define the normalized boundary-adjacent ratio
\[
R := \frac{\Delta(PV)_{d+1}}{-\Delta(IU)_{d+1}},
\]
whenever \(-\Delta(IU)_{d+1}>0\).

From `results/shannon_transition_exhaustive_n17_profile_ratio.json`
(exhaustive \(n\le 17\)):

- `dplus1_ratio_max = 0.125`,
- `dplus1_zero_den_positive_num = 0`.

From `results/shannon_transition_samples_n18_24_profile_ratio.json`
(stratified samples for \(n=18,\dots,24\)):

- each sampled \(n\) also has `dplus1_ratio_max = 0.125`,
- again `dplus1_zero_den_positive_num = 0`.

So the observed bound is
\[
\Delta(PV)_{d+1} \le \frac18\,\bigl(-\Delta(IU)_{d+1}\bigr),
\]
with equality witness pattern at a leaf-child local configuration
(`child_size=1`, `child_arity=0`, `d=3`, `num=3`, `den=24`).

Working quantitative conjecture:
\[
\Delta(PV)_{d+1} \le \frac18\,\bigl(-\Delta(IU)_{d+1}\bigr)
\quad\text{for all rooted-tree Lamport steps.}
\]
If proved together with \(\Delta(PV)_k\le 0\) for \(k\ge d+2\), interior tail
control is immediate with a uniform slack factor.

Cross-check via dedicated slack artifact
`results/shannon_slack_scan_n17_plus_samples_18_24.json`:

- exhaustive \(n\le 17\): leaf checks `2,019,762`, violations `0`, equalities `438`,
- sampled \(n=18,\dots,24\): leaf checks all had violations `0` and
  `max_ratio = 0.125` in every sampled \(n\).

### Local profile (where boundary failures live)

From `results/shannon_transition_exhaustive_n17_profile.json` (exhaustive
\(n\le 17\)), boundary failures are strongly concentrated by child-root arity:

- arity 0: 1,131,007 / 3,801,720 (`29.75%`)
- arity 1: 239,677 / 5,852,496 (`4.10%`)
- arity 2: 37,894 / 2,785,833 (`1.36%`)
- arity 3: 5,697 / 1,121,704 (`0.51%`)
- arity 4: 310 / 424,008 (`0.073%`)
- arity 5: 6 / 157,838 (`0.0038%`)
- arity >= 6: `0` observed boundary failures

`PV` overshoots after \(d\) are even more concentrated:
- 13,133 / 13,169 occur at arity 0,
- all occur at offset \(k=d+1\), none for \(k\ge d+2\).

This suggests a proof split:
1) prove a generic high-arity monotonicity lemma (\( \text{arity} \ge 4 \) or
   at least \( \ge 6 \)),
2) handle low-arity cases (\(0,1,2,3\)) by explicit rooted recurrences.

Additional profiled samples (`results/shannon_transition_samples_n18_24_profile.json`)
support the same direction: for each sampled \(n \in \{18,\dots,24\}\),
`arity_ge6_boundary_fail = 0`; combined, this is 250,816 sampled
arity-\(\ge 6\) steps with zero observed boundary failures.

Dedicated slack scan (`results/shannon_slack_scan_n17_plus_samples_18_24.json`)
gives stronger margin data for \(r\ge 6\):

- exhaustive \(n\le 17\): `min_s0 = 4`, `min_s1 = 14`,
- sampled \(n=18,\dots,24\): `ge6_min_s0` in `{4,5}`, `ge6_min_s1 = 14`.

Here
\[
s0 = (-\Delta(IU)_d)-\Delta(PV)_d,\qquad
s1 = (-\Delta(IU)_{d+1})-\Delta(PV)_{d+1}.
\]

### Arity-threshold closure scan (new)

We ran a dedicated thresholded sign scan
(`results/shannon_arity_sign_n17_plus_samples_18_24_r4.json`,
`results/shannon_arity_sign_n17_plus_samples_18_24_r5.json`,
`results/shannon_arity_sign_n17_plus_samples_18_24.json`) that records
\(\Delta(PV)_d\), \(\Delta(PV)_{d+1}\), and \(s0,s1\) failures by arity floor.

Let \(r\) be the child-root arity in the rooted orientation.

- \(r\ge 4\): 1,185,156 scanned steps.
  `dPV_dplus1.positive = 0`, `s1_fail_count = 0`, but `s0_fail_count = 535`.
- \(r\ge 5\): 594,581 scanned steps.
  `dPV_dplus1.positive = 0`, `s1_fail_count = 0`, and `s0_fail_count = 6`.
- \(r\ge 6\): 344,732 scanned steps.
  `dPV_d.positive = 0`, `dPV_dplus1.positive = 0`,
  `s0_fail_count = 0`, `s1_fail_count = 0`.

So empirically:
1. the \(k=d+1\) obstruction disappears already at \(r\ge 4\),
2. the \(k=d\) obstruction appears confined to \(r\in\{4,5\}\),
3. no \(r\ge 6\) obstruction was seen.

### Arity-tier lemma candidates (updated)

**Candidate A4 (adjacent index safe):** if \(r\ge 4\), then
\[
\Delta(PV)_{d+1} \le -\Delta(IU)_{d+1}.
\]
Data support: zero \(s1\) failures for \(r\ge 4\) in the threshold scan above.

**Candidate H6 (high-arity boundary safe):** if \(r\ge 6\), then
\[
\Delta(PV)_d \le -\Delta(IU)_d,\qquad
\Delta(PV)_d\le 0,\qquad
\Delta(PV)_{d+1}\le 0.
\]
Data support: zero \(s0,s1\) failures and zero positive \(\Delta(PV)\) at both
critical indices for \(r\ge 6\).

**Residual finite obstruction set (if A4 + H6 hold):**
\[
r\in\{0,1,2,3,4,5\}\ \text{at }k=d,\quad
r\in\{0,1,2,3\}\ \text{at }k=d+1.
\]
This is the current closure target.

**Candidate Q1/8 (quantitative slack):**
\[
\Delta(PV)_{d+1}\le \frac18\bigl(-\Delta(IU)_{d+1}\bigr).
\]
If proved, Q1/8 implies the \(k=d+1\) inequality with uniform margin and gives
robustness against perturbative errors in any approximate argument.

### Boundary-failure catalog (exhaustive \(n\le 17\))

Using `results/shannon_boundary_fail_catalog_n17.json`, boundary failures
(\(k=d\), i.e. \(\Delta(PV)_d>-\Delta(IU)_d\)) are distributed as:

\[
r=0:1{,}131{,}007,\;
r=1:239{,}677,\;
r=2:37{,}894,\;
r=3:5{,}697,\;
r=4:310,\;
r=5:6,\;
r\ge 6:0.
\]

Focused obstruction shape:
- arity \(5\): only `6` failures total, all in one tree
  (`n=17`, `tree_index=26052`, graph6 `P???????C?G?G?E??o?B_w?[`).
- arity \(4\): `310` failures over `69` unique trees, concentrated in
  \(n\in\{14,15,16,17\}\) with `255` at \(n=17\).

So after proving \(r\ge 6\) safety, the unresolved \(k=d\) obstruction appears
to be a finite low-arity catalog, with \(r=5\) already nearly atomic.

### Exhaustive \(n=18\) extension (new)

We extended the threshold scan exhaustively to \(n=18\):

- `results/shannon_arity_sign_n18_exhaustive_r5.json`:
  \(r\ge 5\), `steps_thr = 471274`, `s0_fail = 16`, `s1_fail = 0`.
- `results/shannon_arity_sign_n18_exhaustive_r6.json`:
  \(r\ge 6\), `steps_thr = 174369`, `s0_fail = 0`, `s1_fail = 0`.

So at \(n=18\), all threshold-\(r\ge 5\) boundary failures are still confined to
arity \(5\), while arity \(\ge 6\) remains clean.

Focused boundary catalog through \(n\le 18\):
`results/shannon_boundary_fail_catalog_n18_focus5.json`.

- total arity-5 boundary failures: `22` (`n=17:6`, `n=18:16`);
- arity-5 failures occur in `4` trees total:
  - `P???????C?G?G?E??o?B_w?[` (count `6`, n=17),
  - `Q?????????O?O?G?B??M??MB_?w` (count `6`, n=18),
  - `Q?????????O?O?K?@_?E??MB?@w` (count `6`, n=18),
  - `Q?????????O?O?G?B??K??]B_?w` (count `4`, n=18).

This keeps the same closure story:
prove \(r\ge 6\) safety generically, then treat \(r=5\) as a very small explicit
exception family, and handle \(r\le 4\) separately.

### Exhaustive \(n=19\) extension (new)

We further extended threshold scans exhaustively to \(n=19\):

- `results/shannon_arity_sign_n19_exhaustive_r6.json`:
  \(r\ge 6\), `steps_thr = 499536`, `s0_fail = 0`, `s1_fail = 0`.
- `results/shannon_arity_sign_n19_exhaustive_r5.json`:
  \(r\ge 5\), `steps_thr = 1358388`, `s0_fail = 36`, `s1_fail = 0`.

So the same split persists at \(n=19\): arity \(\ge 6\) remains clean, while
remaining boundary failures at threshold \(r\ge 5\) are still at \(r=5\).

Focused arity-5 catalog through \(n\le 19\):
`results/shannon_boundary_fail_catalog_n19_focus5.json`.

- arity-5 boundary failures total: `58`, with by-\(n\) counts
  `n=17:6`, `n=18:16`, `n=19:36`;
- arity-5 failures occur in `12` witness trees (graph6 listed in artifact);
- max observed arity-5 boundary gap: `96` at
  `R???????????_?O?C??o?@_?@{F??w` (n=19).

Empirical status of Candidate H6 is now:
exhaustive through \(n=19\) and sampled for \(n=20,\dots,24\), with no
\(r\ge 6\) boundary failures observed.

### Repro command

```bash
python notes/shannon_transition_scan.py \
  --min-n 1 --max-n 17 --backend geng \
  --sample-n 18,19,20 --sample-mod 16 --sample-per-partition 400 \
  --out results/shannon_transition_scan_n17_plus_samples_18_20.json
```

Extended sampled check:
```bash
python notes/shannon_transition_scan.py \
  --min-n 1 --max-n 0 --backend geng \
  --sample-n 21,22,23,24 --sample-mod 16 --sample-per-partition 200 \
  --out results/shannon_transition_samples_n21_24.json
```

Profile command:
```bash
python notes/shannon_transition_scan.py \
  --min-n 1 --max-n 17 --backend geng \
  --out results/shannon_transition_exhaustive_n17_profile.json
```

Profiled larger-\(n\) sample command:
```bash
python notes/shannon_transition_scan.py \
  --min-n 1 --max-n 0 --backend geng \
  --sample-n 18,19,20,21,22,23,24 \
  --sample-mod 16 --sample-per-partition 200 \
  --out results/shannon_transition_samples_n18_24_profile.json
```

Ratio-profile exhaustive command:
```bash
python notes/shannon_transition_scan.py \
  --min-n 1 --max-n 17 --backend geng \
  --out results/shannon_transition_exhaustive_n17_profile_ratio.json
```

Ratio-profile larger-\(n\) sample command:
```bash
python notes/shannon_transition_scan.py \
  --min-n 1 --max-n 0 --backend geng \
  --sample-n 18,19,20,21,22,23,24 \
  --sample-mod 16 --sample-per-partition 200 \
  --out results/shannon_transition_samples_n18_24_profile_ratio.json
```

Dedicated slack scan command:
```bash
python notes/shannon_slack_scan.py \
  --min-n 1 --max-n 17 --backend geng \
  --sample-n 18,19,20,21,22,23,24 \
  --sample-mod 16 --sample-per-partition 200 \
  --out results/shannon_slack_scan_n17_plus_samples_18_24.json
```

Arity-threshold sign scans:
```bash
python notes/shannon_arity_sign_scan.py \
  --arity-threshold 4 \
  --min-n 1 --max-n 17 --backend geng \
  --sample-n 18,19,20,21,22,23,24 \
  --sample-mod 16 --sample-per-partition 200 \
  --out results/shannon_arity_sign_n17_plus_samples_18_24_r4.json

python notes/shannon_arity_sign_scan.py \
  --arity-threshold 5 \
  --min-n 1 --max-n 17 --backend geng \
  --sample-n 18,19,20,21,22,23,24 \
  --sample-mod 16 --sample-per-partition 200 \
  --out results/shannon_arity_sign_n17_plus_samples_18_24_r5.json

python notes/shannon_arity_sign_scan.py \
  --arity-threshold 6 \
  --min-n 1 --max-n 17 --backend geng \
  --sample-n 18,19,20,21,22,23,24 \
  --sample-mod 16 --sample-per-partition 200 \
  --out results/shannon_arity_sign_n17_plus_samples_18_24.json
```

Boundary-failure catalog:
```bash
python notes/shannon_boundary_fail_catalog.py \
  --max-n 17 --backend geng \
  --focus-arities 4,5 \
  --out results/shannon_boundary_fail_catalog_n17.json
```

Exhaustive \(n=18\) threshold extension:
```bash
python notes/shannon_arity_sign_scan.py \
  --arity-threshold 5 \
  --min-n 18 --max-n 18 --backend geng \
  --out results/shannon_arity_sign_n18_exhaustive_r5.json

python notes/shannon_arity_sign_scan.py \
  --arity-threshold 6 \
  --min-n 18 --max-n 18 --backend geng \
  --out results/shannon_arity_sign_n18_exhaustive_r6.json
```

Focused arity-5 catalog through \(n\le 18\):
```bash
python notes/shannon_boundary_fail_catalog.py \
  --max-n 18 --backend geng \
  --focus-arities 5 \
  --out results/shannon_boundary_fail_catalog_n18_focus5.json
```

Exhaustive \(n=19\) threshold extension:
```bash
python notes/shannon_arity_sign_scan.py \
  --arity-threshold 6 \
  --min-n 19 --max-n 19 --backend geng \
  --out results/shannon_arity_sign_n19_exhaustive_r6.json

python notes/shannon_arity_sign_scan.py \
  --arity-threshold 5 \
  --min-n 19 --max-n 19 --backend geng \
  --out results/shannon_arity_sign_n19_exhaustive_r5.json
```

Focused arity-5 catalog through \(n\le 19\):
```bash
python notes/shannon_boundary_fail_catalog.py \
  --max-n 19 --backend geng \
  --focus-arities 5 \
  --out results/shannon_boundary_fail_catalog_n19_focus5.json
```

### Order-sensitivity stress test

Because Lamport composition is sequential, intermediate states depend on child
attachment order. To stress this, we also scanned child-order permutations:

- exact permutations for parent degree `<=7`,
- `120` random distinct orders for higher degrees,
- all trees up to `n=12`, all roots, all parent states.

Artifact:
- `results/shannon_order_sensitivity_n12.json`

Totals:
- `rooted_states_tested`: 28,733
- `orders_tested`: 1,346,584
- `interior_violations`: 0

So within this range, the boundary-only phenomenon is robust to tested
child-order variation.

Repro command:
```bash
python notes/shannon_order_sensitivity_scan.py \
  --max-n 12 --backend geng --exact-degree 7 \
  --sample-orders 120 --seed 0 \
  --out results/shannon_order_sensitivity_n12.json
```
