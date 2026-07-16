# Frozen target: coefficient monotonicity in the window from zero geometry

## Setting and known facts

Tree T on n vertices, independence number alpha, independence
polynomial I(T;x) = sum_k i_k x^k, all i_k > 0. Write
r_k = i_{k+1}/i_k and L = ceil((2 alpha - 1)/3).

Known (cite precisely; sources in notes/literature/):

(K1) Tail: (i_k) is weakly decreasing for k >= L (Levit--Mandrescu).
(K2) Exhaustive: every tree with n <= 29 is unimodal (this repo,
     8.69e9 trees, 0 failures).
(K3) The smallest-modulus zero beta(T) is unique and real
     (Csikvari, CPC 2013).
(K4) Degree-d+1 zero-free geometry: zeros avoid the Shearer disk
     |z| <= lambda*(d) = d^d/(d+1)^{d+1}; avoid a neighborhood of
     [0, d^d/(d-1)^{d+1}) (Peters--Regts); avoid the semi-disk of
     radius (7/8) tan(pi/(2d)) in the right half-plane (Bencs--
     Csikvari 2018); improved angular profiles in Bencs--Csikvari
     arXiv:2204.04868 Thm 8.1. Rescaled limit: d*U_d -> U_infty inside
     the cardioid C_infty = {-u e^{-u} : |u| < 1}, with explicit
     boundary curve Gamma near the positive endpoint e
     (Bencs--Buys--Peters arXiv:2111.06451, Cor 5.4).
(K5) Certified empirical laws at sampled scale (data.md): non-real
     zeros obey the cusp envelope angle ~ (2(r-1))^{3/2}/6 in the
     dominance ratio r = |z|/beta (matches the C_infty cusp), and a
     positive-axis sector |arg z| >= 0.90 for n <= 160 (known to
     degrade asymptotically: BBP Remark 5.5).
(K6) Cautionary control: the split graph K_20 v E_6 has a valley at
     k = 2 with all non-real zeros at dominance ratio >= 30. Any
     bridge must therefore be asymptotic in alpha; pointwise
     "zeros far => unimodal" is FALSE.

## Target (main)

Formulate and prove an effective statement of the shape:

  BRIDGE. There exist explicit functions alpha_0(d) and a zero-set
  hypothesis H(d) — satisfied by every tree of maximum degree <= d+1
  by (K3)/(K4), or by the quantified cusp/sector form of (K5) — such
  that every tree T with max degree <= d+1 and alpha(T) >= alpha_0(d)
  has NO valley in the window: there are no indices
  a < b < c <= L with i_a > i_b < i_c.

Together with (K1) and (K2) this proves unimodality for all trees of
bounded degree with alpha large; degree-uniformity may be attempted
but is NOT required for success.

## Suggested route (not mandatory)

Saddle-point/Darboux: i_k = (2 pi i)^{-1} oint I(x) x^{-k-1} dx with
positive saddle x_k solving x I'(x)/I(x) = k. Writing
log I(x) = sum_j log(1 - x/z_j), real zeros contribute log-concave
structure along the positive axis; non-real pairs contribute
oscillation controlled by their distance to the saddle contour and
their angle. The window k <= 2 alpha/3 keeps x_k in a bounded range
(quantify!); (K4)-type sector profiles then damp oscillations. The
known danger regime is high-degree hubs where x_k may pass near the
critical activity ~ e/d: the repo's hub-bouquet families live exactly
there and are empirically unimodal with rebound deficit >= ~11/n
(data.md), so the truth is there to be proven.

## Graded sub-targets

T1 (pilot). Prove BRIDGE for spiders S(a_1,...,a_m) or for trees with
    d_leaf <= 1, with explicit alpha_0. Even a re-derivation of Li--
    Yang--Zhang--Zhang's spider unimodality by zero-geometry methods
    counts: it would be the first working instance of the bridge.

T2 (main). BRIDGE for max degree <= d+1 with explicit alpha_0(d).

T3 (refutation mode). Construct a sequence of polynomials with all
    positive coefficients, satisfying the zero-geometry hypotheses
    H(d) (including cusp and sector), whose coefficient sequences have
    window valleys with alpha -> infinity. This proves the stated
    hypotheses insufficient and identifies the missing tree-only
    input. Full success.

T4 (stretch). Degree-uniform version via beta-normalization: replace
    (K4) by the per-tree cusp envelope (K5) and prove BRIDGE with
    alpha_0 absolute.

## What is NOT in scope

- The prefix k < epsilon * alpha (repo has separate machinery:
  prefix-GSB theorems for r <= 3).
- Lean formalization (later; see issue #7 pattern).
- Any claim about general graphs.
