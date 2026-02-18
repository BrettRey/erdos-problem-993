# WHNC Submodularity Progress (2026-02-17)

This note records the two analytic pieces requested after the submodularity
impasse.

## Setup

For heavy set `H = {v : P(v) > 1/3}` define, for non-empty `S ⊆ H`:

- `demand(h) = P(h) - 1/3`
- `supply(u) = 1/3 - P(u)` for `u ∈ N(H)`
- `F(S) = supply(N(S)) - demand(S)`

Removal marginal:

`M(h,S) = F(S) - F(S\{h}) = supply(N_priv(h,S)) - demand(h)`,

where `N_priv(h,S) = N(h) \ N(S\{h})`.

We use Theorem 6 (edge bound):

`P(x) + P(y) < 2/3` for every edge `xy`.

Equivalently on any edge `h-u` with `h ∈ H`:

`supply(u) - demand(h) = 2/3 - (P(u)+P(h)) > 0`.

## Lemma A (Non-leaf private-neighbor peeling)

Assume:

1. `S0 ⊆ H_nonleaf := {h ∈ H : deg(h) ≥ 2}`, `S0 ≠ ∅`.
2. Every non-singleton `S ⊆ H_nonleaf` contains some `h ∈ S` with a private
   neighbor `u ∈ N_priv(h,S)`.

Then iterative peeling by private-neighbor vertices terminates at a singleton,
and for the final singleton `{h*}`:

`F(S0) > F({h*}) ≥ min_{h∈S0} F({h})`.

In particular, singleton dominance holds on `H_nonleaf` subsets.

### Proof

If `|S_t|=1`, stop.

If `|S_t|≥2`, assumption (2) gives `h_t ∈ S_t` with a private neighbor
`u_t ∈ N_priv(h_t,S_t)`. Then:

`M(h_t,S_t) = supply(N_priv(h_t,S_t)) - demand(h_t)`

`≥ supply(u_t) - demand(h_t)`

`= 2/3 - (P(u_t)+P(h_t)) > 0`  (Theorem 6 on edge `u_t h_t`).

So:

`F(S_t) = F(S_t\{h_t}) + M(h_t,S_t) > F(S_t\{h_t})`.

Set `S_{t+1}=S_t\{h_t}`. Size decreases by 1 each step, so in at most
`|S0|-1` steps we terminate at a singleton `{h*}`.

Summing strict inequalities along the chain:

`F(S0) > F({h*})`.

Since `{h*} ⊆ S0`, also `F({h*}) ≥ min_{h∈S0} F({h})`.

Hence `F(S0) > min_{h∈S0} F({h})`, proving singleton dominance on this class.
`□`

### Corollary (verified frontier)

The scan `results/whnc_nonleaf_private_n23.json` reports zero failures through
all `931,596` `d_leaf<=1` trees with `n<=23` for the private-neighbor premise
in Lemma A. Therefore, on this frontier, non-leaf peeling always reaches a
singleton and enforces singleton dominance on `H_nonleaf` subsets.

## Lemma B (Dominant signature `(k,l,r,degs)=(3,2,1,(2,2))`)

Let `S={h,a,b} ⊆ H` with:

- `a,b` heavy leaves,
- `h` heavy non-leaf,
- `N(S)={u,v}`,
- adjacency in `S∪N(S)`: `u~{a,h}`, `v~{b,h}`.

Then `F(S)>0` (strict WHNC slack).

### Cavity derivation

Use directed cavity ratios `R_{x→y}` at `λ=1`.
For side `u-a-h`, define:

- `x := R_{u→h}`,
- `y := R_{h→u}`.

Direct tree-cavity formulas give:

`P(u) = x / (1 + x + y)`,

`P(a) = (1 + y) / (2(1 + x + y))`.

So:

`P(a) + P(u) = (1+y+2x)/(2(1+x+y)) = 1/2 + P(u)/2`.

Similarly on side `v-b-h`:

`P(b) + P(v) = 1/2 + P(v)/2`.

Now expand:

`F(S) = (1/3-P(u)) + (1/3-P(v))`
`      - (P(h)-1/3) - (P(a)-1/3) - (P(b)-1/3)`.

Using the two side identities:

`F(S) = 5/3 - [P(h)+P(u)+P(v)+P(a)+P(b)]`

`     = 2/3 - P(h) - (P(u)+P(v))/2`

`     = (1/2)[(2/3-P(h)-P(u)) + (2/3-P(h)-P(v))]`.

Apply Theorem 6 to edges `hu` and `hv`:

`2/3 - P(h) - P(u) > 0`, `2/3 - P(h) - P(v) > 0`.

Their average is positive, hence `F(S)>0`.
`□`

### Consequence

The computationally dominant all-negative signature
`(3,2,1,(2,2))` has a direct analytic positivity proof; it cannot host a WHNC
counterexample.
