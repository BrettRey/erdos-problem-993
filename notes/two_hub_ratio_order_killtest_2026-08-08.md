# Two-hub ratio-order kill test — 2026-08-08

## Verdict

The corrected partial ratio order does not directly prove the conjectured
log-concavity of trees with at most two branch vertices. Two natural
sufficient-condition decompositions already fail on double stars of orders
eight and six. Both witness polynomials nevertheless remain log-concave and
unimodal. These are exact no-go results for proof routes, not counterexamples
to the tree conjecture.

## Obstruction 1: grouping the adjacent-hub polynomial

For the double star with three leaves at each of its adjacent hubs,

```text
A = (1+x)^6,
D = 2x(1+x)^3,
I = A+D = 1+8x+21x^2+26x^3+17x^4+6x^5+x^6.
```

Both `A` and `D` are log-concave. However, they are neither synchronized nor
partially synchronized in the corrected senses needed by the standard sum
theorems. The exact witness is stored in the JSON certificate. Direct integer
gaps verify that `I` itself is log-concave.

Thus a proof cannot simply group the three adjacent-hub state terms into
`A+D` and invoke either synchronization closure theorem. The missing
positivity is aggregate: the self log-concavity gaps compensate for negative
cross terms.

## Obstruction 2: subdivision-contraction summands

Let `T` be the double star with two leaves at each hub and subdivide the hub
edge. Contracting that edge gives the four-leaf star. Exact tree DP gives

```text
I(T)   = 1+6x+10x^2+6x^3+x^4,
I(T/e) = 1+5x+6x^2+4x^3+x^4,
I(T_e) = I(T)+xI(T/e).
```

The two summands `I(T)` and `xI(T/e)` are not synchronized, although
`I(T)`, `I(T/e)`, and `I(T_e)` are all log-concave and unimodal. Therefore
the corrected subdivision identity plus the synchronization sum theorem does
not supply even a class-restricted subdivision proof for two-hub trees.

## Consequence

The `C_2` conjecture (every tree with at most two branch vertices is
log-concave) survives, but its proof needs an aggregate inequality that keeps
the positive self-gaps instead of requiring every cross-gap to be
nonnegative. The smallest useful algebraic target is the full Turán gap of

```text
Q_1Q_2 + x(R_1Q_2 + Q_1R_2)
```

for adjacent hubs, followed by a connector-path transfer. Pairwise or grouped
synchronization is too strong.

The factorization in `notes/double_star_log_concavity_2026-08-08.md` gives a
short self-contained proof for all double stars, although this base theorem
was already known from Galvin--Hilyard. One connector subdivision is also now
proved. Ordinary synchronization fails again inside the longer-connector
recurrence, but product-binomial-basis inequalities establish the weaker
partial-synchronization relations needed to prove every connector length;
see `notes/connector_partial_sync_route_2026-08-08.md`.

## Replay

```text
python3 verify_two_hub_ratio_order_obstructions_20260808.py
```

Certificate:

```text
results/two_hub_ratio_order_obstructions_20260808.json
```
