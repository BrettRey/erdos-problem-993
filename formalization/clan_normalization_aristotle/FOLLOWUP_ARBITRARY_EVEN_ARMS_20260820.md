# Aristotle continuation: close the arbitrary-even-arm `p = 2` weight theorem

Continue the existing `ClanAudit` Lean project. Preserve all proved declarations and
use the existing definitions wherever possible. Close the first exact unfilled
obligation in `RESULT.md`; do not begin the global adjacent-two-hub partition until
this layer builds cleanly.

## Target

Generalize `localMapP2_normalized_weight_two_arms` from exactly two odd active arms to
a hub with:

- two distinguished odd positive prefixes of lengths `L ≤ M`, as in the existing
  theorem; and
- an arbitrary finite indexed family of `e` additional active arms, each having a
  positive even prefix length.

Formalize the image clan graph after `localMapP2`: the cloned `K₂` components cancel
their factorials; each additional even-prefix arm becomes a separate balanced path
component with normalized imbalance weight `2`; the untouched odd-arm remainder
contributes `z + z⁻¹`. Derive, rather than assume,

```text
W(image state) = 2^e * (z + z⁻¹)
W(active state) + W(image state) = Ablock 1 (2^e).
```

In particular, expose a named theorem giving the derived scalar `c = 2^e` and hence
`1 ≤ c`. It is acceptable—and preferable—to prove a reusable finite-family
disjoint-union/product theorem for `imbalanceGF` or `Wpoly` if that is the cleanest
route.

## Soundness and traceability

- No `sorry`, `admit`, `axiom`, `implemented_by`, `native_decide`, or weakened target.
- Do not postulate component decompositions or normalized weights as hypotheses when
  they can be constructed and proved from the explicit graph models already present.
- Add source-labelled restatements to `RequestProject/ClanAudit.lean`.
- Update `README.md` and `RESULT.md` so this obligation is marked proved and the first
  remaining obligation is the exhaustive global clan-map block partition.
- Run `lake build`, search for escape hatches, and report `#print axioms` for the new
  principal theorems.

If the target is false, return the smallest explicit counterexample and identify the
first failed graph decomposition or weight identity. An honest refutation is success.
