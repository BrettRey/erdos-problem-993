import RequestProject.ClanAudit.Bridge
import RequestProject.ClanAudit.ConnectorPartition

/-!
# Log-concavity of the independence polynomial of the arbitrary-connector two-hub tree

Combining the connector degreewise theorem `conn_sum_normalizedTwoRowCoeff_nonneg` with
the coefficient bridge `indepCount_logConcave_of_degreewise_nonneg` gives the
*unconditional* log-concavity of the independence-count sequence — and hence of the
independence polynomial — of `connectorGraph ell m a n b` for every `ell ≥ 1` and
arbitrary finite ordered families of positive-length pendant paths at the two hubs.
-/

namespace ClanAudit

open Finset

variable (t m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ)

/-- The clan-map space of the connector tree is the generic clan-map space. -/
theorem connMapsOfOrder_eq_multMaps (N : ℕ) :
    connMapsOfOrder t m a n b N = multMaps (ConnV t m a n b) N := by
  ext α
  rw [mem_connMapsOfOrder, mem_multMaps]

/-- The connector degreewise theorem, restated over the generic clan-map space. -/
theorem connGraph_degreewise_nonneg (N k l : ℕ) (hl : 1 ≤ l) (hkl : l ≤ k) (hN : k + l = N) :
    0 ≤ ∑ α ∈ multMaps (ConnV t m a n b) N,
      normalizedTwoRowCoeff (connGraph t m a n b) α k l := by
  rw [← connMapsOfOrder_eq_multMaps]
  exact conn_sum_normalizedTwoRowCoeff_nonneg t m a n b N k l hl hkl hN

/-- **The connector degreewise theorem in terms of the number `ell ≥ 1` of connector
edges.**

The hypothesis `hell : 1 ≤ ell` was requested explicitly and is kept; it is in fact not
needed, since `connectorGraph 0 m a n b` unfolds to the adjacent two-hub tree
`connGraph 0 m a n b`, for which the statement also holds. -/
theorem connectorGraph_degreewise_nonneg (ell : ℕ) (hell : 1 ≤ ell) (N k l : ℕ) (hl : 1 ≤ l)
    (hkl : l ≤ k) (hN : k + l = N) :
    0 ≤ ∑ α ∈ multMaps (ConnV (ell - 1) m a n b) N,
      normalizedTwoRowCoeff (connectorGraph ell m a n b) α k l := by
  obtain ⟨e, rfl⟩ : ∃ e, ell = e + 1 := ⟨ell - 1, by omega⟩
  exact connGraph_degreewise_nonneg e m a n b N k l hl hkl hN

/-- **Log-concavity of the independence-count sequence of the connector tree.** -/
theorem connGraph_indepCount_logConcave (j : ℕ) :
    indepCount (connGraph t m a n b) j * indepCount (connGraph t m a n b) (j + 2)
      ≤ indepCount (connGraph t m a n b) (j + 1) * indepCount (connGraph t m a n b) (j + 1) :=
  indepCount_logConcave_of_degreewise_nonneg (connGraph t m a n b)
    (fun N k l hl hkl hN => connGraph_degreewise_nonneg t m a n b N k l hl hkl hN) j

/-- **Log-concavity of the independence polynomial of the connector tree.** -/
theorem connGraph_indepPoly_logConcave (j : ℕ) :
    (indepPoly (connGraph t m a n b)).coeff j * (indepPoly (connGraph t m a n b)).coeff (j + 2)
      ≤ (indepPoly (connGraph t m a n b)).coeff (j + 1)
        * (indepPoly (connGraph t m a n b)).coeff (j + 1) :=
  indepPoly_logConcave_of_degreewise_nonneg (connGraph t m a n b)
    (fun N k l hl hkl hN => connGraph_degreewise_nonneg t m a n b N k l hl hkl hN) j

/-- **Log-concavity of the independence polynomial of `connectorGraph ell m a n b`** for
every number `ell ≥ 1` of connector edges.

As above, the requested hypothesis `hell : 1 ≤ ell` is kept but is not needed. -/
theorem connectorGraph_indepPoly_logConcave (ell : ℕ) (hell : 1 ≤ ell) (j : ℕ) :
    (indepPoly (connectorGraph ell m a n b)).coeff j
        * (indepPoly (connectorGraph ell m a n b)).coeff (j + 2)
      ≤ (indepPoly (connectorGraph ell m a n b)).coeff (j + 1)
        * (indepPoly (connectorGraph ell m a n b)).coeff (j + 1) := by
  obtain ⟨e, rfl⟩ : ∃ e, ell = e + 1 := ⟨ell - 1, by omega⟩
  exact connGraph_indepPoly_logConcave e m a n b j

end ClanAudit
