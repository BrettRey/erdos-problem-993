import RequestProject.ClanAudit.Bridge
import RequestProject.ClanAudit.GlobalPartition

/-!
# Log-concavity of the independence polynomial of the adjacent two-hub tree

Combining

* the degreewise theorem `sum_normalizedTwoRowCoeff_nonneg` of
  `RequestProject.ClanAudit.GlobalPartition` — for the adjacent two-hub tree
  `dgraph m a n b` with arbitrary finite ordered families of positive-length pendant
  paths and arbitrary total order — with
* the coefficient bridge `indepCount_logConcave_of_degreewise_nonneg` of
  `RequestProject.ClanAudit.Bridge`,

gives the *unconditional* log-concavity of the independence-count sequence of the
adjacent two-hub tree.  This closes the gap recorded in `RESULT.md`.
-/

namespace ClanAudit

open Finset

variable (m : ℕ) (a : Fin m → ℕ) (n : ℕ) (b : Fin n → ℕ)

/-- The clan-map space of the adjacent two-hub tree is the generic clan-map space. -/
theorem mapsOfOrder_eq_multMaps (N : ℕ) :
    mapsOfOrder m a n b N = multMaps (DV m a n b) N := by
  ext α
  rw [mem_mapsOfOrder, mem_multMaps]

/-- The degreewise theorem, restated over the generic clan-map space. -/
theorem dgraph_degreewise_nonneg (N k l : ℕ) (hl : 1 ≤ l) (hkl : l ≤ k) (hN : k + l = N) :
    0 ≤ ∑ α ∈ multMaps (DV m a n b) N, normalizedTwoRowCoeff (dgraph m a n b) α k l := by
  rw [← mapsOfOrder_eq_multMaps]
  exact sum_normalizedTwoRowCoeff_nonneg m a n b N k l hl hkl hN

/-- **Log-concavity of the independence-count sequence of the adjacent two-hub tree.** -/
theorem dgraph_indepCount_logConcave (j : ℕ) :
    indepCount (dgraph m a n b) j * indepCount (dgraph m a n b) (j + 2)
      ≤ indepCount (dgraph m a n b) (j + 1) * indepCount (dgraph m a n b) (j + 1) :=
  indepCount_logConcave_of_degreewise_nonneg (dgraph m a n b)
    (fun N k l hl hkl hN => dgraph_degreewise_nonneg m a n b N k l hl hkl hN) j

/-- **Log-concavity of the independence polynomial of the adjacent two-hub tree.** -/
theorem dgraph_indepPoly_logConcave (j : ℕ) :
    (indepPoly (dgraph m a n b)).coeff j * (indepPoly (dgraph m a n b)).coeff (j + 2)
      ≤ (indepPoly (dgraph m a n b)).coeff (j + 1)
        * (indepPoly (dgraph m a n b)).coeff (j + 1) :=
  indepPoly_logConcave_of_degreewise_nonneg (dgraph m a n b)
    (fun N k l hl hkl hN => dgraph_degreewise_nonneg m a n b N k l hl hkl hN) j

end ClanAudit
