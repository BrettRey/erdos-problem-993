/-
# The spider polynomials `C`, `B` and the quantity `V`

```
C = Q + x*R          (independence polynomial of the spider: hub plus arms)
B = C + x*Q
V k = c k * b k − c (k+1) * b (k-1)
```

Main results:

* `Rof_LR_Cof` : `R ≤lr C`.
* `Cof_LR_X_X_Qof` / `cross_term_nonneg` : **layer 5**, the shifted cross term
  `c k * q (k-1) − c (k+1) * q (k-2) ≥ 0`.
* `V_nonneg_of_logConcave` : the LAW-V inequality `0 ≤ V k`, given log-concavity of `C`.

Source label: agent (layer 5 and the final assembly).
-/
import RequestProject.Products

namespace LawV

open Polynomial

noncomputable section

/-- The spider polynomial `C = Q + x*R`. -/
def Cof (arms : List ℕ) : ℤ[X] := Qof arms + X * Rof arms

/-- `B = C + x*Q`. -/
def Bof (arms : List ℕ) : ℤ[X] := Cof arms + X * Qof arms

/-- `V k = c k * b k − c (k+1) * b (k-1)`, with coefficients read off at integer indices
(so that all boundary cases are covered: coefficients vanish outside the support). -/
def Vco (arms : List ℕ) (k : ℤ) : ℤ :=
  cf (Cof arms) k * cf (Bof arms) k - cf (Cof arms) (k+1) * cf (Bof arms) (k-1)

lemma Cof_niceD (arms : List ℕ) :
    NiceD (max (Qdeg arms) (Rdeg arms + 1)) (Cof arms) :=
  (Qof_niceD arms).add_X_mul (Rof_niceD arms)

lemma Cof_nice (arms : List ℕ) : Nice (Cof arms) := (Cof_niceD arms).nice

lemma Cof_cf_nonneg (arms : List ℕ) (i : ℤ) : 0 ≤ cf (Cof arms) i := (Cof_nice arms).cf_nonneg i

lemma Bof_niceD (arms : List ℕ) :
    NiceD (max (max (Qdeg arms) (Rdeg arms + 1)) (Qdeg arms + 1)) (Bof arms) :=
  (Cof_niceD arms).add_X_mul (Qof_niceD arms)

lemma Bof_nice (arms : List ℕ) : Nice (Bof arms) := (Bof_niceD arms).nice

lemma cf_X_mul_nonneg {p : ℤ[X]} (hp : ∀ i, 0 ≤ cf p i) (i : ℤ) : 0 ≤ cf (X * p) i := by
  rw [cf_X_mul]; exact hp _

/-! ### Elementary comparisons for the spider -/

/-- `R ≤lr Q` upgrades to `R ≤lr C`. -/
lemma Rof_LR_Cof (arms : List ℕ) : LR (Rof arms) (Cof arms) := by
  refine LR.add_right (Rof_LR_Qof arms) ?_
  exact logConcave_iff_LR.1 (Rof_logConcave arms)

/-- `R ≤lr X * Q`. -/
lemma Rof_LR_X_Qof (arms : List ℕ) : LR (Rof arms) (X * Qof arms) :=
  LR.trans_X (Rof_nice arms) (Rof_niceD arms) (Qof_cf_nonneg arms)
    (logConcave_iff_LR.1 (Rof_logConcave arms)) (Rof_LR_Qof arms)

/-- `Q ≤lr X^2 * Q`. -/
lemma Qof_LR_X_X_Qof (arms : List ℕ) : LR (Qof arms) (X * (X * Qof arms)) :=
  LR.trans_X (Qof_nice arms) (Qof_niceD arms) (cf_X_mul_nonneg (Qof_cf_nonneg arms))
    (logConcave_iff_LR.1 (Qof_logConcave arms)) (logConcave_iff_LR.1 (Qof_logConcave arms))

/-- **Layer 5.**  The shifted comparison `C ≤lr X^2 * Q`. -/
theorem Cof_LR_X_X_Qof (arms : List ℕ) : LR (Cof arms) (X * (X * Qof arms)) := by
  rw [Cof]
  refine LR.add_left (Qof_LR_X_X_Qof arms) ?_
  exact LR_X_mul (Rof_LR_X_Qof arms)

/-- **Layer 5, coefficient form.**  `c k * q (k-1) − c (k+1) * q (k-2) ≥ 0`. -/
theorem cross_term_nonneg (arms : List ℕ) (k : ℤ) :
    0 ≤ cf (Cof arms) k * cf (Qof arms) (k-1) - cf (Cof arms) (k+1) * cf (Qof arms) (k-2) := by
  have h := Cof_LR_X_X_Qof arms k
  simp only [cf_X_mul] at h
  have e1 : k - 1 - 1 = k - 2 := by ring
  have e2 : k + 1 - 1 = k := by ring
  rw [e1, e2] at h
  linarith [h]

/-! ### The LAW-V inequality -/

/-- `0 ≤ V k` for all `k` is exactly the statement `C ≤lr X * B`. -/
lemma Vco_nonneg_iff (arms : List ℕ) : (∀ k : ℤ, 0 ≤ Vco arms k) ↔ LR (Cof arms) (X * Bof arms) := by
  constructor
  · intro h k
    have := h k
    simp only [Vco] at this
    simp only [cf_X_mul, add_sub_cancel_right]
    linarith [this]
  · intro h k
    have := h k
    simp only [cf_X_mul, add_sub_cancel_right] at this
    simp only [Vco]
    linarith [this]

/-- **The LAW-V one-hub inequality, given the spider log-concavity input.**
If `C` is log-concave then `0 ≤ V k` for every integer index `k`. -/
theorem V_nonneg_of_logConcave (arms : List ℕ) (hLC : LogConcave (Cof arms)) (k : ℤ) :
    0 ≤ Vco arms k := by
  refine (Vco_nonneg_iff arms).2 ?_ k
  have hB : X * Bof arms = X * Cof arms + X * (X * Qof arms) := by
    rw [Bof]; ring
  rw [hB]
  exact LR.add_right (logConcave_iff_LR.1 hLC) (Cof_LR_X_X_Qof arms)

end

end LawV
