/-
# Layer 3: products of path polynomials

For a list `arms : List ℕ` of arm lengths we set

```
Qof arms = ∏ᵢ P (aᵢ)     = ∏ᵢ W (aᵢ + 1)
Rof arms = ∏ᵢ P (aᵢ - 1) = ∏ᵢ W aᵢ
```

Main results:

* `Qof_niceD`, `Rof_niceD` : the products have positive coefficients on an initial segment
  and vanish afterwards (log-concavity "with no internal zeros").
* `Qof_logConcave`, `Rof_logConcave` : both products are log-concave.
* `Rof_LR_Qof` : `R ≤lr Q`, obtained by replacing the factors one at a time.

Source label: agent (layer 3).
-/
import RequestProject.PathPoly

namespace LawV

open Polynomial

noncomputable section

/-- `Q = ∏ᵢ P (aᵢ)`, the independence polynomial of the disjoint arms. -/
def Qof : List ℕ → ℤ[X]
  | [] => 1
  | a :: t => W (a+1) * Qof t

/-- `R = ∏ᵢ P (aᵢ - 1)`, the arms with their hub-adjacent vertices deleted. -/
def Rof : List ℕ → ℤ[X]
  | [] => 1
  | a :: t => W a * Rof t

@[simp] lemma Qof_nil : Qof [] = 1 := rfl
@[simp] lemma Rof_nil : Rof [] = 1 := rfl
lemma Qof_cons (a : ℕ) (t : List ℕ) : Qof (a :: t) = W (a+1) * Qof t := rfl
lemma Rof_cons (a : ℕ) (t : List ℕ) : Rof (a :: t) = W a * Rof t := rfl

/-- The degree of `Qof arms`: the independence number of the disjoint arms. -/
def Qdeg : List ℕ → ℕ
  | [] => 0
  | a :: t => (a+1)/2 + Qdeg t

/-- The degree of `Rof arms`. -/
def Rdeg : List ℕ → ℕ
  | [] => 0
  | a :: t => a/2 + Rdeg t

lemma Qof_niceD : ∀ arms : List ℕ, NiceD (Qdeg arms) (Qof arms)
  | [] => NiceD.one
  | (a :: t) => (W_niceD (a+1)).mul (Qof_niceD t)

lemma Rof_niceD : ∀ arms : List ℕ, NiceD (Rdeg arms) (Rof arms)
  | [] => NiceD.one
  | (a :: t) => (W_niceD a).mul (Rof_niceD t)

lemma Qof_nice (arms : List ℕ) : Nice (Qof arms) := (Qof_niceD arms).nice
lemma Rof_nice (arms : List ℕ) : Nice (Rof arms) := (Rof_niceD arms).nice

lemma Qof_cf_nonneg (arms : List ℕ) (i : ℤ) : 0 ≤ cf (Qof arms) i := (Qof_nice arms).cf_nonneg i
lemma Rof_cf_nonneg (arms : List ℕ) (i : ℤ) : 0 ≤ cf (Rof arms) i := (Rof_nice arms).cf_nonneg i

/-- **Layer 3.**  Products of path polynomials are log-concave. -/
lemma Qof_logConcave : ∀ arms : List ℕ, LogConcave (Qof arms)
  | [] => by rw [Qof_nil]; exact NiceD.one.logConcave_zero
  | (a :: t) => by
      rw [Qof_cons]
      exact LogConcave.mul (W_niceD (a+1)) (W_logConcave (a+1)) (Qof_niceD t) (Qof_logConcave t)

lemma Rof_logConcave : ∀ arms : List ℕ, LogConcave (Rof arms)
  | [] => by rw [Rof_nil]; exact NiceD.one.logConcave_zero
  | (a :: t) => by
      rw [Rof_cons]
      exact LogConcave.mul (W_niceD a) (W_logConcave a) (Rof_niceD t) (Rof_logConcave t)

/-- **Layer 3.**  Replacing the factors one at a time gives `R ≤lr Q`. -/
lemma Rof_LR_Qof : ∀ arms : List ℕ, LR (Rof arms) (Qof arms)
  | [] => by rw [Rof_nil, Qof_nil]; exact LR.refl _
  | (a :: t) => by
      rw [Rof_cons, Qof_cons]
      -- first replace the tail, then the head factor
      have h1 : LR (W a * Rof t) (W a * Qof t) :=
        LR.mul_left (N := max (Rdeg t) (Qdeg t)) (W_niceD a) (W_logConcave a)
          (fun j hj => (Rof_niceD t).2 j (by omega)) (fun j hj => (Qof_niceD t).2 j (by omega))
          (LR.gen (Rof_LR_Qof t) (Rof_niceD t) (Qof_niceD t))
      have h2 : LR (W a * Qof t) (W (a+1) * Qof t) := by
        have := LR.mul_left (N := (a+1)/2 + 1) (Qof_niceD t) (Qof_logConcave t)
          (fun j hj => (W_niceD a).2 j (by omega)) (fun j hj => (W_niceD (a+1)).2 j (by omega))
          (W_LRgen_succ a)
        simpa only [mul_comm] using this
      exact LR.trans' (Nice.mul (W_nice a) (Rof_nice t))
        ((W_niceD a).mul (Qof_niceD t)) (fun i => (Nice.mul (W_nice (a+1)) (Qof_nice t)).cf_nonneg i)
        h1 h2

/-- All `2×2` minors for the pair `(R, Q)`. -/
lemma Rof_LRgen_Qof (arms : List ℕ) : LRgen (Rof arms) (Qof arms) :=
  LR.gen (Rof_LR_Qof arms) (Rof_niceD arms) (Qof_niceD arms)

end

end LawV
