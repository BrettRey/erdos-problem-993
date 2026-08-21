import RequestProject.ClanAudit.Normalized

/-!
# The Li–Li–Yang–Zhang coefficient bridge

This file turns the *degreewise* normalized two-row nonnegativity statement of the clan
machinery into ordinary log-concavity of the independence-count sequence of a finite
simple graph.  Nothing is postulated: the coefficient identity

```text
∑ α with ∑ v α v = k + l, normalizedTwoRowCoeff G α k l = i_k * i_l - i_(k+1) * i_(l-1)
```

is *derived* from the clan/multicolouring definitions of
`RequestProject.ClanAudit.Clan` and the stable-pair bookkeeping of
`RequestProject.ClanAudit.Stable`.

The mechanism is the following.  For a fixed multiplicity map `α`, an ordered stable
two-block partition `(A,B)` of the clan graph `clan G α` has at most one clone of each
vertex in each block (clones of a vertex form a clique), so its *supports*
`S = π(A)`, `T = π(B)` are independent sets of `G` with `α = 1_S + 1_T` and
`|S| = |A|`, `|T| = |B|`.  Conversely each such pair of supports has exactly
`∏ v (α v)!` preimages — this is exactly the clan normalization `clan_fiber_card`,
applied to the two-colour multicolouring `colSets S T`.  Summing over `α` of a fixed
total order therefore counts *all* ordered pairs of independent sets of sizes `k` and
`l`, i.e. `i_k * i_l`.

Main results:

* `sum_normalizedTwoRowCoeff_eq` — the `(k,l)` coefficient identity;
* `sum_normalizedTwoRowCoeff_diag` — its diagonal case `l = k`, i.e.
  `i_k ^ 2 - i_(k-1) * i_(k+1)`;
* `sum_normalizedTwoRowCoeff_zero` — the `k = 0` boundary case;
* `indepCount_logConcave_of_degreewise_nonneg` — degreewise normalized two-row
  nonnegativity implies log-concavity of the independence-count sequence;
* `indepPoly_logConcave_of_degreewise_nonneg` — the same, phrased for the coefficients
  of the independence polynomial.
-/

namespace ClanAudit

open Finset

open scoped Classical

variable {V : Type*} [Fintype V] [DecidableEq V]

/-! ### Independent sets and the independence polynomial -/

/-- The independent sets of `G` with exactly `j` vertices. -/
noncomputable def indepSets (G : SimpleGraph V) (j : ℕ) : Finset (Finset V) :=
  Finset.univ.filter (fun S => IsIndep G S ∧ S.card = j)

omit [DecidableEq V] in
theorem mem_indepSets {G : SimpleGraph V} {j : ℕ} {S : Finset V} :
    S ∈ indepSets G j ↔ IsIndep G S ∧ S.card = j := by
  simp [indepSets]

/-- `indepCount G j` is the number `i(G, j)` of independent sets of `G` of size `j`;
these are the coefficients of the independence polynomial. -/
noncomputable def indepCount (G : SimpleGraph V) (j : ℕ) : ℕ := (indepSets G j).card

omit [DecidableEq V] in
theorem indepCount_zero (G : SimpleGraph V) : indepCount G 0 = 1 := by
  have : indepSets G 0 = {∅} := by
    ext S
    simp only [mem_indepSets, Finset.mem_singleton, Finset.card_eq_zero]
    constructor
    · rintro ⟨-, h⟩; exact h
    · rintro rfl; exact ⟨by simp [IsIndep], rfl⟩
  rw [indepCount, this, Finset.card_singleton]

omit [DecidableEq V] in
theorem indepCount_eq_zero_of_gt {G : SimpleGraph V} {j : ℕ} (hj : Fintype.card V < j) :
    indepCount G j = 0 := by
  have : indepSets G j = ∅ := by
    rw [Finset.eq_empty_iff_forall_notMem]
    intro S hS
    have := (mem_indepSets.mp hS).2
    have hle : S.card ≤ Fintype.card V := Finset.card_le_univ S
    omega
  rw [indepCount, this, Finset.card_empty]

/-- The independence polynomial `I(G; x) = ∑ j i(G,j) x^j`. -/
noncomputable def indepPoly (G : SimpleGraph V) : Polynomial ℕ :=
  ∑ j ∈ Finset.range (Fintype.card V + 1), Polynomial.monomial j (indepCount G j)

omit [DecidableEq V] in
theorem coeff_indepPoly (G : SimpleGraph V) (j : ℕ) :
    (indepPoly G).coeff j = indepCount G j := by
  rw [indepPoly, Polynomial.finset_sum_coeff]
  by_cases hj : j ∈ Finset.range (Fintype.card V + 1)
  · rw [Finset.sum_eq_single j]
    · simp
    · intro i _ hij; simp [Polynomial.coeff_monomial, hij]
    · intro h; exact absurd hj h
  · rw [Finset.sum_eq_zero, (indepCount_eq_zero_of_gt (by simpa using hj)).symm]
    intro i hi
    have : i ≠ j := by
      intro h; exact hj (h ▸ hi)
    simp [Polynomial.coeff_monomial, this]

/-! ### The multiplicity map of a pair of independent sets -/

/-- The multiplicity map `1_S + 1_T` of an ordered pair of sets. -/
def pairMult (S T : Finset V) : V → ℕ :=
  fun v => (if v ∈ S then 1 else 0) + (if v ∈ T then 1 else 0)

theorem sum_pairMult (S T : Finset V) : ∑ v, pairMult S T v = S.card + T.card := by
  simp [pairMult, Finset.sum_add_distrib, Finset.sum_ite_mem]

/-- The finite set of all clan maps of `V` of total order `N`. -/
noncomputable def multMaps (V : Type*) [Fintype V] [DecidableEq V] (N : ℕ) : Finset (V → ℕ) :=
  (Fintype.piFinset fun _ => Finset.range (N + 1)).filter (fun α => ∑ v, α v = N)

theorem mem_multMaps {N : ℕ} {α : V → ℕ} : α ∈ multMaps V N ↔ ∑ v, α v = N := by
  rw [multMaps, Finset.mem_filter, Fintype.mem_piFinset]
  refine ⟨fun h => h.2, fun h => ⟨fun x => Finset.mem_range.mpr ?_, h⟩⟩
  have hx : α x ≤ ∑ y, α y :=
    Finset.single_le_sum (f := α) (fun i _ => Nat.zero_le _) (Finset.mem_univ x)
  omega

/-! ### Slices of a block of the clan graph -/

/-- The clones of `v` belonging to a set `A` of clan vertices. -/
noncomputable def sliceOf {α : V → ℕ} (A : Finset (Σ v : V, Fin (α v))) (v : V) :
    Finset (Fin (α v)) :=
  Finset.univ.filter (fun i => (⟨v, i⟩ : Σ v : V, Fin (α v)) ∈ A)

omit [Fintype V] in
theorem card_sliceOf_le_one {G : SimpleGraph V} {α : V → ℕ}
    {A : Finset (Σ v : V, Fin (α v))} (hA : IsIndep (clan G α) A) (v : V) :
    (sliceOf A v).card ≤ 1 := by
  rw [Finset.card_le_one]
  intro i hi j hj
  simp only [sliceOf, Finset.mem_filter, Finset.mem_univ, true_and] at hi hj
  by_contra hne
  exact hA _ hi _ hj (Or.inl ⟨rfl, by simp [Sigma.mk.injEq, hne]⟩)

omit [Fintype V] in
theorem sliceOf_nonempty_iff {α : V → ℕ} {A : Finset (Σ v : V, Fin (α v))} (v : V) :
    (sliceOf A v).Nonempty ↔ v ∈ A.image Sigma.fst := by
  constructor
  · rintro ⟨i, hi⟩
    simp only [sliceOf, Finset.mem_filter, Finset.mem_univ, true_and] at hi
    exact Finset.mem_image.mpr ⟨_, hi, rfl⟩
  · intro hv
    obtain ⟨x, hx, hxv⟩ := Finset.mem_image.mp hv
    subst hxv
    exact ⟨x.2, by simpa [sliceOf] using hx⟩

omit [Fintype V] in
theorem card_sliceOf {G : SimpleGraph V} {α : V → ℕ} {A : Finset (Σ v : V, Fin (α v))}
    (hA : IsIndep (clan G α) A) (v : V) :
    (sliceOf A v).card = if v ∈ A.image Sigma.fst then 1 else 0 := by
  have h1 := card_sliceOf_le_one hA v
  by_cases hv : v ∈ A.image Sigma.fst
  · rw [if_pos hv]
    have := (sliceOf_nonempty_iff (A := A) v).mpr hv
    have h2 : 1 ≤ (sliceOf A v).card := Finset.card_pos.mpr this
    omega
  · rw [if_neg hv]
    rw [Finset.card_eq_zero, ← Finset.not_nonempty_iff_eq_empty]
    intro hc
    exact hv ((sliceOf_nonempty_iff (A := A) v).mp hc)

theorem card_sliceOf_add {G : SimpleGraph V} {α : V → ℕ}
    {A B : Finset (Σ v : V, Fin (α v))} (h : IsStablePair (clan G α) A B) (v : V) :
    (sliceOf A v).card + (sliceOf B v).card = α v := by
  have hdisj : Disjoint (sliceOf A v) (sliceOf B v) := by
    rw [Finset.disjoint_left]
    intro i hi hi'
    simp only [sliceOf, Finset.mem_filter, Finset.mem_univ, true_and] at hi hi'
    exact Finset.disjoint_left.mp h.2.2.1 hi hi'
  have hunion : sliceOf A v ∪ sliceOf B v = Finset.univ := by
    ext i
    simp only [sliceOf, Finset.mem_union, Finset.mem_filter, Finset.mem_univ, true_and,
      iff_true]
    have : (⟨v, i⟩ : Σ v : V, Fin (α v)) ∈ A ∪ B := by rw [h.2.2.2]; exact Finset.mem_univ _
    exact Finset.mem_union.mp this
  rw [← Finset.card_union_of_disjoint hdisj, hunion, Finset.card_univ, Fintype.card_fin]

omit [Fintype V] in
theorem card_image_fst {G : SimpleGraph V} {α : V → ℕ} {A : Finset (Σ v : V, Fin (α v))}
    (hA : IsIndep (clan G α) A) : (A.image Sigma.fst).card = A.card := by
  refine Finset.card_image_of_injOn ?_
  intro x hx y hy hxy
  by_contra hne
  exact hA _ (Finset.mem_coe.mp hx) _ (Finset.mem_coe.mp hy) (Or.inl ⟨hxy, hne⟩)

omit [Fintype V] in
theorem indep_image_fst {G : SimpleGraph V} {α : V → ℕ} {A : Finset (Σ v : V, Fin (α v))}
    (hA : IsIndep (clan G α) A) : IsIndep G (A.image Sigma.fst) := by
  intro u hu w hw huw
  obtain ⟨x, hx, rfl⟩ := Finset.mem_image.mp hu
  obtain ⟨y, hy, rfl⟩ := Finset.mem_image.mp hw
  exact hA _ hx _ hy (Or.inr huw)

/-! ### The supports of a stable pair of the clan graph -/

/-- The pair of supports of an ordered pair of sets of clan vertices. -/
noncomputable def suppPair {α : V → ℕ} (p : Finset (Σ v : V, Fin (α v)) ×
    Finset (Σ v : V, Fin (α v))) : Finset V × Finset V :=
  (p.1.image Sigma.fst, p.2.image Sigma.fst)

/-- The ordered stable two-block partitions of `H` with block sizes `k` and `l`. -/
noncomputable def stablePairsOfSize {W : Type*} [Fintype W] [DecidableEq W]
    (H : SimpleGraph W) (k l : ℕ) : Finset (Finset W × Finset W) :=
  (stablePairs H).filter (fun p => p.1.card = k ∧ p.2.card = l)

theorem card_stablePairsOfSize {W : Type*} [Fintype W] [DecidableEq W]
    (H : SimpleGraph W) (k l : ℕ) : (stablePairsOfSize H k l).card = stableCount H k l := by
  rfl

/-- Ordered pairs of independent sets of sizes `k` and `l`. -/
noncomputable def indepPairs (G : SimpleGraph V) (k l : ℕ) : Finset (Finset V × Finset V) :=
  (indepSets G k) ×ˢ (indepSets G l)

omit [DecidableEq V] in
theorem card_indepPairs (G : SimpleGraph V) (k l : ℕ) :
    (indepPairs G k l).card = indepCount G k * indepCount G l := by
  rw [indepPairs, Finset.card_product, indepCount, indepCount]

/-- **The support of a stable pair of the clan graph.**  Both supports are independent,
they have the sizes of the blocks, and the multiplicity map is recovered as `1_S + 1_T`. -/
theorem suppPair_mem {G : SimpleGraph V} {α : V → ℕ} {k l : ℕ}
    {p : Finset (Σ v : V, Fin (α v)) × Finset (Σ v : V, Fin (α v))}
    (hp : p ∈ stablePairsOfSize (clan G α) k l) :
    suppPair p ∈ (indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α) := by
  rw [stablePairsOfSize, Finset.mem_filter, mem_stablePairs] at hp
  obtain ⟨hst, hk, hl⟩ := hp
  refine Finset.mem_filter.mpr ⟨?_, ?_⟩
  · rw [indepPairs, Finset.mem_product]
    constructor
    · refine mem_indepSets.mpr ⟨indep_image_fst hst.1, ?_⟩
      show (p.1.image Sigma.fst).card = k
      rw [card_image_fst hst.1, hk]
    · refine mem_indepSets.mpr ⟨indep_image_fst hst.2.1, ?_⟩
      show (p.2.image Sigma.fst).card = l
      rw [card_image_fst hst.2.1, hl]
  · funext v
    have h1 := card_sliceOf hst.1 v
    have h2 := card_sliceOf hst.2.1 v
    have h3 := card_sliceOf_add hst v
    simp only [suppPair, pairMult]
    omega

/-! ### The two-colour multicolouring attached to a pair of independent sets -/

/-- The two-colour multicolouring of type `1_S + 1_T` attached to a pair of independent
sets: `v` gets the colour `0` if `v ∈ S` and the colour `1` if `v ∈ T`. -/
noncomputable def colSets (S T : Finset V) : V → Finset ℕ :=
  fun v => (if v ∈ S then ({0} : Finset ℕ) else ∅) ∪ (if v ∈ T then ({1} : Finset ℕ) else ∅)

omit [Fintype V] in
theorem mem_colSets {S T : Finset V} {v : V} {c : ℕ} :
    c ∈ colSets S T v ↔ (c = 0 ∧ v ∈ S) ∨ (c = 1 ∧ v ∈ T) := by
  unfold colSets
  by_cases hS : v ∈ S <;> by_cases hT : v ∈ T <;> simp [hS, hT]

omit [Fintype V] in
theorem card_colSets (S T : Finset V) (v : V) :
    (colSets S T v).card = pairMult S T v := by
  unfold colSets pairMult
  by_cases hS : v ∈ S <;> by_cases hT : v ∈ T <;> simp [hS, hT]

omit [Fintype V] in
theorem isMulticoloring_colSets {G : SimpleGraph V} {α : V → ℕ} {S T : Finset V}
    (hS : IsIndep G S) (hT : IsIndep G T) (hα : pairMult S T = α) :
    IsMulticoloring G α (colSets S T) where
  card v := by rw [card_colSets, hα]
  disj u v huv := by
    rw [Finset.disjoint_left]
    intro c hcu hcv
    rcases mem_colSets.mp hcu with ⟨rfl, hu⟩ | ⟨rfl, hu⟩
    · rcases mem_colSets.mp hcv with ⟨-, hv⟩ | ⟨hc, -⟩
      · exact hS _ hu _ hv huv
      · exact absurd hc (by norm_num)
    · rcases mem_colSets.mp hcv with ⟨hc, -⟩ | ⟨-, hv⟩
      · exact absurd hc.symm (by norm_num)
      · exact hT _ hu _ hv huv

/-! ### The fibre of the support map -/

section Fiber

variable {G : SimpleGraph V} {α : V → ℕ} {k l : ℕ} {S T : Finset V}

omit [Fintype V] [DecidableEq V] in
theorem mem_toMulti_iff {f : (Σ v : V, Fin (α v)) → ℕ} (v : V) (c : ℕ) :
    c ∈ toMulti α f v ↔ ∃ i : Fin (α v), f ⟨v, i⟩ = c := by
  simp [toMulti, eq_comm]

theorem mem_image_fst_filter {f : (Σ v : V, Fin (α v)) → ℕ} (v : V) (c : ℕ) :
    v ∈ (Finset.univ.filter (fun x => f x = c)).image Sigma.fst ↔ c ∈ toMulti α f v := by
  rw [mem_toMulti_iff]
  constructor
  · intro h
    obtain ⟨x, hx, rfl⟩ := Finset.mem_image.mp h
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx
    exact ⟨x.2, hx⟩
  · rintro ⟨i, hi⟩
    exact Finset.mem_image.mpr ⟨⟨v, i⟩, by simp [hi], rfl⟩

/-- The pair of blocks attached to a colouring of the clan graph by the colours `0`, `1`. -/
noncomputable def pairOfColoring (α : V → ℕ) (f : (Σ v : V, Fin (α v)) → ℕ) :
    Finset (Σ v : V, Fin (α v)) × Finset (Σ v : V, Fin (α v)) :=
  (Finset.univ.filter (fun x => f x = 0), Finset.univ.filter (fun x => f x = 1))

/-- The colouring of the clan graph attached to an ordered pair of blocks. -/
noncomputable def coloringOfPair (α : V → ℕ)
    (p : Finset (Σ v : V, Fin (α v)) × Finset (Σ v : V, Fin (α v))) :
    (Σ v : V, Fin (α v)) → ℕ := fun x => if x ∈ p.1 then 0 else 1

omit [Fintype V] in
/-- A colouring lying over the two-colour multicolouring only uses the colours `0`, `1`. -/
theorem coloring_val {f : (Σ v : V, Fin (α v)) → ℕ} (hm : toMulti α f = colSets S T)
    (x : Σ v : V, Fin (α v)) : f x = 0 ∨ f x = 1 := by
  have hx : f x ∈ toMulti α f x.1 := (mem_toMulti_iff _ _).mpr ⟨x.2, rfl⟩
  rw [hm] at hx
  rcases mem_colSets.mp hx with ⟨h, -⟩ | ⟨h, -⟩
  · exact Or.inl h
  · exact Or.inr h

theorem pairOfColoring_mem (hSc : S.card = k) (hTc : T.card = l)
    {f : (Σ v : V, Fin (α v)) → ℕ} (hf : IsProperColoring (clan G α) f)
    (hm : toMulti α f = colSets S T) :
    pairOfColoring α f ∈
      (stablePairsOfSize (clan G α) k l).filter (fun p => suppPair p = (S, T)) := by
  have hval := coloring_val (S := S) (T := T) hm
  have himg0 : (Finset.univ.filter (fun x => f x = 0)).image Sigma.fst = S := by
    ext v
    rw [mem_image_fst_filter, hm]
    simp [mem_colSets]
  have himg1 : (Finset.univ.filter (fun x => f x = 1)).image Sigma.fst = T := by
    ext v
    rw [mem_image_fst_filter, hm]
    simp [mem_colSets]
  have hst : IsStablePair (clan G α) (pairOfColoring α f).1 (pairOfColoring α f).2 := by
    refine ⟨?_, ?_, ?_, ?_⟩
    · intro x hx y hy hxy
      simp only [pairOfColoring, Finset.mem_filter, Finset.mem_univ, true_and] at hx hy
      exact hf x y hxy (by rw [hx, hy])
    · intro x hx y hy hxy
      simp only [pairOfColoring, Finset.mem_filter, Finset.mem_univ, true_and] at hx hy
      exact hf x y hxy (by rw [hx, hy])
    · rw [Finset.disjoint_left]
      intro x hx hx'
      simp only [pairOfColoring, Finset.mem_filter, Finset.mem_univ, true_and] at hx hx'
      rw [hx] at hx'
      exact absurd hx' (by norm_num)
    · ext x
      simp only [pairOfColoring, Finset.mem_union, Finset.mem_filter, Finset.mem_univ,
        true_and, iff_true]
      exact hval x
  refine Finset.mem_filter.mpr ⟨Finset.mem_filter.mpr ⟨mem_stablePairs.mpr hst, ?_, ?_⟩, ?_⟩
  · rw [← card_image_fst hst.1]
    show ((Finset.univ.filter (fun x => f x = 0)).image Sigma.fst).card = k
    rw [himg0, hSc]
  · rw [← card_image_fst hst.2.1]
    show ((Finset.univ.filter (fun x => f x = 1)).image Sigma.fst).card = l
    rw [himg1, hTc]
  · show ((Finset.univ.filter (fun x => f x = 0)).image Sigma.fst,
      (Finset.univ.filter (fun x => f x = 1)).image Sigma.fst) = (S, T)
    rw [himg0, himg1]

theorem coloringOfPair_spec
    {p : Finset (Σ v : V, Fin (α v)) × Finset (Σ v : V, Fin (α v))}
    (hp : p ∈ (stablePairsOfSize (clan G α) k l).filter (fun p => suppPair p = (S, T))) :
    IsProperColoring (clan G α) (coloringOfPair α p) ∧
      toMulti α (coloringOfPair α p) = colSets S T := by
  rw [Finset.mem_filter, stablePairsOfSize, Finset.mem_filter, mem_stablePairs] at hp
  obtain ⟨⟨hst, -, -⟩, hsupp⟩ := hp
  have hS' : p.1.image Sigma.fst = S := congrArg Prod.fst hsupp
  have hT' : p.2.image Sigma.fst = T := congrArg Prod.snd hsupp
  have hfilter0 : Finset.univ.filter (fun x => coloringOfPair α p x = 0) = p.1 := by
    ext x
    by_cases hx : x ∈ p.1 <;> simp [coloringOfPair, hx]
  have hfilter1 : Finset.univ.filter (fun x => coloringOfPair α p x = 1) = p.2 := by
    ext x
    by_cases hx : x ∈ p.1
    · have hx2 : x ∉ p.2 := fun hc => Finset.disjoint_left.mp hst.2.2.1 hx hc
      simp [coloringOfPair, hx, hx2]
    · have hx2 : x ∈ p.2 := (mem_snd_iff hst x).mpr hx
      simp [coloringOfPair, hx, hx2]
  refine ⟨?_, ?_⟩
  · intro x y hxy hcol
    by_cases hx : x ∈ p.1 <;> by_cases hy : y ∈ p.1
    · exact hst.1 x hx y hy hxy
    · simp [coloringOfPair, hx, hy] at hcol
    · simp [coloringOfPair, hx, hy] at hcol
    · exact hst.2.1 x ((mem_snd_iff hst x).mpr hx) y ((mem_snd_iff hst y).mpr hy) hxy
  · funext v
    ext c
    rw [← mem_image_fst_filter]
    match c with
    | 0 => rw [hfilter0, hS']; simp [mem_colSets]
    | 1 => rw [hfilter1, hT']; simp [mem_colSets]
    | (n + 2) =>
        have hempty : Finset.univ.filter (fun x => coloringOfPair α p x = n + 2) = ∅ := by
          ext x
          by_cases hx : x ∈ p.1 <;> simp [coloringOfPair, hx]
        rw [hempty]
        simp [mem_colSets]

theorem coloringOfPair_pairOfColoring {f : (Σ v : V, Fin (α v)) → ℕ}
    (hval : ∀ x, f x = 0 ∨ f x = 1) : coloringOfPair α (pairOfColoring α f) = f := by
  funext x
  by_cases hx : f x = 0
  · simp [coloringOfPair, pairOfColoring, hx]
  · have h1 : f x = 1 := (hval x).resolve_left hx
    simp [coloringOfPair, pairOfColoring, h1]

theorem pairOfColoring_coloringOfPair
    {p : Finset (Σ v : V, Fin (α v)) × Finset (Σ v : V, Fin (α v))}
    (hst : IsStablePair (clan G α) p.1 p.2) : pairOfColoring α (coloringOfPair α p) = p := by
  refine Prod.ext ?_ ?_
  · show Finset.univ.filter (fun x => coloringOfPair α p x = 0) = p.1
    ext x
    by_cases hx : x ∈ p.1 <;> simp [coloringOfPair, hx]
  · show Finset.univ.filter (fun x => coloringOfPair α p x = 1) = p.2
    ext x
    by_cases hx : x ∈ p.1
    · have hx2 : x ∉ p.2 := fun hc => Finset.disjoint_left.mp hst.2.2.1 hx hc
      simp [coloringOfPair, hx, hx2]
    · have hx2 : x ∈ p.2 := (mem_snd_iff hst x).mpr hx
      simp [coloringOfPair, hx, hx2]

/-- The fibre of the support map over `(S,T)` is in bijection with the fibre of the clan
normalization over the two-colour multicolouring `colSets S T`. -/
noncomputable def fiberEquiv (hSc : S.card = k) (hTc : T.card = l) :
    {f : (Σ v : V, Fin (α v)) → ℕ //
        IsProperColoring (clan G α) f ∧ toMulti α f = colSets S T}
      ≃ {p : Finset (Σ v : V, Fin (α v)) × Finset (Σ v : V, Fin (α v)) //
        p ∈ (stablePairsOfSize (clan G α) k l).filter (fun p => suppPair p = (S, T))} where
  toFun f := ⟨pairOfColoring α f.1, pairOfColoring_mem hSc hTc f.2.1 f.2.2⟩
  invFun p := ⟨coloringOfPair α p.1, coloringOfPair_spec p.2⟩
  left_inv f := Subtype.ext (coloringOfPair_pairOfColoring (coloring_val f.2.2))
  right_inv p := by
    refine Subtype.ext (pairOfColoring_coloringOfPair (G := G) ?_)
    exact mem_stablePairs.mp (Finset.mem_filter.mp (Finset.mem_filter.mp p.2).1).1

/-- **Each pair of supports has exactly `∏ v (α v)!` lifts.**  This is the clan
normalization `clan_fiber_card` applied to the two-colour multicolouring. -/
theorem card_fiber_suppPair (hS : IsIndep G S) (hT : IsIndep G T) (hSc : S.card = k)
    (hTc : T.card = l) (hα : pairMult S T = α) :
    ((stablePairsOfSize (clan G α) k l).filter (fun p => suppPair p = (S, T))).card
      = ∏ v : V, Nat.factorial (α v) := by
  rw [← clan_fiber_card (isMulticoloring_colSets hS hT hα),
    Nat.card_congr (fiberEquiv (G := G) hSc hTc), Nat.card_eq_finsetCard]

end Fiber

/-! ### The coefficient identity -/

/-- **Fibering the stable pairs of a fixed clan map over their supports.**  Every ordered
stable two-block partition of `clan G α` has a pair of independent supports inducing `α`,
and every such pair of supports has exactly `∏ v (α v)!` lifts. -/
theorem stableCount_clan_eq (G : SimpleGraph V) (α : V → ℕ) (k l : ℕ) :
    stableCount (clan G α) k l
      = ((indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α)).card
        * ∏ v : V, Nat.factorial (α v) := by
  have hfib : ∀ q ∈ (indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α),
      ((stablePairsOfSize (clan G α) k l).filter (fun p => suppPair p = q)).card
        = ∏ v : V, Nat.factorial (α v) := by
    rintro ⟨S, T⟩ hq
    rw [Finset.mem_filter, indepPairs, Finset.mem_product] at hq
    obtain ⟨⟨hq1, hq2⟩, hqa⟩ := hq
    obtain ⟨hS, hSc⟩ := mem_indepSets.mp hq1
    obtain ⟨hT, hTc⟩ := mem_indepSets.mp hq2
    exact card_fiber_suppPair hS hT hSc hTc hqa
  rw [← card_stablePairsOfSize,
    Finset.card_eq_sum_card_fiberwise (fun p hp => suppPair_mem hp),
    Finset.sum_congr rfl hfib, Finset.sum_const, smul_eq_mul]

/-- The normalized count of stable pairs of a fixed clan map is exactly the number of
ordered pairs of independent sets inducing it. -/
theorem stableCount_div_clanFactorial (G : SimpleGraph V) (α : V → ℕ) (k l : ℕ) :
    (stableCount (clan G α) k l : ℚ) / clanFactorial α
      = (((indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α)).card : ℚ) := by
  rw [stableCount_clan_eq]
  have h1 : (((((indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α)).card
        * ∏ v : V, Nat.factorial (α v)) : ℕ) : ℚ)
      = ((((indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α)).card : ℚ))
        * clanFactorial α := by
    rw [clanFactorial]
    push_cast
    ring
  rw [h1, mul_div_assoc, div_self (clanFactorial_ne_zero α), mul_one]

/-- Summing the support fibres over all clan maps of the right total order counts *all*
ordered pairs of independent sets of sizes `k` and `l`. -/
theorem sum_card_fiber_multMaps (G : SimpleGraph V) (k l N : ℕ) (hN : k + l = N) :
    ∑ α ∈ multMaps V N, ((indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α)).card
      = indepCount G k * indepCount G l := by
  rw [← card_indepPairs]
  refine (Finset.card_eq_sum_card_fiberwise ?_).symm
  intro q hq
  rw [Finset.mem_coe, indepPairs, Finset.mem_product] at hq
  refine Finset.mem_coe.mpr (mem_multMaps.mpr ?_)
  rw [sum_pairMult, (mem_indepSets.mp hq.1).2, (mem_indepSets.mp hq.2).2, hN]

/-- If the sizes do not add up to the total order, all support fibres are empty. -/
theorem sum_card_fiber_multMaps_of_ne (G : SimpleGraph V) (k l N : ℕ) (hN : k + l ≠ N) :
    ∑ α ∈ multMaps V N, ((indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α)).card
      = 0 := by
  refine Finset.sum_eq_zero fun α hα => ?_
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro q hq hc
  rw [indepPairs, Finset.mem_product] at hq
  have h2 : ∑ v, pairMult q.1 q.2 v = N := by rw [hc]; exact mem_multMaps.mp hα
  rw [sum_pairMult, (mem_indepSets.mp hq.1).2, (mem_indepSets.mp hq.2).2] at h2
  exact hN h2

/-- The normalized two-row coefficient of a fixed clan map, counted by pairs of
independent sets. -/
theorem normalizedTwoRowCoeff_eq_card_sub (G : SimpleGraph V) (α : V → ℕ) (k l : ℕ) :
    normalizedTwoRowCoeff G α k l
      = (((indepPairs G k l).filter (fun q => pairMult q.1 q.2 = α)).card : ℚ)
        - (((indepPairs G (k + 1) (l - 1)).filter (fun q => pairMult q.1 q.2 = α)).card : ℚ) := by
  rw [normalizedTwoRowCoeff, twoRowCoeff]
  push_cast
  rw [sub_div, stableCount_div_clanFactorial, stableCount_div_clanFactorial]

/-- **The Li–Li–Yang–Zhang coefficient bridge.**  Summing the normalized two-row
coefficients over all clan maps of total order `k + l` gives the two-row minor
`i_k i_l - i_(k+1) i_(l-1)` of the independence-count sequence. -/
theorem sum_normalizedTwoRowCoeff_eq (G : SimpleGraph V) (k l N : ℕ) (hl : 1 ≤ l)
    (hN : k + l = N) :
    ∑ α ∈ multMaps V N, normalizedTwoRowCoeff G α k l
      = (indepCount G k : ℚ) * (indepCount G l : ℚ)
        - (indepCount G (k + 1) : ℚ) * (indepCount G (l - 1) : ℚ) := by
  rw [Finset.sum_congr rfl (fun α _ => normalizedTwoRowCoeff_eq_card_sub G α k l),
    Finset.sum_sub_distrib, ← Nat.cast_sum, ← Nat.cast_sum,
    sum_card_fiber_multMaps G k l N hN,
    sum_card_fiber_multMaps G (k + 1) (l - 1) N (by omega)]
  push_cast
  ring

/-- The diagonal case of the bridge:
`∑ α normalizedTwoRowCoeff G α k k = i_k ^ 2 - i_(k-1) i_(k+1)`. -/
theorem sum_normalizedTwoRowCoeff_diag (G : SimpleGraph V) (k : ℕ) (hk : 1 ≤ k) :
    ∑ α ∈ multMaps V (2 * k), normalizedTwoRowCoeff G α k k
      = (indepCount G k : ℚ) ^ 2 - (indepCount G (k - 1) : ℚ) * (indepCount G (k + 1) : ℚ) := by
  rw [sum_normalizedTwoRowCoeff_eq G k k (2 * k) hk (by omega)]
  ring

/-- The `k = 0` boundary case, handled explicitly: there is a single clan map of order
`0`, and the corresponding sum is `i_0 ^ 2 = 1` (the term `i_1 i_(-1)` is absent). -/
theorem sum_normalizedTwoRowCoeff_zero (G : SimpleGraph V) :
    ∑ α ∈ multMaps V 0, normalizedTwoRowCoeff G α 0 0 = (indepCount G 0 : ℚ) ^ 2 := by
  rw [Finset.sum_congr rfl (fun α _ => normalizedTwoRowCoeff_eq_card_sub G α 0 0),
    Finset.sum_sub_distrib, ← Nat.cast_sum, ← Nat.cast_sum,
    sum_card_fiber_multMaps G 0 0 0 rfl,
    sum_card_fiber_multMaps_of_ne G 1 0 0 (by norm_num)]
  push_cast
  ring

/-! ### From degreewise nonnegativity to log-concavity -/

/-- **The reusable transfer theorem.**  Degreewise nonnegativity of the normalized
two-row coefficients implies log-concavity of the independence-count sequence. -/
theorem indepCount_logConcave_of_degreewise_nonneg (G : SimpleGraph V)
    (h : ∀ N k l : ℕ, 1 ≤ l → l ≤ k → k + l = N →
      0 ≤ ∑ α ∈ multMaps V N, normalizedTwoRowCoeff G α k l) (j : ℕ) :
    indepCount G j * indepCount G (j + 2) ≤ indepCount G (j + 1) * indepCount G (j + 1) := by
  have hnn := h (2 * (j + 1)) (j + 1) (j + 1) (by omega) le_rfl (by omega)
  rw [sum_normalizedTwoRowCoeff_eq G (j + 1) (j + 1) (2 * (j + 1)) (by omega) (by omega)] at hnn
  simp only [Nat.add_sub_cancel] at hnn
  have hQ : (indepCount G j * indepCount G (j + 2) : ℚ)
      ≤ (indepCount G (j + 1) * indepCount G (j + 1) : ℚ) := by nlinarith [hnn]
  exact_mod_cast hQ

/-- The same conclusion phrased for the coefficients of the independence polynomial. -/
theorem indepPoly_logConcave_of_degreewise_nonneg (G : SimpleGraph V)
    (h : ∀ N k l : ℕ, 1 ≤ l → l ≤ k → k + l = N →
      0 ≤ ∑ α ∈ multMaps V N, normalizedTwoRowCoeff G α k l) (j : ℕ) :
    (indepPoly G).coeff j * (indepPoly G).coeff (j + 2)
      ≤ (indepPoly G).coeff (j + 1) * (indepPoly G).coeff (j + 1) := by
  simpa [coeff_indepPoly] using indepCount_logConcave_of_degreewise_nonneg G h j

end ClanAudit
