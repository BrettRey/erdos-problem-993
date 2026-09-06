import Mathlib.Tactic

namespace DepthThree.RootedForestCertificate

/-- A plane rooted forest: the first tree's child forest and the remaining trees.
    Choosing an order only enlarges the 1,842 unordered types to 16,796 cases. -/
inductive Forest where
  | nil
  | cons (children rest : Forest)
  deriving DecidableEq, Repr

def size : Forest → ℕ
  | .nil => 0
  | .cons children rest => size children + size rest + 1

def enumerateAux : ℕ → ℕ → List Forest
  | 0, _ => []
  | _ + 1, 0 => [.nil]
  | fuel + 1, n + 1 => (List.range (n + 1)).flatMap fun i =>
      (enumerateAux fuel i).flatMap fun children =>
        (enumerateAux fuel (n - i)).map fun rest => .cons children rest

def enumerate (n : ℕ) : List Forest := enumerateAux (n + 1) n

theorem mem_enumerateAux_size (f : Forest) (fuel : ℕ) (h : size f < fuel) :
    f ∈ enumerateAux fuel (size f) := by
  induction f generalizing fuel with
  | nil => cases fuel <;> simp_all [size, enumerateAux]
  | cons children rest ihc ihr =>
    cases fuel with
    | zero => omega
    | succ fuel =>
      have hc : size children < fuel := by simp only [size] at h; omega
      have hr : size rest < fuel := by simp only [size] at h; omega
      simp only [size, enumerateAux]
      apply List.mem_flatMap.mpr
      refine ⟨size children, List.mem_range.mpr (by omega), ?_⟩
      apply List.mem_flatMap.mpr
      refine ⟨children, ihc fuel hc, ?_⟩
      apply List.mem_map.mpr
      refine ⟨rest, ?_, rfl⟩
      simpa using ihr fuel hr

theorem mem_enumerate_size (f : Forest) : f ∈ enumerate (size f) :=
  mem_enumerateAux_size f (size f + 1) (Nat.lt_succ_self _)

theorem enumeration_complete {f : Forest} {n : ℕ} (h : size f = n) :
    f ∈ enumerate n := by simpa [h] using mem_enumerate_size f

structure TopFive where
  c0 : ℕ
  c1 : ℕ
  c2 : ℕ
  c3 : ℕ
  c4 : ℕ
  deriving DecidableEq, Repr

def add (p q : TopFive) : TopFive :=
  ⟨p.c0 + q.c0, p.c1 + q.c1, p.c2 + q.c2, p.c3 + q.c3, p.c4 + q.c4⟩

def shift (p : TopFive) : TopFive := ⟨0, p.c0, p.c1, p.c2, p.c3⟩

def mul (p q : TopFive) : TopFive :=
  ⟨p.c0*q.c0,
   p.c0*q.c1 + p.c1*q.c0,
   p.c0*q.c2 + p.c1*q.c1 + p.c2*q.c0,
   p.c0*q.c3 + p.c1*q.c2 + p.c2*q.c1 + p.c3*q.c0,
   p.c0*q.c4 + p.c1*q.c3 + p.c2*q.c2 + p.c3*q.c1 + p.c4*q.c0⟩

/-- Top five reversed coefficients of corona P and root-deleted Q, both of
    degree `size f`. Graph/evaluator correctness remains a separate bridge. -/
def profile : Forest → TopFive × TopFive
  | .nil => (⟨1,0,0,0,0⟩, ⟨1,0,0,0,0⟩)
  | .cons children rest =>
      let (pc, qc) := profile children
      let (pr, qr) := profile rest
      let rootOut := add pc (shift pc)
      (mul (add rootOut qc) pr, mul rootOut qr)

def numerator (f : Forest) : ℕ :=
  let (p, q) := profile f
  58 * p.c0 + 21 * p.c1 + 3 * p.c2 + 30 * q.c0 + 9 * q.c1 + q.c2

def denominator (f : Forest) : ℕ :=
  let p := (profile f).1
  126 * p.c0 + 84 * p.c1 + 36 * p.c2 + 9 * p.c3 + p.c4

def withinBound (f : Forest) : Bool :=
  decide (217404 * numerator f ≤ 52513 * denominator f)

set_option maxRecDepth 100000 in
set_option maxHeartbeats 100000000 in
theorem all_codes_checked : (enumerate 10).all withinBound = true := by decide +kernel

theorem rooted_forest_bound (f : Forest) (h : size f = 10) :
    217404 * numerator f ≤ 52513 * denominator f := by
  have hmem := enumeration_complete h
  have hf := List.all_eq_true.mp all_codes_checked f hmem
  exact of_decide_eq_true hf

set_option maxRecDepth 100000 in
set_option maxHeartbeats 100000000 in
theorem enumerated_count : (enumerate 10).length = 16796 := by decide +kernel

#print axioms enumeration_complete
#print axioms rooted_forest_bound
#print axioms enumerated_count

end DepthThree.RootedForestCertificate
