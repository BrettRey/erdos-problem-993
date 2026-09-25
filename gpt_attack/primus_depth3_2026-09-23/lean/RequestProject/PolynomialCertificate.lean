import RequestProject.RootedForestCertificate
import Mathlib.Algebra.Polynomial.Coeff

namespace DepthThree.PolynomialCertificate

open Polynomial RootedForestCertificate

/-- The attachment-concentration identity in a form using only addition.
Here the full factors are q1+r1 and q2+r2, so r1,r2 count excluded sets. -/
theorem concentration_identity {R : Type*} [CommSemiring R]
    (a x q1 q2 r1 r2 : R) :
    (a * (q1+r1) * (q2+r2) + x*q1*q2) * (a+x) =
      (a*(q1+r1)+x*q1) * (a*(q2+r2)+x*q2) + a*x*r1*r2 := by ring

theorem concentration_coeff_le (a x q1 q2 r1 r2 : ℕ[X]) (k : ℕ) :
    ((a*(q1+r1)+x*q1) * (a*(q2+r2)+x*q2)).coeff k ≤
      ((a*(q1+r1)*(q2+r2)+x*q1*q2)*(a+x)).coeff k := by
  rw [concentration_identity, coeff_add]
  exact Nat.le_add_right _ _

def firstFive (p : ℕ[X]) : TopFive :=
  ⟨p.coeff 0, p.coeff 1, p.coeff 2, p.coeff 3, p.coeff 4⟩

theorem firstFive_add (p q : ℕ[X]) :
    firstFive (p+q) = RootedForestCertificate.add (firstFive p) (firstFive q) := by
  simp [firstFive, RootedForestCertificate.add]

theorem firstFive_X_mul (p : ℕ[X]) :
    firstFive (X*p) = shift (firstFive p) := by
  simp [firstFive, shift, coeff_X_mul]

theorem firstFive_mul (p q : ℕ[X]) :
    firstFive (p*q) = RootedForestCertificate.mul (firstFive p) (firstFive q) := by
  simp [firstFive, RootedForestCertificate.mul, coeff_mul,
    Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk, Finset.sum_range_succ]

/-- Full reversed polynomial recurrence corresponding to the finite evaluator.
Its interpretation as an independence polynomial of a graph is separate. -/
noncomputable def reversedProfile : Forest → ℕ[X] × ℕ[X]
  | .nil => (1,1)
  | .cons children rest =>
      let (pc,qc) := reversedProfile children
      let (pr,qr) := reversedProfile rest
      let rootOut := pc + X*pc
      ((rootOut+qc)*pr, rootOut*qr)

theorem profile_eq_firstFive (f : Forest) :
    profile f = (firstFive (reversedProfile f).1, firstFive (reversedProfile f).2) := by
  induction f with
  | nil => simp [profile, reversedProfile, firstFive, coeff_one]
  | cons children rest ihc ihr =>
      simp only [profile, reversedProfile, ihc, ihr]
      simp only [firstFive_mul, firstFive_add]
      simp [firstFive, coeff_X, shift, RootedForestCertificate.mul]

#print axioms concentration_coeff_le
#print axioms profile_eq_firstFive

end DepthThree.PolynomialCertificate
