import Mathlib
import RequestProject.Defs
import RequestProject.ZBasic
import RequestProject.CauchyBinet
import RequestProject.Conv
import RequestProject.Factors
import RequestProject.Sector
import RequestProject.Sharpness

open scoped BigOperators
open scoped Real
open scoped Nat
open scoped Classical
open scoped Pointwise

set_option maxHeartbeats 8000000
set_option maxRecDepth 4000
set_option synthInstance.maxHeartbeats 20000
set_option synthInstance.maxSize 128

set_option relaxedAutoImplicit false
set_option autoImplicit false

set_option pp.fullNames true
set_option pp.structureInstances true
set_option pp.coercions.types true
set_option pp.funBinderTypes true
set_option pp.letVarTypes true
set_option pp.piBinderTypes true

set_option grind.warning false

/-! ## Axiom audit for every named result of this project -/

#print axioms LC.inner_mul_ge
#print axioms LC.step_mul_ge
#print axioms LC.cauchyBinet_two
#print axioms LC.cauchyBinet_two_nonneg
#print axioms LC.logConcave_of_conv
#print axioms LC.noInternalZeros_of_conv
#print axioms LC.logConcave_conv
#print axioms LC.logConcave_coeff_mul
#print axioms LC.logConcave_linear
#print axioms LC.noInternalZeros_linear
#print axioms LC.logConcave_quadratic
#print axioms LC.nonneg_quadratic
#print axioms LC.noInternalZeros_quadratic
#print axioms LC.logConcave_quadratic_of_abs_le_pi_div_three
#print axioms LC.sector_aux
#print axioms LC.logConcave_of_roots_in_sector
#print axioms LC.conv_not_logConcave_without_noInternalZeros
