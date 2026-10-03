"""Generator for the cell files PBScalar/PaperCells/Cells*.lean (see regen_paper_cells.sh).

Usage: python3 scripts/gen_paper_cells.py m1 m2 ...   (cells m in 4..15)
"""
import json,sys
d=json.load(open('data/universal_pb_finite_bernstein_full_certificate_2026-07-16.json'))
cells=d['payload']['cells']
which=[int(x) for x in sys.argv[1:]] if len(sys.argv)>1 else list(range(4,16))
def q(s):
    if '/' in s:
        p,r=s.split('/'); return f"({p}/{r} : ℝ)"
    return f"({s} : ℝ)"
from math import comb
SEP=" ∧\n    "
out=[]
for c in cells[1:]:
    m=c['left_endpoint']
    if m not in which: continue
    dg=c['degree']; a=c['numerator_power_coefficients_low_to_high']; b=c['bernstein_coefficients_low_to_high']
    poly=" +\n        ".join(f"{q(x)} * t ^ {k}" for k,x in enumerate(a))
    pos=SEP.join(f"0 < {q(x)}" for x in b)
    bern=" +\n        ".join(f"{q(x)} *\n          {comb(dg,i)} * t ^ {i} * (1 - t) ^ {dg-i}" for i,x in enumerate(b))
    out.append(f"""
/-- Cell `[{m}, {m+1}]`: numerator identity `4 H^{2*m} (S_{m} T_{m} - Q(H)) = P_{m}(H)` at `H = {m} + t`
(coefficients from the certificate file). -/
theorem cell{m}_numerator (t : ℝ) (ht : 0 ≤ t) :
    4 * ({m} + t) ^ {2*m} * (Sm ({m} + t) {m} * Tm ({m} + t) {m} - Qf ({m} + t)) =
      {poly} := by
  have hH : ({m} : ℝ) + t ≠ 0 := by positivity
  simp only [Sm, Tm, bH, Qf, Finset.sum_range_succ, Finset.prod_range_succ, Finset.sum_range_zero,
    Finset.prod_range_zero]
  field_simp
  ring

/-- Cell `[{m}, {m+1}]`: Bernstein expansion of `P_{m}({m} + t)` (coefficients from the certificate). -/
theorem cell{m}_bernstein (t : ℝ) :
    {poly} =
      {bern} := by
  ring

/-- Cell `[{m}, {m+1}]`: all Bernstein coefficients are positive. -/
theorem cell{m}_coeffs_pos :
    {pos} := by
  norm_num

set_option maxRecDepth 20000 in
/-- Cell `[{m}, {m+1}]` by the paper's route: `Q(H) ≤ S_{m} T_{m}` for `H ∈ [{m}, {m+1}]`. -/
theorem cell{m}_paper (t : ℝ) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    Qf ({m} + t) ≤ Sm ({m} + t) {m} * Tm ({m} + t) {m} := by
  have h := cell{m}_numerator t ht0
  rw [cell{m}_bernstein] at h
  have hs : 0 ≤ 1 - t := by linarith
  have hpos : 0 ≤ 4 * ({m} + t) ^ {2*m} * (Sm ({m} + t) {m} * Tm ({m} + t) {m} - Qf ({m} + t)) := by
    rw [h]; generalize 1 - t = s at hs ⊢; positivity
  have hH : (0 : ℝ) < 4 * ({m} + t) ^ {2*m} := by positivity
  have := (mul_nonneg_iff_of_pos_left hH).mp hpos
  linarith
""")
print("\n".join(out))
