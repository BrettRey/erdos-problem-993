# Summary of changes for run e6301287-ecb0-4d75-8b65-b6711945038f
Formalized `pascal_smoothing_shadow_lemma_2026-08-20.md` in Lean 4 / Mathlib, in the new file `RequestProject/PascalSmoothing.lean`. The file builds with no `sorry` and the main results depend only on the standard axioms (`propext`, `Classical.choice`, `Quot.sound`).

Setup (definitions):
- `faceCount Δ j` — the number of `j`-faces `a_j` of a family `Δ` of subsets of `Fin M`;
- `erasureShadow Δ d` — the profile `E_d = ∑_j a_j C(M-j, d)` of equation (1);
- `layerRatio` / `pascalTransform` — the normalized face numbers `b_j = a_j / C(M,j)` and their Pascal transform `S_m = ∑_j C(m,j) b_j`;
- `shadowGenFun` — the generating polynomial `F(t) = (1+t)^M A_Δ(t/(1+t)) = ∑_j a_j t^j (1+t)^{M-j}`.

Results proved (downward-closedness of `Δ` is the only structural hypothesis; no graph, poset or forest assumption is used):
- `pascal_smoothing_shadow_lemma` — the Pascal-smoothing lemma (2): for `1 ≤ d ≤ M-1` and `m = M-d`, `8(d+1)(m+1) E_{d-1} E_{d+1} ≤ 9 d m E_d²` (stated over ℕ with denominators cleared, which is equivalent to `E_d²/(E_{d-1}E_{d+1}) ≥ (8/9)(d+1)(m+1)/(dm)`);
- `erasureShadow_log_concave_of` — the sufficient condition (8): log-concavity `E_{d-1}E_{d+1} ≤ E_d²` whenever `d(M-d) ≤ 8(M+1)`;
- `erasureShadow_log_concave_depth_le_eight` — consequence 1: log-concavity through defect depth eight (`d ≤ 8`);
- `erasureShadow_log_concave_of_M_le_33` — consequence 2: log-concavity of the whole sequence when `M ≤ 33`;
- `erasureShadow_depth_three` and `erasureShadow_depth_three_strict` — consequence (3): `32(M-2) E_2 E_4 ≤ 27(M-3) E_3²` for `M ≥ 4`, and the strict form `E_2 E_4 < E_3²` under the (necessary) assumption `E_2 E_4 > 0`;
- supporting results matching the write-up: the local-LYM monotonicity `(j+1)a_{j+1} ≤ (M-j)a_j` and `b_{j+1} ≤ b_j` (4), the absorption identity `E_d = C(M,d) S_{M-d}` (5), the shifted-transform Pascal recurrences (6), the core inequality `8 S_{m-1}S_{m+1} ≤ 9 S_m²` (7), and the generating-function reading `E_d = [t^{M-d}] F(t)`.

Two notes on faithfulness: the nonemptiness hypothesis of the source statement is kept in the user-facing theorems though the inequality does not require it (documented in the docstring), and the strict depth-three inequality needs the extra positivity hypothesis since equality `0 = 0` occurs when `E_2 E_4 = 0`. All work is committed and pushed.