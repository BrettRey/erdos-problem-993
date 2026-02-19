# Sub-claim A helper lower bound: global status (2026-02-18)

Helper inequality used in the lane-2 route:

`lb_u2(k,j) := mu_{k,j+2}(u2_{k,j}) - mu_{k,j}(lambda_{k,j}) >= 1`,

where

- `lambda_{k,j}` is the tie fugacity at mode `m_{k,j}`,
- `u2_{k,j} = a'_{m_{k,j}}/a'_{m_{k,j}+1}`,
- `a'_t = [x^t](1+2x)^k(1+x)^(j+2)`.

---

## 1) New checker

- `prove_subclaim_A_helper_lb.py`

This script checks:

1. mode-step prerequisite `m_{k,j+2}=m_{k,j}+1`,
2. helper inequality `lb_u2 >= 1`,
3. exact-rational values on the failure set.

---

## 2) Result: strict global form is false

Run A:

`python3 prove_subclaim_A_helper_lb.py --k-min 6 --k-max 1200 --j-max 120 --out results/whnc_subclaim_A_helper_lb_k6_1200_j120.json`

Run B:

`python3 prove_subclaim_A_helper_lb.py --k-min 6 --k-max 4000 --j-max 80 --exact-cap 10 --out results/whnc_subclaim_A_helper_lb_k6_4000_j80.json`

Observed in both runs:

- mode-step failures: `0`,
- helper failures: exactly `5`.

Failure set:

- `(k,j)=(6,0)`,
- `(k,j)=(6,2)`,
- `(k,j)=(7,1)`,
- `(k,j)=(8,0)`,
- `(k,j)=(8,2)`.

Exact-rational checks confirm all five are genuinely `<1`:

- `0.979591142954440...`,
- `0.994381961240462...`,
- `0.990794685234226...`,
- `0.986831750231117...`,
- `0.999457876207099...`.

So the strict statement "`lb_u2>=1` for all `k>=6,j>=0`" cannot be proved as-is.

---

## 3) Corrected global candidate

In the same scans, the corrected form holds:

- `lb_u2(k,j) >= 1` for all checked `k>=9`, `j>=0`.

Min witness in run B (`k<=4000,j<=80`):

- `min lb_u2 (k>=9) = 1.000089495361863` at `(k,j)=(3998,80)`.

Additional larger run (same day):

`python3 prove_subclaim_A_helper_lb.py --k-min 9 --k-max 4000 --j-max 120 --exact-cap 0 --out results/whnc_subclaim_A_helper_lb_k9_4000_j120.json`

Observed:

- mode-step failures: `0`,
- helper failures: `0`,
- `min lb_u2 (k>=9) = 1.000087559139956` at `(k,j)=(3998,120)`.

This is compatible with a finite-base correction:

- treat the five `(k,j)` exceptions explicitly,
- apply the helper bound on all remaining pairs.

---

## 4) Relation to Sub-claim A itself

Even on the five helper-failure pairs, the actual target mean-step remains `>1`
(checked exactly in `prove_subclaim_A_meanstep_lane2.py`).

So what fails is only the surrogate helper inequality, not Sub-claim A’s target
step inequality.
