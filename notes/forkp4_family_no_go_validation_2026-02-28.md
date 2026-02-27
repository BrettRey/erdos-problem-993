# Fork2/P4 Family No-Go Validation (2026-02-28)

## Scope
Validated the new family claim:
- core `T0` + `(X/Y)` kernel + `(Fork2/P4)` scaffolds
- symbolic lock identity at degree 5:
  - `9*i4 - 7*i5 = -S/120`
  - `S` coefficientwise positive, no constant term.

## Source artifacts (downloaded)
- `/Users/brettreynolds/Downloads/family_attack_T0_XY_ForkP4_no_go.py`
- `/Users/brettreynolds/Downloads/family_attack_T0_XY_ForkP4_no_go_identity.txt`
- `/Users/brettreynolds/Downloads/family_attack_T0_XY_ForkP4_no_go_identity_terms.json`

## In-repo reproducible script
- `scripts/family_attack_t0_xy_forkp4_no_go.py`

## Reproduction command
This machine does not expose `sympy` globally, so run in a local venv:

```bash
python3 -m venv .venv_sympy
. .venv_sympy/bin/activate
python -m pip install sympy
python scripts/family_attack_t0_xy_forkp4_no_go.py \
  --out-txt results/family_attack_t0_xy_forkp4_no_go_identity.txt \
  --out-json results/family_attack_t0_xy_forkp4_no_go_identity_terms.json
```

## Produced local artifacts
- `results/family_attack_t0_xy_forkp4_no_go_identity.txt`
- `results/family_attack_t0_xy_forkp4_no_go_identity_terms.json`

## Equality check (downloaded vs local)
Term-set comparison by exponent tuple `(a,b,c,d)`:
- term count: `125` vs `125`
- exact coefficient match: `true`

Local summary:
- `num_terms = 125`
- `min_coeff = 7168`
- `all_positive_nonconstant = true`
- `has_constant_term = false`

## Consequence
Within this family, any nontrivial attachment `(a,b,c,d) != (0,0,0,0)` yields:
- `9*i4 - 7*i5 < 0`
- therefore `i4/i5 < 7/9`
- so lock `(m,lambda) = (5,7/9)` is impossible outside the trivial base point.

This validates the family-level no-go mechanism, but remains family-specific (not a global injectivity proof).
