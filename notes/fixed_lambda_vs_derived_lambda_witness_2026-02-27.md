# Fixed-`lambda` vs Derived-Canonical-`lambda` Check (2026-02-27)

C supplied a realizable pair with same `rho` at fixed `lambda=7/9`:
- `g6 = JgcO_OC_??_`
- `g6 = LgcOOC@@??o??@`

Both are valid min-`u` canonical gated trees with canonical triplet `(2,1,0)`.

## Verified canonical records

Using `build_record_from_triplet` under min-`u`:

1) `JgcO_OC_??_` (`n=11`)
- `m=3`
- derived canonical `lambda = i_{m-1}/i_m = 15/29`
- `N=[x]P=8`
- `P=[1,8,22,25,10]`
- `Q=[0,1,6,12,10,3]`
- `I=[1,11,45,87,82,33,3]`

2) `LgcOOC@@??o??@` (`n=13`)
- `m=4`
- derived canonical `lambda = i_{m-1}/i_m = 167/221`
- `N=[x]P=10`
- `P=[1,10,37,62,45,10]`
- `Q=[0,1,8,23,29,16,3]`
- `I=[1,13,66,167,221,145,39,3]`

## Fixed-`lambda` check
At externally fixed `lambda=7/9`, both satisfy
- `rho(7/9) = 17920/44229`.

## Conclusion
This is **not** a canonical `(m,lambda,rho)` split witness because the derived canonical
`lambda` values differ (`15/29` vs `167/221`) and `m` differs (`3` vs `4`).

It is instead a valid fixed-activity collision example:
- same `Q(7/9)/P(7/9)`, different `N`.

So it does not refute injectivity of the canonical key `(m,lambda,rho) -> N`.
