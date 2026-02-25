# Canonical Pair: Same P, Different Derived Lambda (2026-02-25)

## Purpose

This note records an explicit **canonical-valid** pair showing that canonical-derived
`lambda` is not determined by `P` alone (because `Q` can vary while `P` stays fixed).

That directly supports the bridge-2 obstruction:

```text
I(T;x) = (1+2x)P(x) + (1+x)Q(x)
```

and canonical `lambda = i_{m-1}/i_m` depends on `Q`.

## Gate-verified pair

Both trees pass:

- tree property
- `is_dleaf_le_1(...) == True`
- `bridge_decomposition(..., require_dleaf=True) != None`

### Tree A

- `g6 = H?AA@bg`
- canonical triplet from decomposition: `(leaf,support,u) = (3,7,2)`
- `P = [1, 6, 10, 5]`
- `Q = [0, 1, 5, 8, 4]`
- `I = [1, 9, 28, 38, 22, 4]`
- leftmost mode: `m = 3`
- derived `lambda = i_2/i_3 = 28/38 = 14/19`

### Tree B

- `g6 = H?ABAag`
- canonical triplet from decomposition: `(leaf,support,u) = (3,7,1)`
- `P = [1, 6, 10, 5]`
- `Q = [0, 1, 5, 6, 2]`
- `I = [1, 9, 28, 36, 18, 2]`
- leftmost mode: `m = 3`
- derived `lambda = i_2/i_3 = 28/36 = 7/9`

## Conclusion

`P` is identical across the pair, but canonical-derived `lambda` differs:

```text
14/19 != 7/9.
```

So preserving `P` (or any `P`-only jet/cumulant data) does not by itself preserve
canonical `(m,lambda)`. Any canonical lift from fixed-`lambda` no-go must include a
`Q`-sensitive control.
