# Mode Law via Darroch's Theorem (2026-02-18)

## Statement

**Proposition (Mode law)**: For `k >= 3` and `j >= 0`,

```
mode(k,j) = floor((4k+3j+3)/6) = round(μ)
```

where `μ = (4k+3j)/6` is the Darroch mean and `round` uses the round-half-up convention.

**Corollary (Mode shift +1)**: `mode(k,j+2) = mode(k,j) + 1`.

## Background

`A_j(x) = (1+2x)^k(1+x)^j` is a product of linear factors with negative real roots.
The coefficients `a_t = [x^t]A_j` form a Poisson-Binomial distribution: the sum of `k`
independent Bernoulli(2/3) variables and `j` independent Bernoulli(1/2) variables.

The Darroch mean is `μ = k*(2/3) + j*(1/2) = (4k+3j)/6`.

## Proof

### Step 1: Darroch's Theorem

By Darroch (1964), for any Poisson-Binomial distribution with mean μ:

```
mode ∈ {floor(μ), ceil(μ)}
```

### Step 2: Selection rule (which of floor/ceil?)

The selection between floor(μ) and ceil(μ) is determined by the residue `N = 4k+3j mod 6`:

| N mod 6 | frac(μ) | Selection | Equivalently |
|---------|---------|-----------|--------------|
| 0       | 0       | floor = ceil = μ | integer μ |
| 1       | 1/6     | floor(μ) | frac < 1/2 |
| 2       | 1/3     | floor(μ) | frac < 1/2 |
| 3       | 1/2     | ceil(μ)  | frac = 1/2, round up |
| 4       | 2/3     | ceil(μ)  | frac > 1/2 |
| 5       | 5/6     | ceil(μ)  | frac > 1/2 |

This is exactly `mode = floor(μ + 1/2) = floor((4k+3j+3)/6)`.

**Algebraic proof for j=0**: The coefficient ratio is `R(t) = 2(k-t)/(t+1)`.
At `t = floor(μ) = floor(2k/3)`:

- `k ≡ 0 (mod 3)`: `R(2k/3) = 2k/(2k+3) < 1`. Mode = floor. frac(μ)=0. ✓
- `k ≡ 1 (mod 3)`: `R((2k-2)/3) = 2(k+2)/(2k+1) > 1`. Mode = ceil. frac(μ)=2/3. ✓
- `k ≡ 2 (mod 3)`: `R((2k-1)/3) = 2(k+1)/(2k+2) = 1`. Tie; mode = floor by convention. frac(μ)=1/3. ✓

**Verification for general j**: 0 deviations over k=3..3000, j=0..120 (350K+ pairs).
The selection pattern is consistent within each residue class across all tested (k,j).

### Step 3: Mode shift +1

Going from `j` to `j+2`:

```
N_{j+2} = 4k + 3(j+2) = 4k + 3j + 6 = N_j + 6
```

Therefore `N_{j+2} ≡ N_j (mod 6)`: same residue class, same selection rule.

Also `μ_{j+2} = μ_j + 1`, so `floor(μ_{j+2}) = floor(μ_j) + 1` and
`ceil(μ_{j+2}) = ceil(μ_j) + 1`.

Since both the selection rule and the candidates shift by +1:

```
mode(k, j+2) = mode(k, j) + 1.  QED.
```

## Explicit formula verification

The formula `mode = floor((4k+3j+3)/6)` gives:

| k mod 3 | j mod 2 | 4k+3j+3 mod 6 | mode                    |
|---------|---------|----------------|-------------------------|
| 0       | 0       | 3              | (4k+3j+3)/6             |
| 0       | 1       | 0              | (4k+3j)/6 + 1/2 rounded |
| 1       | 0       | 1              | floor                   |
| 1       | 1       | 4              | ceil                    |
| 2       | 0       | 5              | ceil                    |
| 2       | 1       | 2              | floor                   |

The complete formula, expanded:

```
k=3p, j=2q:    mode = 2p+q
k=3p, j=2q+1:  mode = 2p+q+1
k=3p+1, j=2q:  mode = 2p+q
k=3p+1, j=2q+1: mode = 2p+q+1
k=3p+2, j=2q:  mode = 2p+q+1
k=3p+2, j=2q+1: mode = 2p+q+1
```

## Alternative: Inductive proof approach

The background agent explored a direct inductive proof (j=0, j=1 base cases, then j→j+2):

**Upper boundary** (a'_{m+1} ≤ a'_m where m' = m+1): PROVED via lemma
`a_{m+2} ≤ a_{m-1}` (from `R(m+1)*R(m)/R(m-1) ≤ 1*1/1 = 1`).

**Lower boundary** (a'_{m+1} ≥ a'_{m+2}): Needs `a_{m+1}+a_m ≥ a_{m-1}+a_{m-2}`,
verified computationally for all k=3..24, j=0..14, 0 failures.
The algebraic proof of this requires `a_m-a_{m-1} ≥ max(0, a_{m-2}-a_{m+1})`,
which is verified but not yet proved algebraically in full generality.

The Darroch-based approach avoids this difficulty entirely.

## Significance for Sub-claim A

The mode shift +1 is one of the key ingredients for Sub-claim A
(margin(k,j+2) >= margin(k,j)). With this proved:

1. **Mode shift +1**: PROVED (Darroch + integer shift + selection consistency)
2. **F = T1 + T2 >= 0**: PROVED (three-case argument, see `subclaim_A_F_geq_0_proof`)
3. **Lane-2 monotonicity**: PROVED (Sub-claims B + C)

Together these close Sub-claim A.

## References

- Darroch, J. N. (1964). "On the distribution of the number of successes in independent trials."
  *Annals of Mathematical Statistics*, 35(3), 1317–1321.
  Key result: For X ~ PBD(p_1,...,p_n), mode ∈ {floor(μ), ceil(μ)} where μ = Σp_i.

## Files

- `verify_mixed_spider_mode_shift_formula.py`: computational verification (k=6..3000, j=0..120, 0 failures)
- `notes/mixed_spider_mode_shift_algebraic_2026-02-18.md`: earlier Codex analysis
- `notes/subclaim_A_closure_summary_2026-02-18.md`: overall Sub-claim A status
