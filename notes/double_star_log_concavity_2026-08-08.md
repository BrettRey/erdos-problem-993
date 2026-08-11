# Double stars are log-concave — 2026-08-08

## Result

Let `D_{p,q}` be the tree formed from an edge by adjoining `p` leaves to one
endpoint and `q` leaves to the other, where `p,q >= 0`. Its independence
polynomial is log-concave, and hence unimodal.

This is a self-contained reproof of a known base lemma for the two-hub (`C_2`)
program. Galvin--Hilyard's Theorem 1.9 already proves all double stars
log-concave as the `n=1` case of a stronger alternating-caterpillar theorem.
The argument below supplies a short hand proof of the base inequality that
their paper discharged computationally. It does not prove the full `C_2`
conjecture: subdividing pendant arms or the path between the hubs is the
remaining transfer problem.

## Proof

Conditioning on the two adjacent hubs gives

$$
I(D_{p,q};x)=(1+x)^{p+q}+x(1+x)^p+x(1+x)^q. \tag{1}
$$

By symmetry, assume `p <= q` and put `d=q-p`. Factor (1) as

$$
I(D_{p,q};x)=(1+x)^p H(x),
$$

where

$$
\begin{aligned}
H(x)&=G(x)+x,\\
G(x)&=(1+x)^q+x(1+x)^d\\
    &=(1+x)^d\bigl((1+x)^p+x\bigr).
\end{aligned}
$$

The coefficient sequence of the star polynomial
`(1+x)^p+x` is log-concave. Away from indices 1 and 2 this is ordinary
binomial log-concavity; the only exceptional checks are

$$
(p+1)^2\geq \binom p2,
$$

and

$$
\binom p2^2-(p+1)\binom p3
=\frac{p(p-1)(p^2-p+4)}{12}\geq0.
$$

The standard convolution closure theorem for nonnegative log-concave
sequences without internal zeros therefore makes `G` log-concave.

Passing from `G` to `H=G+x` increases only the coefficient of `x`. This
improves the Turan inequality at index 1, leaves every inequality from index
3 onward unchanged, and can only make the index-2 inequality harder. The
three relevant coefficients are

$$
h_1=q+2,\qquad
h_2=\binom q2+d,\qquad
h_3=\binom q3+\binom d2.
$$

Thus it remains to prove

$$
\Delta_q(d):=
\left(\binom q2+d\right)^2
-(q+2)\left(\binom q3+\binom d2\right)\geq0
\quad(0\leq d\leq q). \tag{2}
$$

For fixed `q`, `Delta_q(d)` is a concave quadratic in `d`, with quadratic
coefficient `-q/2`. Its minimum on `[0,q]` is therefore attained at an
endpoint. Direct simplification gives

$$
\Delta_q(0)=
\frac{q(q-1)(q^2-3q+8)}{12}\geq0
$$

and

$$
\Delta_q(q)=
\frac{q(q+1)(q^2+q+4)}{12}\geq0.
$$

This proves (2), hence `H` is log-concave. A final application of convolution
closure to `(1+x)^p H(x)` proves the result.

## Exact replay

The proof is algebraic rather than computational. The replay script checks
the graph-DP formula, the factorization, both endpoint identities, and a
finite parameter box in exact integer arithmetic:

```text
python3 verify_double_star_log_concavity_20260808.py
python3 test_double_star_log_concavity.py
```

## Boundary of the result

The earlier synchronization obstructions remain valid: the natural summands
in (1), and the two summands in the subdivision-contraction identity, need not
be synchronized. The proof succeeds because it retains the aggregate Turan
gap after a different factorization. The next nontrivial target is a transfer
lemma for subdividing a double star, not another pairwise-synchronization
criterion.

Prior-art source: D. Galvin and J. Hilyard, *The independent set sequence of
some families of trees*, Australasian Journal of Combinatorics 70 (2018),
236--252, Theorem 1.9(1).
