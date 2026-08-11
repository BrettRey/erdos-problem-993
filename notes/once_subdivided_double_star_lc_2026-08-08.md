# Once-subdivided double stars are log-concave — 2026-08-08

## Result

Let `S_{p,q}` be obtained from the double star `D_{p,q}` by subdividing its
hub edge exactly once. For every `p,q >= 0`, the independence polynomial of
`S_{p,q}` is log-concave, and hence unimodal.

This is the first connector-transfer step beyond the unsubdivided theorem in
`notes/double_star_log_concavity_2026-08-08.md`. It does not prove the result
for two or more connector subdivisions, nor for subdivided pendant arms.

## Polynomial and log-concave core

Conditioning on the two hubs and the intervening vertex gives

$$
F_{p,q}(x)=(1+x)^{p+q+1}+x(1+x)^p+x(1+x)^q+x^2. \tag{1}
$$

Assume by symmetry that `p <= q`, and put `d=q-p`. Removing the final `x^2`
term leaves

$$
\begin{aligned}
K_{p,q}(x)
&=F_{p,q}(x)-x^2\\
&=(1+x)^p
  \left((1+x)^{q+1}+x(1+x)^d+x\right). \tag{2}
\end{aligned}
$$

The bracket in (2) is exactly the auxiliary polynomial proved log-concave in
the double-star argument, with its parameter `q` there replaced by `q+1`.
Multiplication by `(1+x)^p` preserves log-concavity by convolution closure.
Thus `K_{p,q}` is log-concave.

## The two exceptional inequalities

Restoring `x^2` increases only coefficient 2. It improves the Turan inequality
at index 2, leaves every inequality from index 4 onward unchanged, and can
only make the inequalities at indices 1 and 3 harder.

Put `s=p+q`. At index 1, the required inequality is

$$
(s+3)^2\geq\binom{s+2}{2},
$$

whose doubled gap is `s^2+9s+16 > 0`.

For index 3, write `f_k=[x^k]F_{p,q}` and
`z=(p-q)^2`. Direct symmetric expansion gives

$$
\begin{aligned}
144(f_3^2-f_2f_4)
={}&s(s^5+6s^4+s^3-9s^2+16s-60)\\
&+3\bigl(s(s-1)(s+4)+12\bigr)z+9z^2. \tag{3}
\end{aligned}
$$

For `s=0`, both sides of the desired inequality vanish. For `s=1`, the
integer constraints force `z=1`, and the right side of (3) again vanishes.
For `s >= 2`, set `u=s-2`. The only factor in (3) not manifestly
nonnegative satisfies

$$
s^5+6s^4+s^3-9s^2+16s-60
=u^5+16u^4+89u^3+221u^2+264u+72>0.
$$

The coefficient of `z` is positive because
`s(s-1)(s+4)+12 >= 12`. Hence (3) is nonnegative for all integer
`p,q >= 0`. The two exceptional inequalities hold, proving that (1) is
log-concave.

## Exact replay

The proof is algebraic. The verifier independently checks the tree-DP
formula, factorization (2), symmetric identity (3), and a finite parameter
box using exact integers:

```text
python3 verify_once_subdivided_double_star_lc_20260808.py
python3 test_once_subdivided_double_star_lc.py
```

## Remaining frontier

For `t >= 0` internal connector vertices, with path polynomials `P_t` and the
conventions `P_{-2}=0`, `P_{-1}=P_0=1`, the full connector family has

$$
(1+x)^{p+q}P_t
+x\bigl((1+x)^p+(1+x)^q\bigr)P_{t-1}
+x^2P_{t-2}. \tag{4}
$$

The present sparse-perturbation proof covers `t=0` and `t=1`. The subsequent
product-binomial-basis proof in
`notes/connector_partial_sync_route_2026-08-08.md` establishes the two base
partial-synchronization relations and iterates formula (4) to every `t`.
