# Log-concavity of the independence polynomials of double brooms

## Status and scope

A **double broom** `T_{p,q,t}` consists of two vertices joined by a path with
`t` internal vertices, with `p` leaves attached to one endpoint and `q` leaves
attached to the other. The parameters `p`, `q`, and `t` may be zero.

> **Theorem.** For every `p,q,t >= 0`, the independence polynomial of
> `T_{p,q,t}` is log-concave, and hence unimodal.

This is an all-parameter algebraic result, not a finite search. It covers
arbitrary subdivision of the connector between the two leaf-stars, but not
subdivision of their pendant arms. It therefore does **not** settle the full
class of trees with at most two branch vertices or Erdős Problem #993.

This repository note arose through substantive generative-AI assistance. The
formulas and finite identities have been replayed in exact arithmetic, but the
argument has not received independent human peer review.

## Connector recurrence

Write `F_t=I(T_{p,q,t};x)`, and let `P_j=I(P_j;x)` be the independence
polynomial of the path on `j` vertices, with the convenient conventions
`P_{-2}=0` and `P_{-1}=P_0=1`. Conditioning on the two endpoint vertices gives

$$
F_t=(1+x)^{p+q}P_t
+x\bigl((1+x)^p+(1+x)^q\bigr)P_{t-1}
+x^2P_{t-2}. \tag{1}
$$

The path recurrence immediately yields

$$
F_t=F_{t-1}+xF_{t-2}\qquad(t\geq2). \tag{2}
$$

## Partial-synchronization induction

For coefficient sequences `A=(a_i)` and `B=(b_i)`, write `A ~_p B` when

$$
a_mb_n+a_nb_m\geq
a_{m+1}b_{n-1}+a_{n-1}b_{m+1}\qquad(m\geq n). \tag{3}
$$

Hu, Wang, Zhao, and Zhao prove that partial synchronization is symmetric,
preserved by common convolution, and compatible with nonnegative sums of
pairwise partially synchronized log-concave sequences. In particular, the
sum of two partially synchronized log-concave sequences is log-concave.

Suppose the two base relations

$$
F_1\sim_p F_0,
\qquad
F_1\sim_p xF_0 \tag{4}
$$

hold. Induction on (2) then gives

$$
F_t\sim_p F_{t-1},
\qquad
F_t\sim_p xF_{t-1}\qquad(t\geq1), \tag{5}
$$

and makes every `F_t` log-concave. Indeed, the two summands in (2) are
partially synchronized. For the first relation at the next step, add
`F_{t-1} ~_p F_{t-1}` to `xF_{t-2} ~_p F_{t-1}`. For the second, add
`F_{t-1} ~_p xF_{t-1}` to the common-`x` convolution of
`F_{t-2} ~_p F_{t-1}`. It remains to prove (4).

## Reduction to two base polynomials

Put `s=p+q`, `A=1+x`, and

$$
B=A^s+x(A^p+A^q)=F_0,
\qquad
C=A^s+x.
$$

The subdivision identity gives

$$
F_1=B+xC. \tag{6}
$$

The polynomial `B` is log-concave. This is the known double-star case; for a
short direct verification, assume `p <= q`, put `d=q-p`, and factor

$$
B=A^p\bigl(A^q+xA^d+x\bigr).
$$

The polynomial `A^q+xA^d=A^d(A^p+x)` is log-concave by convolution closure.
Adding the final `x` can make only the index-2 Turán inequality harder. That
gap is

$$
\left(\binom q2+d\right)^2
-(q+2)\left(\binom q3+\binom d2\right).
$$

It is a concave quadratic in `d` on `0 <= d <= q`, so its minimum is at an
endpoint; after multiplication by 12, the endpoint values are

$$
q(q-1)(q^2-3q+8)
\quad\text{and}\quad
q(q+1)(q^2+q+4),
$$

both nonnegative. Thus `B` is log-concave. The polynomial `C` is also
log-concave: only its coefficient `c_1=s+1` differs from the binomial
sequence, and the only nontrivial new check reduces to

$$
3s(s-1)\geq2(s+1)(s-2),
$$

equivalently `s^2-s+4 >= 0`.

We now prove two stronger facts:

1. `B` and `C` are ordinarily synchronized.
2. `B` and `xC` are partially synchronized.

They imply (4). Log-concavity gives `B ~_p B` and `B ~_p xB`. The second
fact and (6) give `F_1 ~_p B`. The first fact implies `C ~_p B`; common
convolution by `x` gives `xC ~_p xB`, and hence `F_1 ~_p xB`.

## Ordinary synchronization of `B` and `C`

Write `b_k=[x^k]B` and `c_k=[x^k]C`, with binomial coefficients outside
their support interpreted as zero. Ordinary synchronization is equivalent to

$$
U_k:=c_kb_k-c_{k-1}b_{k+1}\geq0,
\qquad
V_k:=c_kb_k-c_{k+1}b_{k-1}\geq0. \tag{7}
$$

At the first two indices,

$$
\begin{aligned}
U_1&=\frac{(s+1)(s+4)}2,&
V_1&=\frac{s^2+7s+4}2,\\
U_2&=\frac{(s+1)\bigl(s(s-1)(s-2)+12pq\bigr)}{12},&
V_2&=\frac{s(s-1)(s^2+3s+8)}{12},
\end{aligned} \tag{8}
$$

which are nonnegative. For `k >= 3`, put `S_j=binom(s,j)`. Then

$$
U_k=S_k^2-S_{k-1}S_{k+1}
+\sum_{r\in\{p,q\}}
\left(S_k\binom r{k-1}-S_{k-1}\binom rk\right). \tag{9}
$$

The first term is nonnegative by binomial log-concavity, and each remaining
term is nonnegative because

$$
\frac{\binom rk}{\binom r{k-1}}
=\frac{r-k+1}{k}
\leq\frac{s-k+1}{k}
=\frac{S_k}{S_{k-1}}.
$$

At `k=3`, direct symmetric simplification gives

$$
V_3=
\frac{s(s-1)(s-2)\bigl(s^3+6s^2+5s-24pq\bigr)}{144}\geq0, \tag{10}
$$

because `pq <= s^2/4`.

For `k >= 4`, expand `V_k` in the nonnegative product-binomial basis
`binom(p,i)binom(q,j)`. Vandermonde's identity gives

$$
\binom{p+q}{a}\binom pb
=\sum_{i,j\geq0}
\binom ib\binom b{a+b-i-j}
\binom pi\binom qj, \tag{11}
$$

and

$$
\binom{p+q}{a}\binom{p+q}{b}
=\sum_{i,j\geq0}
\binom{i+j}{a}\binom a{a+b-i-j}
\binom pi\binom qj. \tag{12}
$$

Only indices with `k <= r:=i+j <= 2k` occur. Put `h=2k-r`. The coefficient
of `binom(p,i)binom(q,j)` is

$$
\begin{aligned}
\mathcal V_{i,j}^{(k)}={}&
\frac{\binom rk\binom kh}{k-h+1}\\
&+\frac{\binom{k-1}{h-1}}{k-1}
\left[
\binom i{k-2}(2-j)+\binom j{k-2}(2-i)
\right].
\end{aligned} \tag{13}
$$

The magnitude of the negative part in brackets is at most

$$
(j-2)_+\binom i{k-2}+(i-2)_+\binom j{k-2}
\leq2\binom rk, \tag{14}
$$

by the corresponding Vandermonde terms. Therefore

$$
\mathcal V_{i,j}^{(k)}\geq
\binom rk\left(
\frac{\binom kh}{k-h+1}
-\frac{2\binom{k-1}{h-1}}{k-1}
\right)\geq0. \tag{15}
$$

The last inequality is equivalent to
`k(k-1) >= 2h(k-h+1)`. It is direct for `k=4`; for `k >= 5`, it follows from
`h(k-h+1) <= (k+1)^2/4`. This proves ordinary synchronization.

## Partial synchronization of `B` and `xC`

First consider the diagonal inequalities

$$
W_k:=2c_{k-1}b_k-c_{k-2}b_{k+1}-c_kb_{k-1}\geq0. \tag{16}
$$

The exceptional low indices are

$$
\begin{aligned}
W_1&=s+3,\\
W_2&=\frac{2s^3+9s^2+13s+6pq}{6},\\
W_3&=\frac{s^5+s^4+s^3-s^2-2s-12pq(s^2-s+2)}{24}.
\end{aligned} \tag{17}
$$

The first two are immediate. For `s >= 2`, the bound `pq <= s^2/4` makes
the numerator of `W_3` at least
`s(s-2)(s^3+4s+1)`, while `s=0,1` can be checked directly.

For `k >= 4`, use the same product-binomial basis. Only
`k-1 <= r:=i+j <= 2k-1` occurs; put `h=2k-1-r`. The coefficient is

$$
\begin{aligned}
\mathcal W_{i,j}^{(k)}={}&
\frac{2\binom r{k-1}\binom{k-1}h}{k-h+1}\\
&+\frac{\binom{k-1}{h-1}}{(k-1)(k-h+1)}
\Bigl[
\binom i{k-2}\bigl(2(i-k+2)-j(j-1)\bigr)\\
&\hspace{39mm}+
\binom j{k-2}\bigl(2(j-k+2)-i(i-1)\bigr)
\Bigr].
\end{aligned} \tag{18}
$$

The positive linear terms may be discarded. The magnitude of what remains is
bounded by

$$
j(j-1)\binom i{k-2}+i(i-1)\binom j{k-2}
\leq4\binom rk. \tag{19}
$$

Comparison with the first line of (18) reduces nonnegativity to
`k(k-1) >= 2h`, which holds for `k >= 4` and `0 <= h <= k`. Thus all the
diagonal inequalities hold.

To obtain full partial synchronization, set
`alpha_k=c_k/c_{k-1}` and `beta_k=b_k/b_{k-1}` on the positive supports.
Ordinary synchronization gives

$$
\beta_{k+1}\leq\alpha_k,
\qquad
\alpha_{k+1}\leq\beta_k. \tag{20}
$$

For `m >= n+1`, log-concavity and (20) imply

$$
c_{m-1}b_n\geq c_mb_{n-1},
\qquad
c_{n-1}b_m\geq c_{n-2}b_{m+1}. \tag{21}
$$

Their sum is (3) for `xC` and `B`; when `m=n`, (3) is exactly (16). Outside
the positive supports it is automatic. Hence `xC ~_p B`, completing the base
relations and the proof of the theorem.

## Exact replay

The main verifier independently checks the graph-DP formula, recurrence, base
partial-synchronization relations, product-binomial identities, and their
nonnegative coefficient bounds. It uses exact integer and rational arithmetic:

```text
python3 verify_connector_partial_sync_route_20260808.py
```

Its checked result is
[`results/connector_partial_sync_route_20260808.json`](../results/connector_partial_sync_route_20260808.json).
The finite ranges audit the algebra; they are not the proof of the universal
statement.

## Bounded two-hub companion

A separate exact enumeration extends the recurrence argument beyond double
brooms. Attach pendant paths with positive lengths
`a=(a_1,...,a_r)` and `b=(b_1,...,b_s)` to the two connector endpoints.
If the total pendant weight satisfies

$$
\sum_i a_i+\sum_j b_j\leq24,
$$

then the resulting tree has a log-concave independence polynomial for every
connector length. The connector remains unbounded; only the pendant core is
bounded.

The certificate enumerates all 163,523 unordered pendant-core pairs through
weight 24 and checks four exact base conditions for each pair. Its canonical
enumeration digest is

```text
25c37f1435984d3bbdb275424ea4f6f5766c33a7c6e792a67700a90dae91efb1
```

Replay it with

```text
python3 verify_c2_bounded_pendant_core_20260808.py
python3 -m unittest test_c2_bounded_pendant_core.py
```

The checked report is
[`results/c2_bounded_pendant_core_20260808.json`](../results/c2_bounded_pendant_core_20260808.json).
This is a computer-assisted bounded theorem, not evidence for pendant weight
above 24.

## Prior-art boundary

Galvin and Hilyard prove log-concavity for the unsubdivided double-star case as
part of a stronger alternating-caterpillar theorem. Their stated family does
not include the general endpoint-only bristle sequence `(p,0,...,0,q)` treated
here. A targeted search found no exact prior statement of the all-connector
result, but this has not been authenticated through MathSciNet or zbMATH.

## References

- D. Galvin and J. Hilyard, “The independent set sequence of some families of
  trees,” *Australasian Journal of Combinatorics* 70 (2018), 236–252.
  [PDF](https://academicweb.nd.edu/~dgalvin1/pdf/journal/centipedes_J.pdf)
- H. Hu, D. G. L. Wang, F. Zhao, and T. Zhao, “Convolution preserves partial
  synchronicity of log-concave sequences,” *Mathematical Inequalities &
  Applications* 20 (2017), 91–103.
  [DOI](https://doi.org/10.7153/mia-20-07)
