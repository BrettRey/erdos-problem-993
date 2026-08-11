# All connector-subdivided double stars are log-concave — 2026-08-08

## Theorem

Let `T_{p,q,t}` be the tree with `p` leaves at one hub, `q` leaves at the
other, and `t` internal vertices on the path between the hubs. For every
`p,q,t >= 0`, its independence polynomial is log-concave, and hence unimodal.

This settles arbitrary subdivision of the hub connector when all pendant arms
have length one. It does not settle arbitrary subdivision of the pendant arms,
so it is not a proof of the full two-branch-vertex (`C_2`) conjecture or of
Erdos Problem 993.

## Connector recurrence

Write `F_t=I(T_{p,q,t};x)`. Exact conditioning on the hubs gives

$$
F_t=(1+x)^{p+q}P_t
+x\bigl((1+x)^p+(1+x)^q\bigr)P_{t-1}
+x^2P_{t-2}, \tag{1}
$$

where `P_j=I(P_j;x)` is the independence polynomial of the `j`-vertex path
and `P_{-2}=0`, `P_{-1}=P_0=1`. Hence

$$
F_t=F_{t-1}+xF_{t-2}\qquad(t\geq2). \tag{2}
$$

The cases `t=0` and `t=1` are log-concave by the elementary perturbation
proofs in the preceding two notes.

## Partial-synchronization induction

For coefficient sequences `A=(a_i)` and `B=(b_i)`, write `A ~_p B` when

$$
a_mb_n+a_nb_m\geq a_{m+1}b_{n-1}+a_{n-1}b_{m+1}
\quad(m\geq n). \tag{3}
$$

Hu--Wang--Zhao--Zhao prove that partial synchronization is symmetric,
preserved by common convolution, and preserved under nonnegative sums of
pairwise partially synchronized log-concave sequences. It also implies that
the sum of the two sequences is log-concave.

Suppose

$$
F_1\sim_p F_0
\quad\text{and}\quad
F_1\sim_p xF_0. \tag{4}
$$

Then induction on (2) gives, for every `t >= 1`,

$$
F_t\sim_p F_{t-1},
\qquad
F_t\sim_p xF_{t-1}, \tag{5}
$$

and makes every `F_t` log-concave. Indeed, if (5) holds at `t-1`, then the
two summands in (2) are partially synchronized, so their sum is log-concave.
For the first relation at `t`, add the relations
`F_{t-1} ~_p F_{t-1}` and `xF_{t-2} ~_p F_{t-1}`. For the second, add
`F_{t-1} ~_p xF_{t-1}` to the common-`x` convolution of
`F_{t-2} ~_p F_{t-1}`. Thus it remains only to prove (4).

## Reduction of the base relations

Put `s=p+q`, `A=1+x`, and define

$$
B=A^s+x(A^p+A^q)=F_0,
\qquad
C=A^s+x.
$$

The subdivision identity gives

$$
F_1=B+xC. \tag{6}
$$

We prove two claims:

1. `B` and `C` are ordinarily synchronized.
2. `B` and `xC` are partially synchronized.

Both sequences needed in Claim 1 are log-concave. This is already known for
`B=F_0`. For `C`, the only altered binomial coefficient is
`c_1=s+1`: the log-concavity inequality at index 1 is immediate, while the
one at index 2 reduces to

$$
3s(s-1)\geq 2(s+1)(s-2),
$$

or `s^2-s+4 >= 0`; all later inequalities are ordinary binomial
log-concavity.

The claims imply (4). Since `B` is log-concave, it is partially synchronized
with both itself and `xB`. Claim 2 and sum closure applied to (6) give
`F_1 ~_p B`. Claim 1 implies `C ~_p B`; after common convolution with `x`,
we have `xC ~_p xB`. Sum closure then gives `F_1 ~_p xB`.

## Ordinary synchronization of `B` and `C`

Write `b_k=[x^k]B`, `c_k=[x^k]C`, and use the convention that binomial
coefficients outside their support vanish. Ordinary synchronization amounts to

$$
U_k:=c_kb_k-c_{k-1}b_{k+1}\geq0,
\qquad
V_k:=c_kb_k-c_{k+1}b_{k-1}\geq0. \tag{7}
$$

For the first two indices, direct simplification gives

$$
\begin{aligned}
U_1&=\frac{(s+1)(s+4)}2,\\
V_1&=\frac{s^2+7s+4}2,\\
U_2&=\frac{(s+1)\bigl(s(s-1)(s-2)+12pq\bigr)}{12},\\
V_2&=\frac{s(s-1)(s^2+3s+8)}{12},
\end{aligned} \tag{8}
$$

all nonnegative for integer `s >= 0`.

For `k >= 3`, put `S_j=binom(s,j)` and
`e_j=binom(p,j)+binom(q,j)`. Then

$$
U_k=S_k^2-S_{k-1}S_{k+1}
+\sum_{r\in\{p,q\}}
\left(S_k\binom r{k-1}-S_{k-1}\binom rk\right). \tag{9}
$$

Every term is nonnegative: the first is binomial log-concavity, while

$$
\frac{\binom rk}{\binom r{k-1}}
=\frac{r-k+1}{k}
\leq\frac{s-k+1}{k}
=\frac{S_k}{S_{k-1}}.
$$

For `V_3`, direct symmetric simplification gives

$$
V_3=
\frac{s(s-1)(s-2)\bigl(s^3+6s^2+5s-24pq\bigr)}{144}\geq0, \tag{10}
$$

because `pq <= s^2/4` makes the last factor at least `s^3+5s`.

It remains to prove `V_k >= 0` for `k >= 4`. We use the product-binomial
basis. Vandermonde's identity and the product identity for two binomial
coefficients give

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

Expanding `V_k` with (11)--(12), only indices with
`k <= r:=i+j <= 2k` occur. Set `h=2k-r`. The coefficient of
`binom(p,i)binom(q,j)` is

$$
\begin{aligned}
\mathcal V_{i,j}^{(k)}={}&
\frac{\binom rk\binom kh}{k-h+1}\\
&+\frac{\binom{k-1}{h-1}}{k-1}
\left[
\binom i{k-2}(2-j)+\binom j{k-2}(2-i)
\right]. \tag{13}
\end{aligned}
$$

The negative part of the bracket in (13) is at most

$$
(j-2)_+\binom i{k-2}+(i-2)_+\binom j{k-2}
\leq2\binom rk. \tag{14}
$$

Indeed, each summand is bounded by the corresponding Vandermonde term
`binom(j,2)binom(i,k-2)` or its transpose, and each such term is at most
`binom(r,k)`. Therefore

$$
\mathcal V_{i,j}^{(k)}\geq
\binom rk\left(
\frac{\binom kh}{k-h+1}
-\frac{2\binom{k-1}{h-1}}{k-1}
\right)\geq0. \tag{15}
$$

The last inequality is equivalent to
`k(k-1) >= 2h(k-h+1)`. It is checked directly for `k=4`; for `k >= 5` it
follows from `h(k-h+1) <= (k+1)^2/4`. Thus every coefficient in the
nonnegative basis is nonnegative, proving `V_k >= 0` and Claim 1.

## Partial synchronization of `B` and `xC`

First prove the diagonal, or weak-synchronization, inequalities

$$
W_k:=2c_{k-1}b_k-c_{k-2}b_{k+1}-c_kb_{k-1}\geq0. \tag{16}
$$

The exceptional low indices simplify to

$$
\begin{aligned}
W_1&=s+3,\\
W_2&=\frac{2s^3+9s^2+13s+6pq}{6},\\
W_3&=\frac{s^5+s^4+s^3-s^2-2s-12pq(s^2-s+2)}{24}.
\end{aligned} \tag{17}
$$

The first two are immediate. For `s >= 2`, using `pq <= s^2/4` bounds the
numerator of `W_3` below by

$$
s(s-2)(s^3+4s+1)\geq0;
$$

the cases `s=0,1` give equality directly.

For `k >= 4`, expand (16) in the same product-binomial basis. Only
`k-1 <= r:=i+j <= 2k-1` occurs; put `h=2k-1-r`. The coefficient is

$$
\begin{aligned}
\mathcal W_{i,j}^{(k)}={}&
\frac{2\binom r{k-1}\binom{k-1}h}{k-h+1}\\
&+\frac{\binom{k-1}{h-1}}{(k-1)(k-h+1)}
\Bigl[
\binom i{k-2}\bigl(2(i-k+2)-j(j-1)\bigr)\\
&\hspace{43mm}+
\binom j{k-2}\bigl(2(j-k+2)-i(i-1)\bigr)
\Bigr]. \tag{18}
\end{aligned}
$$

Whenever `binom(i,k-2)` is nonzero, `i-k+2 >= 0`, and similarly for `j`.
The magnitude of the remaining negative part is at most

$$
j(j-1)\binom i{k-2}+i(i-1)\binom j{k-2}
\leq4\binom rk, \tag{19}
$$

again by taking the two relevant Vandermonde terms. Comparing (19) with the
first line of (18) reduces nonnegativity to

$$
k(k-1)\geq2h,
$$

which holds for `k >= 4` and `0 <= h <= k`. Hence every coefficient in (18)
is nonnegative, so (16) holds at every index.

It remains to upgrade the diagonal inequalities to full partial
synchronization. Let
`alpha_k=c_k/c_{k-1}` and `beta_k=b_k/b_{k-1}` on the positive supports.
Claim 1 gives the interlacing inequalities

$$
\beta_{k+1}\leq\alpha_k,
\qquad
\alpha_{k+1}\leq\beta_k. \tag{20}
$$

For `m >= n+1`, log-concavity and (20) give

$$
c_{m-1}b_n\geq c_mb_{n-1},
\qquad
c_{n-1}b_m\geq c_{n-2}b_{m+1}. \tag{21}
$$

Adding (21) is exactly (3) for `xC` and `B`. When `m=n`, it is (16), and
outside the positive supports the inequality is automatic. Thus
`xC ~_p B`, proving Claim 2 and therefore the base relations (4).

Combining the base relations with the induction proves the theorem.

## Exact replay and proof boundary

The proof above is algebraic and all-parameter. The verifier independently
checks:

- formula (1) against direct tree DP;
- recurrence (2);
- the two base partial-synchronization relations on a finite box;
- the product-binomial coefficient identities and nonnegative coefficient
  bounds through a configurable symbolic index range; and
- the exact failure of ordinary synchronization inside the recurrence at
  `p=q=2,t=3`.

Replay:

```text
python3 verify_connector_partial_sync_route_20260808.py
```

Certificate:

```text
results/connector_partial_sync_route_20260808.json
```

The finite replay audits the formulas; it is not the reason the theorem holds.
The universal proof is the coefficientwise argument above.

## Prior-art boundary

Galvin--Hilyard already prove all unsubdivided double stars log-concave as the
`n=1` case of their stronger alternating-caterpillar theorem. Their theorem
does not, as stated, settle the non-alternating caterpillars
`(p,0,...,0,q)` treated here. No novelty claim is made without an authenticated
prior-art search.

References:

- D. Galvin and J. Hilyard, *The independent set sequence of some families of
  trees*, Australasian Journal of Combinatorics 70 (2018), 236--252.
- H. Hu, D. G. L. Wang, F. Zhao, and T. Zhao, *Convolution preserves partial
  synchronicity of log-concave sequences*, Mathematical Inequalities &
  Applications 20 (2017), 91--103.
