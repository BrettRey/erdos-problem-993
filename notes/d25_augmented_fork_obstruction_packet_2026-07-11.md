# D25 GSB-augmented fork: exact all-parameter obstruction

## Verdict

The proposed local family

```text
Q_c <= B_c+N[lambda a_c+(1-lambda)S_c],
0 <= lambda <= 1,                                      (AF_lambda)
```

is false in the required prefix.  One exact tree, rank, and center has a
negative gap already at `lambda=0`, and the gap strictly decreases with
`lambda`.  Thus the witness refutes every fixed or adaptively chosen
`lambda in [0,1]` for this family.

The witness still satisfies the true GSB inequality and has a unimodal
independence sequence.  This is an obstruction to one local allocation, not a
counterexample to Erdős #993.

The exact replay is
`scratch_d25_augmented_fork_obstruction_certificate_20260711.py`; its compact
output is saved in
`results/d25_augmented_fork_obstruction_20260711.json`.

## Candidate and its exact global role

Fix a tree `T`, a rank `r`, and a center `c`.  Put

```text
N   = i_r(T),
a_v = i_r(T-N[v]),
J_c = sum_{u<v in N(c)} i_r(T-(N[u] union N[v])),

Q_c = 2[N J_c-sum_{u<v in N(c)}a_u a_v],
B_c = sum_{u in N(c)}a_u^2/deg(u)
      +a_c sum_{u in N(c)}a_u,
S_c = sum_{u in N(c)}a_u/deg(u).
```

Every distance-two pair has one center, so

```text
sum_c Q_c = 2q2.
```

The two allocations in the right side satisfy

```text
sum_c B_c
  = sum_v a_v^2+2 sum_{uv in E(T)}a_u a_v,

sum_c a_c = sum_c S_c = sum_v a_v.
```

Consequently, summing `(AF_lambda)` would give

```text
2q2 <= sum_v a_v^2+2 sum_{uv in E(T)}a_u a_v
       +N sum_v a_v.                                  (1)
```

This is only the distance-two portion of GSB.  The full coefficient inequality
also contains `2qfar` and would still require the independent aggregate claim
`qfar<=0`.  The augmented-fork route therefore had two separate dependencies
even before the obstruction below.

## Exact witness

Let `G(m,t)=T_(m,t,1)` be Galvin's rooted tree: the root has `m` children,
and every child has `t` pendant paths of length two.  Take five disjoint copies
of `G(60,18)` and join all five roots to one new center `c`.

Each branch has

```text
1+60+2*60*18 = 2221 vertices,
```

so the full tree has

```text
n=11106, edges=11105.
```

Every branch root has full-tree degree `61`.  For one rooted branch, exact
rooted-tree recurrence gives

```text
R=(1+2x)^1080,
E=((1+2x)^18+x(1+x)^18)^60,
P=E+xR.
```

Here `P` is the total branch polynomial, `E` is the root-excluded state, and
`xR` is the root-selected state.  Hence the full independence polynomial is

```text
I_T(x)=P(x)^5+xE(x)^5.
```

Its independence number and prefix boundary are

```text
alpha=5701,
L=ceil((2alpha-1)/3)=3801,
r=L-2=3799.
```

The certificate constructs the displayed polynomials in `ZZ[x]` and
independently compares `P,E,R` with the repository's rooted tree DP on the
explicit 2,221-vertex branch.

## The all-parameter failure

At rank `r=3799`, write

```text
N   = [x^r](P^5+xE^5),
a_c = [x^r]E^5,
a   = [x^r]RP^4,
j   = [x^r]R^2P^3.
```

All five branch roots have the same one-mark count `a`, and every pair has the
same joint count `j`.  Therefore

```text
Q_c = 20(Nj-a^2),
B_c = 5a^2/61+5a_c a,
S_c = 5a/61.
```

Let `D(lambda)` be right side minus left side in `(AF_lambda)`.  Clearing the
positive denominator `61` gives the affine function

```text
61D(lambda)=G0+lambda H,

G0 = 5a^2+305a_c a+5Na-61Q_c,
H  = N(61a_c-5a).
```

Exact integer arithmetic gives

```text
G0 < 0,       H < 0.
```

The equality point is

```text
lambda = -0.02838508823786863213585119225875164172526...
```

and hence lies strictly below the allowed interval.  The compact exact
fingerprints are

```text
G0:  5137 digits, sha256 prefix e6af643fbcbeeb578be9, sign -1,
H:   5138 digits, sha256 prefix cf42ac994945a2c15812, sign -1.
```

Thus `D(lambda)<0` for every `lambda in [0,1]`.

## Global safeguards

The same exact replay verifies all of the following:

- the branch and full constructions are connected trees with the stated
  orders and edge counts;
- `r=3799` is exactly the last required prefix rank;
- the true GSB gap at this rank is positive (5,139 digits, SHA-256 prefix
  `6fc3ed6d0ceaa05e69be`); and
- the complete independence sequence of the witness is unimodal.

The local obstruction therefore cannot be mistaken for a nonunimodal tree or
for a failure of prefix GSB.

## Search context and route decision

The original exact scan covered 184,914 applicable center/rank rows through
every nonisomorphic tree of order 15 and left the whole interval `[0,1]`
feasible.  The first structured hard scan covered 1,936 configurations and
still retained `lambda=1/2`.  A broader scan allowing Galvin states in every
branch role exposed the symmetric five-branch witness above.  This is a useful
warning against treating a large finite interval of feasible parameters as a
structural invariant.

D25 is now frozen.  Replacing the constant by a center-dependent value in
`[0,1]` does not help, since this center fails for the entire interval.  A
future distance-two allocation would need a genuinely different, nonlocal
budget.  Aggregate `qfar<=0` remains a separate theorem-strength proposition;
the present obstruction neither proves nor refutes it.  No broad D28 search is
warranted from this packet.

## Reproduction

Run

```bash
python3 scratch_d25_augmented_fork_obstruction_certificate_20260711.py
```

The run uses exact integer polynomial arithmetic and ends with
`"certificate": "passed"`.
