# Erdős 993 as a statement about the mean free-vertex count; Mason's inequality below the Levit–Mandrescu threshold
<!-- SUMMARY: Exact identities rewrite unimodality as "the mean number of free vertices over uniform k-independent sets never re-crosses the line k+1"; the sufficient condition "that mean is decreasing" is Mason's inequality (b) for the tree's independence complex, verified at every index for all 23.9M trees n<=23, failing first in exactly one tree at n=24 (T_{3,3,4}, at k=alpha-2), and holding below the threshold in every known LC-failing family; the universal slope-1 version is refuted by T_{80,12} · status: bounded evidence and a conjecture, no theorem · updated: 2026-09-02 -->

Date: 2026-09-02. Status: exact identities (verified), one refuted
strengthening, one surviving conjecture with bounded evidence. Nothing here is a
theorem about all trees, and no manuscript claim changes.

Scripts: `scripts/free_count_stats_20260902.py` (trivariate DP for the free-count
moments, asserts the identities), `scripts/mason_slope_20260902.py` and
`scripts/mason_slope_threshold_20260902.py` (polynomial-level diagnostics),
`scripts/mason_slope_search_20260902.py` (adversarial hill-climb),
`scripts/mason_b_census_gentreeg_20260902.py` (exhaustive check on `gentreeg -p`
output). Numbers: `results/mason_free_count_evidence_20260902.json`.

## 1. The identities

Fix a tree $T$ and an independent set $S$ with $|S|=k$. Call a vertex *free*
if it lies outside $N[S]$, write $f(S)$ for the number of free vertices, $F_S$
for the forest they induce, $e(S)$ for its edge count and $c(S)=f-e$ for its
component count. Every $(k+1)$-set arises from exactly $k+1$ of its $k$-subsets
by adding a free vertex, so

\[
(k+1)\,i_{k+1}=\sum_{|S|=k} f(S),
\qquad
(k+1)(k+2)\,i_{k+2}=\sum_{|S|=k}\bigl(f^2-f-2e\bigr)(S). \tag{1}
\]

The second identity counts ordered pairs of distinct non-adjacent free vertices.
Both were asserted by the DP on all 32,507 trees with $n\le16$.

Write $E_k$ for expectation over a uniformly random $k$-independent set and

\[
\mu_k:=\frac{(k+1)\,i_{k+1}}{i_k}=E_k[f].
\]

So $\mu_k$ is the mean number of free vertices, and $i_{k+1}/i_k=\mu_k/(k+1)$.
The sequence descends at $k$ exactly when $\mu_k<k+1$. From (1),

\[
\mu_{k+1}-\mu_k=\frac{\operatorname{Var}_k(f)}{\mu_k}-1-\frac{2E_k[e]}{\mu_k}. \tag{2}
\]

There is also an exact form of the Levit–Mandrescu tail bound. For a forest,
$\alpha(F)=|F|-\nu(F)$, so $f=2\alpha(F_S)-\operatorname{def}(F_S)$ where
$\operatorname{def}$ counts vertices missed by a maximum matching. Setting
$\operatorname{gap}(S)=\alpha(T)-k-\alpha(F_S)\ge0$ (how far $S$ falls short
of extending to a maximum independent set),

\[
\mu_k \;=\; 2(\alpha-k)\;-\;E_k\bigl[\,2\operatorname{gap}+\operatorname{def}\,\bigr]. \tag{3}
\]

Verified by brute force on all 200 trees with $n\le10$. Dropping the
correction gives $\mu_k\le2(\alpha-k)$, which is $\le k+1$ for
$k\ge(2\alpha-1)/3$: that is the Levit–Mandrescu theorem, and (3) says exactly
what its slack is.

## 2. Unimodality as a crossing statement

Unimodality of $(i_k)$ is equivalent to: for every $k$,

\[
\mu_k<k+1\;\Longrightarrow\;\mu_{k+1}<k+2 .
\]

The line $k+1$ rises by one per step, so a sufficient condition is a slope bound:

- **(S1)** $\mu_{k+1}-\mu_k\le1$;
- **(N)** $\mu_{k+1}\le\mu_k$, i.e. $(k+2)\,i_{k+2}\,i_k\le(k+1)\,i_{k+1}^2$,
  i.e. $i_{k+1}^2\ge\bigl(1+\tfrac1{k+1}\bigr)\,i_k\,i_{k+2}$.

(N) is Mason's second inequality for the independence complex of $T$ (the
form proved for matroids by Huh, Schröter and Wang, arXiv:1806.02675). It is
strictly stronger than log-concavity and strictly weaker than the binomial
ultra-log-concavity that Lorentzian methods deliver. Since $k!\,i_k$ counts
ordered independent $k$-tuples, (N) says that sequence is log-concave. By (2),
(N) $\Leftrightarrow \operatorname{Var}_k(f)\le E_k[f]+2E_k[e]$ and
(S1) $\Leftrightarrow \operatorname{Var}_k(f)\le 2E_k[f]+2E_k[e]$: a
sub-Poisson dispersion bound on the free count, with an allowance for free
edges. In pair form, (N) at level $k$ reads

\[
\sum_{u\ne v,\ u\not\sim v} P_k(u,v\text{ both free})\;\le\;\Bigl(\sum_v P_k(v\text{ free})\Bigr)^2 .
\]

Because rises $i_j<i_{j+1}$ are impossible for $j\ge\theta:=\lceil(2\alpha-1)/3\rceil$
(Levit–Mandrescu), a valley needs a slope larger than one at some
$k\le\theta-2$. So the versions that matter are the **prefix** ones,
$k\le\theta-2$. Prefix-(N) implies the hypothesis of Levit–Kadrawi
Corollary 2.20 (log-concavity for $j\le\theta-1$) and implies prefix-(S1);
either implies Erdős 993.

## 3. The universal versions are false

Universal (N) fails at the top of every log-concavity-failing tree, since (N)
is stronger than log-concavity: the first $n=26$ witness has
$\mu_{13}-\mu_{12}=+0.052$.

Universal (S1) also fails, but only deep in the tail. In the bouquet
$T_{m,t}$ (a root with $m$ children, each carrying $t$ pendant $P_2$ legs),
the slope exceeds one once $m$ is large relative to $2^t/t$: $T_{14,8}$ gives
$1.08$ and $T_{80,12}$ gives $2.40$, the latter at $k=961$ with $\alpha=1040$,
$\theta=693$ and $\mu_k=0.45$. The sequence there is descending by a factor of
about $2000$ per step, so the crossing statement is nowhere near violated. Over
the whole grid $t\le14$, $m\le100$ the first index with slope $>1$ satisfies
$k-\theta\ge29$.

## 4. Evidence for the prefix versions

All arithmetic is exact. "First positive slope" is the least $k-\theta$ with
$\mu_{k+1}>\mu_k$; prefix-(N) needs it to be $\ge-1$.

| Data | Trees | First positive slope ($k-\theta$) | Max slope, $k\le\theta-2$ |
|---|---|---|---|
| All trees $4\le n\le23$ | 23,942,356 | never (N holds at every index) | $-2/n$ (the star) |
| All trees $n=24$ (`gentreeg`, 8 shards) | 39,299,897 | 2, in exactly one tree: $T_{3,3,4}$, $\alpha=13$, $\theta=9$, slope $+0.0046$ at $k=11$ (no LC failure; the first LC failure is $T_{3,4,4}$ at $n=26$) | $-2/24$ |
| Census LC failures $n\le32$ (`results/lc_census_20260814/`) | 1,228 | 3 | $-1.54$ |
| $TG_{m,t}$, $m\le20$, $t\le10$ (up to 19 consecutive breaks, LC ratio to 13.7) | 200 | 4 | $-1.15$ |
| $T_{m,t}$, $t\le14$, $m\le100$ | 288 | 4 | $-1.55$ |
| Bautista–Ramos $m$-families (three families, up to 3 consecutive breaks) | 111 | 12 | $-1.82$ |
| Evolutionary champions, 2026-08-14 break-depth campaign | 30 | 3 | $-0.81$ |

Within the descending prefix ($k\le\theta-2$ and $\mu_k<k+1$, the only
region where a valley can start) the largest slope seen anywhere is $-0.77$
($n\le16$ exhaustive). Two adversarial hill-climbs maximising exactly that
quantity (leaf moves, regrafts, subdivisions, pendant $P_2$ grafts;
$n\le50$; 2.15M and 9.09M exact evaluations) ended at $-0.82$ and $-0.78$,
with zero unimodality alarms. A third hill-climb on the plain prefix slope
(1.09M evaluations) converged to $-0.04$, which is the star $K_{1,49}$: the
star's $-2/n$ is extremal for (N) in the prefix and nothing found beats it.

So in every tree and family examined, $\mu_k$ is strictly decreasing up to
$k=\theta+1$, and the first failures of (N) sit at distance $\ge2$ above the
threshold (the unique $n=24$ witness) and $\ge3$ in every larger family, inside
the distance-$\ge4$ wall recorded for log-concavity failures on 2026-08-14.
The smallest tree violating Mason (b) is therefore the Kadrawi–Levit tree
$T_{3,3,4}$ on 24 vertices, two orders before the first log-concavity
failure and at the same top-of-sequence position $k=\alpha-2$. The star shows prefix-(N) has no slack
in the ascending region; prefix-(S1) has slack about $1.8$ everywhere in the
descending region.

## 5. What a proof would need, and why this is not one

By (2), prefix-(S1) is the variance bound
$\operatorname{Var}_k(f)\le2E_k[f+e(F_S)]$ for $k\le\theta-2$. Writing
$f=\sum_v X_v$ with $X_v$ the free indicator,
$\operatorname{Var}(f)=\sum_v p_v(1-p_v)+\sum_{u\ne v}\operatorname{Cov}(X_u,X_v)$,
and adjacent pairs contribute at most $2E[e]$. So (N) would follow from
"non-adjacent free indicators are negatively correlated in aggregate", up to
the slack $\sum_v p_v^2$. Pairwise negative correlation is false: in the star
at $k=1$ two leaves have $\operatorname{Cov}=(s-3)/(s+1)^2>0$, and the
inequality holds there only because $\sum_v p_v^2\approx s$ absorbs it. I do
not have a route to the aggregate bound; it is the same difficulty as the
prefix log-concavity target, reached from the probabilistic side.

Two smaller facts that fall out: the product of a log-concave sequence with no
internal zeros and a unimodal sequence is unimodal (convolution with a
$\mathrm{TP}_2$ kernel is variation-diminishing), so a forest is unimodal
whenever all but one component are log-concave; and by (3) a descent at $k$
means $E_k[2\operatorname{gap}+\operatorname{def}]>2\alpha-3k-1$, so a valley
would need that mean correction to fall by more than three in one step,
whereas pointwise it falls by at most one (adding a free vertex $v$ changes
it by $d_{F}(v)-1$).

## 6. Literature check (narrow, evidenced)

Local: `grep -il mason notes/literature/*.txt notes/*.md paper/main_v2.tex`
hits only the Rota-disproof paper and `notes/matroid_nonunimodality_2026-08-14.md`,
neither about tree independence complexes. arXiv API, run 2026-09-02:
`all:"Mason" AND all:"independence polynomial"` (0 results);
`all:"Mason's conjecture" AND all:"independence complex"` (0);
`abs:"independent set" AND abs:"Mason" AND abs:"tree"` (2: Huh–Schröter–Wang
1806.02675 and Ardila–Cristancho–Denham 2601.02547, both about matroids);
`abs:"ultra log-concave" AND abs:"independence polynomial"` (0). OpenAlex
`"independence polynomial" "Mason" trees` returned four unrelated titles.
Semantic Scholar returned HTTP 429. This is not a novelty claim; a proper
outward search (hyperresearch light tier) is the cheap next step if the
conjecture is worth stating publicly.

## 7. Verdict

No proof and no counterexample. The reformulation moves the problem onto the
dispersion of the free-vertex count and identifies Mason's inequality below
the Levit–Mandrescu threshold as the clean intermediate statement, with the
star as its extremal case. Every mechanism I can see for a valley runs into
the same entropy-versus-shift trade-off the 2026-07 barrier law
($1-V\sim4/n$) already quantifies.
