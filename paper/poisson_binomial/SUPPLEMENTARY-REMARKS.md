# Supplementary remark: Pitman's ratio bound and the adjacent-ratio drop (1.6)

This remark was in Section 5.3 of earlier drafts of *Variance and local log-concavity of Poisson–binomial laws*. It was moved here to save space. No proof in the paper depends on it. Equation numbers refer to the paper; the notation for exponentially tilted laws is that of the proof in Section 5.2.

For every Bernoulli sum with $V\geq1$, the bound (1.6) gives
$$
\frac{f_{D+1}}{f_D}\leq1-\frac1{4V}.
$$

Pitman's ratio bound [Pitman 1997, equation (21)], combined with the tilted variance bound $V(t)\leq e^tV$ for $t\geq0$, gives
$$
\frac{f_{k+1}}{f_k}<\frac{V}{V+k-\mathbb EW}\qquad\text{for integers } \mathbb EW<k<n.
$$
Indeed, $f_{k+1}/f_k<e^{-t_k}$, where $\mu(t_k)=k$, and
$$
k-\mathbb EW=\int_0^{t_k}V(t)\,dt\leq V\bigl(e^{t_k}-1\bigr).
$$

At $k=D$ this bound tends to one when $D-\mathbb EW\downarrow0$. An example is $\operatorname{Bin}(n,p)$ with $np=\ell+1-\varepsilon$ for an integer $\ell\geq1$ and $p<\varepsilon$, where $D=\ell+1$. So Pitman's ratio bound does not yield the drop (1.6).
