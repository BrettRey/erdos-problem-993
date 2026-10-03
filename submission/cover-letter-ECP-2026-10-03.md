Dear Editor,

I submit "Variance and local log-concavity of Poisson–binomial laws" for consideration in *Electronic Communications in Probability*.

Let $W$ be a sum of independent Bernoulli variables with pmf $f$ and variance $V\geq1$, and let $D$ be the first index at which $f$ decreases. The paper asks how tightly the variance alone constrains the local log-concavity of $f$, measured by the normalized Turán deficit $\delta_k=1-f_{k-1}f_{k+1}/f_k^2$. It proves $V\delta_D\geq1/4$, with no constant above $1/3$ possible, and propagates this to $\delta_k\geq1/(4V+|k-D|)$ at every support index. At the mode it shows $1/(4V+1)\leq\delta_c<2/V$, so $1/V$ is the right scale. The lower bounds on $\delta_k$ that I know of scale with the number of summands or with the success probabilities, and they degenerate at fixed variance. An explicit family shows that ultra-log-concavity alone gives no variance-scaled bound. For eigenvalue counts of Haar unitary matrices in a half circle, the new bound is of order $1/\log N$, where degree-based bounds give $1/N$.

The proof reduces the theorem to a one-variable inequality. On one range, a monotonicity argument reduces it to 32 exact rational inequalities, checked by a program that uses only the Python standard library, with the margins printed in the paper. On the rest, it is proved by Bernstein expansions whose coefficients are printed. The whole proof of the main theorem, including the cited inequalities of Hillion and Johnson and of Bobkov, Marsiglietti, and Melbourne, is also formally verified in the Lean proof assistant, starting from the definition of a Poisson–binomial law. The supplementary archive contains the programs, the Lean projects, an independent second certificate for the first range, and replay instructions. It is deposited at [ZENODO DOI].

Generative AI systems assisted with this work. The title-page footnote and the acknowledgements say which systems were used and for what. I am responsible for the mathematics and for every released file.

The manuscript is not under consideration elsewhere. Its source is visible in a public GitHub repository, but it has not been posted to a preprint server. I have no competing interests to declare.

Sincerely,

Brett Reynolds
Humber College, Toronto
brett.reynolds@humber.ca
ORCID 0000-0003-2407-9448
