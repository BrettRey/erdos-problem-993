# Proof Context: Endpoint-Aware Curvature Propagation

This packet isolates the first Lean milestone for the universal finite
Poisson-binomial effective-drop theorem proved on 2026-07-10.

For a positive Poisson-binomial mass sequence $f_k$, define

\[
\delta_k=1-\frac{f_{k-1}f_{k+1}}{f_k^2},
\qquad
\delta_0=\delta_n=1.
\]

Hillion--Johnson's cubic inequalities imply the normalized neighboring law

\[
\delta_{k\pm1}(1-\delta_k)\le\delta_k.
\]

The Lean file takes this inequality as the hypothesis `hstep`; formalizing the
Hillion--Johnson theorem itself is a later milestone. Iterating the map
$x\mapsto x/(1-x)$ should give

\[
\delta_{D\pm r}\le\frac{d}{1-rd}
\qquad ((r+1)d<1),
\]

where $d=\delta_D$. If an endpoint occurred in this window, its curvature one
would satisfy

\[
1\le\frac{d}{1-rd}<1,
\]

a contradiction. This endpoint exclusion is mathematically essential.

At the modal crossing, the same propagation gives

\[
q_D\ge\frac{1-2d}{1-d}.
\]

Finally, if $D$ is a strict descent and

\[
R_- = \frac{f_D}{f_{D-1}}<1,
\qquad
R_+ = \frac{f_{D+1}}{f_D},
\]

then

\[
1-R_+\ge 1-\frac{R_+}{R_-}
=1-\frac{f_{D-1}f_{D+1}}{f_D^2}.
\]

That last implication converts the effective quarter bound into the raw
quarter reserve used downstream.

The full informal proof and exact symbolic certificates live outside this
minimal upload packet. This job must not claim to formalize the complete
Poisson-binomial theorem.
