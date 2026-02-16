**To:** dgalvin1@nd.edu
**From:** Brett Reynolds
**Subject:** A subdivision-contraction identity for tree independence polynomials

Dear Professor Galvin,

I've been working computationally on Erdős Problem #993 (unimodality of tree independence sequences) and have a preprint with some results that connect to your work on log-concavity failures in subdivided stars. I'd be grateful for any feedback.

The main result is a subdivision-contraction identity: for any tree T and edge e,

  I(T_e; x) = I(T; x) + x · I(T/e; x)

where T_e is the subdivided tree and T/e is the edge contraction. This decomposes the effect of subdivision into the original polynomial plus a shifted copy of the contracted tree's polynomial.

The identity reduces the question of whether subdivision preserves unimodality to a conjectured mode stability property for edge contraction: |mode(I(T)) - mode(I(T/e))| ≤ 1. I've verified this for all 24.7 million edges in trees on up to 20 vertices. If it holds, any minimal counterexample to the Alavi–Malde–Schwenk–Erdős conjecture is homeomorphically irreducible.

Separately, I prove a chain of structural lemmas (Hub Exclusion, Transfer) that reduce the mode bound to a single conjecture about trees where every vertex has at most one leaf-child. The paper also includes exhaustive verification through n = 26 (447 million trees, independently reproducing your and Kadrawi and Levit's two LC failures) and an asymptotic analysis of the near-miss ratio.

I should be upfront: I'm a linguist, not a mathematician. I came to this problem through an interest in combinatorial structure and did the work with substantial help from LLMs (acknowledged in the paper). I'm confident in the computational results and the identity proof, but I'd value a mathematician's eye on the overall framing and the conjectures.

The preprint is attached. Code and data are at https://github.com/BrettRey/erdos-problem-993.

Thank you for your time.

Best regards,
Brett Reynolds
