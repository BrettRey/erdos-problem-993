# Level/category and terminology pass
Target: `paper/poisson_binomial/main.tex`
## Scope
This pass will enforce type discipline at the manuscript's mathematical interfaces: random variables, laws, probability mass functions, indices, scalar masses, derived statistics, auxiliary bounds, proof objects, programs, and data files. “No synonymy and no homonymy” will govern the paper's own technical vocabulary. Ordinary English polysemy and explicitly attributed source terminology will remain where eliminating them would make the prose unnatural or misreport a source.

No proof strategy, numerical constant, or cited result will change. One current literature-scope sentence does require correction: its universal novelty claim omits the theorem's hypothesis (V\ge 1) and is false without it.
## Typed ontology
| Level | Objects | Canonical language |
| --- | --- | --- |
| Stochastic construction | (B_i,X,Y) and their constituents | Bernoulli **summand**; success probability (p_i); finite Bernoulli sum |
| Distribution | law of (W) | Poisson–binomial **law**; binomial law; extended Bernoulli-sum law only when discussing Gnedin |
| Discrete representation | (f,h,g) | probability mass function (**pmf**) |
| Support and position | support, (D,c,d,u) | support index; first-descent index; rightmost modal index |
| Scalar mass | (f_k), (p=f_c) | mass; maximal mass |
| Derived quantity | (V,q_k,\delta_k) | variance; adjacent mass ratio; normalized Turán deficit |
| Auxiliary proof object | (L_r,R_r,b_r,c_r,w_i) | normalized-mass lower bound; auxiliary weight |
| Exact proof device | Bernstein expansion | positive rational Bernstein coefficients; compact-range Bernstein verification |
| Computational object | Python source / JSON / hashes | generator or checker **program**; certificate **data file**; digest |
| Conditional formal proof | Lean development | conditional formalization under the stated recurrence assumption |
## Canonical terminology decisions
| Preferred term | Terms to retire in this sense | Boundary that must remain visible |
| --- | --- | --- |
| Bernoulli summand | factor, parameter (when (B_i) is meant) | A polynomial factor is not a random summand and cannot be probabilistically independent. |
| pmf | mass function (after the initial expansion), distribution (when the function is meant) | A random variable has a law; a law has a pmf; (f_k) is one mass. |
| rightmost modal index | largest mode, rightmost mode | The index is (c=D-1); the maximal mass is (p=f_c). Gnedin's attributed “leading mode” remains source-specific. |
| first-descent index | first strict descent, postmodal index | (D) is the destination index of the first strict decrease and is not the modal index. |
| normalized Turán deficit | curvature, three-mass deficit |     |
| (\delta_D=1-f_{D-1}f_{D+1}/f_D^2) is distinct from the raw adjacent-ratio drop (1-f_{D+1}/f_D). |     |     |
| adjacent mass ratio | coefficient ratio | Dümgen–Wellner's index-weighted ratio remains separately named and displayed. |
| normalized-mass lower bound | expression, mass lower bound, weight (for (L_r,R_r)) | (L_r,R_r) bound (f_{c\pm r}/p), not the unnormalized masses. Once placed in (A), (w_i) are auxiliary weights. |
| variance-scaled bound | variance-normalized / variance-scale bound | The theorem multiplies the deficit by (V); it does not divide by (V). |
| positive Bernstein expansion / verification | Bernstein certificate (for the mathematical argument) | “Certificate” will be reserved for the released JSON data files. |
## High-impact prose changes
1. Retitle the paper as **“A variance-scaled bound for the Turán deficit at the first descent of Poisson–binomial laws”** and use **“Turán deficit at first descent of Poisson–binomial laws”** as the short title. This removes the curvature/deficit synonym pair and corrects “variance-normalized.”
  
2. Define “probability mass function (pmf)” and thereafter use `pmf` for (f), `mass` for (f_k), and `law` for the distribution.
  
3. State explicitly that (c:=D-1) is the rightmost modal index and that (p=f_c) is the maximal mass. Rewrite the support sentence so indices lie in the support and their corresponding masses are positive.
  
4. Replace every stochastic use of “factor” by “Bernoulli summand.” In the signed corollary, say that all Bernoulli summands occurring in (X) and (Y) are mutually independent.
  
5. Use “first-descent index” after the definition. Attribute symmetry and strict log-concavity to the pmf, not vaguely to the law.
  
6. Use “normalized Turán deficit” throughout. Preserve “raw adjacent-ratio drop” as a separate quantity and only claim the proved inequality between it and (\delta_D).
  
7. Use “normalized-mass lower bound” for (L_r,R_r). Replace “their first moment vanishes” with the literal identity (\sum r w_r=0), since the auxiliary weights are not a probability law.
  
8. Separate programs, certificate data, digests, the positive-Bernstein argument, and the conditional Lean formalization. Say exactly which finite cells the programs check and repeat that none gives end-to-end machine verification of the theorem.
  
9. Correct the novelty sentence to restrict the claim to finite Poisson–binomial laws with (V\ge1).
  
10. Remove the accidentally duplicated sentence “Write (Y=\sum_jY_j).”
  
## Symbol hygiene
Repeated conventional symbols are harmless only when their scopes are disjoint. Two local collisions are worth removing because they cross representation levels in the same proof:

- Rename the generic continuous smoothing density (h) to (\varphi), keeping (h_k) solely for the earlier discrete pmf.
  
- Rename the center (a) in the smoothing argument to (x_0), keeping (a=(1-2\delta)/(1-\delta)) solely for the propagation parameter.
  

I will not globally rename locally scoped counters such as (m), nor the conventional (p_i) and maximal-mass abbreviation (p), because a mechanical rewrite would increase proof risk without improving a reader's local type judgments.
## Preservation gates
- (D) must remain the first-descent index and (c=D-1) the rightmost modal index.
  
- (\delta_D) must remain the normalized Turán deficit; the raw drop remains only a corollary.
  
- (L_r,R_r) must remain lower bounds for normalized masses (f_{c\pm r}/p).
  
- The novelty claim must include (V\ge1).
  
- The certificate/checker language must not imply independent human audit or end-to-end theorem verification.
  
- Gnedin's “leading mode,” max-atom, ULC, and (c)-log-concavity may remain as marked source terms rather than being silently assimilated to the paper's vocabulary.
  
## Verification after editing
1. Run the terminology search and house-style audit.
  
2. Run a fresh read-only proofreading and level/category audit.
  
3. Build with XeLaTeX and biber; inspect warnings and visually inspect affected pages.
  
4. Run the unit tests and the independent exact-certificate checker to confirm that prose cleanup did not disturb formulas or supporting artifacts.
