# Authenticated MathSciNet pass, 2026-10-04
<!-- SUMMARY: MathSciNet keyword searches plus forward-citation lists of the four papers a prior result would most likely cite; nothing states or implies Theorem 1.1 · status: bounded negative, venue-record condition met · updated: 2026-10-04 -->

**Access.** Brett's U of T library session in Chrome, through the proxy
`mathscinet-ams-org.myaccess.library.utoronto.ca`. Brett said "I'm signed in
in Chrome. Drive it." The session passed through without a password. Claude
Code (Opus 5.5) drove the browser and read the results pages. The search was
kept to normal manual scale: 15 searches or record views.

**Interface notes, so the counts can be read correctly.**
- The default search is literal and covers all fields. Quoted phrases
  combined with wildcards, and wildcards on hyphenated words, return
  nothing.
- A results page read before it has finished loading reports "0 results".
  Every count below was taken after a wait, except the one flagged.

## Keyword searches (all fields)

| Results | Query | Relevant |
|---:|---|---|
| 17 | `Poisson binomial log-concave` | none. See MR3873905 below |
| 19 | `Bernoulli log-concavity` | none. Many hits match the journal name *Bernoulli* |
| 1 | `Turán inequality Poisson binomial` | none (a 2002 Festschrift) |
| 18 | `Newton inequalities real zeros coefficients` | none. Katkova–Shapiro–Vishnyakova 2024 (MR4725153) asks which Newton-type coefficient inequalities hold for real-rooted polynomials. It gives no variance-scaled deficit bound |
| 23 | `Poisson binomial mode` | none |
| (0) | `"Poisson binomial" AND log-concav*`, `Bernoulli log-concav*` | syntax failures, superseded by the plain-word searches above |
| (0) | `Bernoulli sums log-concave` | read before the page loaded, so unreliable. Covered by `Bernoulli log-concavity` |

## Forward citations (`cp:` lists)

These are the papers a prior result of this kind would most likely cite.

| Citing papers | Cited paper | Relevant |
|---:|---|---|
| 6 | Hillion–Johnson 2016, MR3456338 (the cubic inequalities the proof uses) | none. Entropy papers, plus Johnson 2017, which is cited |
| 17 | Johnson 2017, MR3729642 (c-log-concavity, discrete log-Sobolev) | none. The nearest are Aravinda–Marsiglietti–Melbourne 2022 (concentration for ULC laws) and Aryan–López-Rivera–Shenfeld 2025, ECP (Wu's log-Sobolev via Poisson–Föllmer) |
| 63 | Pitman 1997, MR1429082 (real-rooted coefficients as Bernoulli sums) | none. See MR4700316 below. Marsiglietti–Melbourne 2026 (MR5033003) is already cited |
| 42 | Darroch 1964, MR0164359 (the Poisson–binomial mode) | none. Baillon–Cominetti–Vaisman 2016 is already cited. Singhal 2022 bounds maximum point probabilities |

## Records read

- **MR3873905**, Hillion–Johnson–Saumard, *Statist. Probab. Lett.* 145 (2019).
  - The review states that strong log-concavity with parameter σ implies $\mathrm{Var}\le\sigma$, with discrete analogues characterizing the Poisson and binomial laws.
  - That is the curvature-to-variance direction, the same as Johnson's Lemma 5.3, $c\le 1/\mathbb EW$. Theorem 1.1 goes from variance to curvature at the single point $D$, for Bernoulli sums. So it does not imply the theorem.
  - It is a natural companion citation if Brett wants one.
- **MR4700316**, Marsiglietti–Melbourne, *Discrete Comput. Geom.* 71 (2024). Dilation, Prékopa–Leindler, deviation and convolution inequalities for log-concave sequences. No deficit bound.

## Conclusion (bounded)

Taken together with the zbMATH Open and arXiv sweep in `REPORT.md`, nothing
found states or implies a lower bound of order $1/V$ on the normalized Turán
deficit of a Poisson–binomial law at its first descent or its mode. The paper's
"we are not aware of an earlier proof" stands on a recorded search.

The venue record's condition, an authenticated MathSciNet pass plus a delta
search since July, is met. The delta search is the arXiv filter in `REPORT.md`;
MathSciNet's newest relevant items are listed above. This is still not
exhaustive: free-text matching is literal, and an expert could know of a result
stated in other terms.

Possible citations surfaced, all optional:
- Guo 2021 (adjacent Turán differences);
- Hillion–Johnson–Saumard 2019 (the converse direction);
- Chen–Wang 2026 (quantitative ULC).
