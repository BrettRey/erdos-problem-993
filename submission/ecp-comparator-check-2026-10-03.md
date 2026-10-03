# ECP comparator check
<!-- SUMMARY: PB paper compared with ECP's template and three recent ECP papers · status: no required changes; abstract and keyword length recommended · updated: 2026-10-03 -->

**Question (Brett, 2026-10-03):** compare the paper with something similar in
the same journal to see whether anything needs to change.

## Sources actually read

| Source | What it gave | Access |
|---|---|---|
| ECP/EJP class sample, `ejpecp` v1.11.1, `texmf-dist/doc/latex/ejpecp/sample.tex` (TeX Live 2026) | Official front matter, source rules, back-matter order | local |
| Pinelis, "Preservation of the Bernstein property for sums of independent random variables", ECP 30 (2025) 1–8, doi:10.1214/25-ECP712 | Metadata only, from the Project Euclid landing page in Brett's browser: two-sentence abstract, 3 keywords, primary MSC 60E15 only, 8 pp. This is the venue record's closest format and subject comparator | Full PDF not downloaded (a download needs Brett's go-ahead) |
| Niles-Weed, "Reweighted information inequalities", ECP 2026, pp. 1–13 (arXiv:2603.13135v1) | Abstract, structure, acknowledgements | arXiv version, author formatting |
| Chewi, Pooladian, Zhang, "Stability of the Kim–Milman flow map", ECP 2026, pp. 1–13 (arXiv:2511.01154v2) | Abstract, structure, acknowledgements | arXiv version, author formatting |

How I found them: OpenAlex query on ECP (source S87494164), 2023–2026, with
the terms "log-concave", "Bernoulli", "Poisson binomial", "unimodal" and "sums
of independent". The two 2026 papers are the closest in topic (inequalities for
log-concave or near-log-concave laws) and in length (13 pp.) that have open
arXiv copies. Project Euclid and the Michigan Tech repository both returned
bot challenges to automated requests; I didn't try to get around them.

## Comparison

| Feature | Ours | Comparators | Verdict |
|---|---|---|---|
| Length | 13 pp. | 8, 13, 13 pp. | Fine (the IMS limit recorded in the venue record is ≤ 13) |
| Abstract | 210 words, 8 sentences | about 25 (Pinelis), 29, 89 words | **Long.** Recommend roughly 110–130 words |
| Keywords | 9 | 3 (Pinelis); 4 in the class sample | **Long.** Recommend about 6 |
| MSC | primary 60E15; secondary 60C05; 05A20 | Pinelis: primary 60E15 only | Fine; secondary codes are allowed by the class |
| Main result placement | Theorem 1.1 on p. 2 | main theorem in §1 in both 2026 papers | Fine |
| §1 size | pp. 1–5 of 13; includes Props 1.2, 1.4, 1.5 and Cor 1.3 with proofs, Table 1, Example 1.6 | the 2026 papers state results in §1 and prove them later | Optional: move the proofs of the secondary results to a short later section |
| AI acknowledgement | 270-word acknowledgement crediting each system specifically; title footnote on p. 1 | Niles-Weed (ECP 2026) credits "the technical assistance of ChatGPT-5.2-Pro and Claude Opus 4.6" in one clause | ECP has published an AI credit. Ours is much longer because the contributions were substantive. Keep the specifics (Brett's credit rule); tightening the wording is optional |
| Funding | none stated | both 2026 papers state funding in the acknowledgements | Confirm "no funding" (portal open item) |
| Back-matter order | supplement, acknowledgements, bibliography | class sample: supplement, bibliography, acknowledgements | Optional: match the sample. The class prints `acks` wherever it is placed |
| Bibliography | `amsplain` with MR numbers (all 11), `.bbl` inlined by the source-bundle builder | sample requires inline `thebibliography` and MR links via `amsplain` | Compliant |
| Macros, labels | `\Pp`, `\Var`, `\Ber`; prefixed labels | sample: no one-letter macros, prefixed labels | Compliant |

## Recommendation

Nothing is required. Two changes would bring the paper into line with ECP
practice: a shorter abstract and fewer keywords. The rest is optional.
