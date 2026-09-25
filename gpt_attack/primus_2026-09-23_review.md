# Primus handoff: proposed brief and ZIP contents

Draft for Brett's review, 23 September 2026. This prepares a local upload
package; it does not submit a job, spend credits, publish material, or restart
our own research campaign.

## Proposed assignment

Ask Primus for a proof or exact counterexample to the remaining bounded
depth-three inequality for tree independence sequences. Give it the established
graph-level results and counterexamples to failed auxiliary claims, so it can
start at the current boundary rather than repeat old work.

The broader motivation is Erdős Problem #993, but completing this assignment
would establish only one bounded coefficient window, not solve #993.

For a finite tree T, let n be its number of vertices, alpha its independence
number, and s_d the number of independent sets of size alpha-d. The primary
target is

\[
s_3^2\ge s_2s_4
\]

for every tree with 33 <= n <= 38, alpha in {17,18,19}, and
delta = 2 alpha - n <= 5.

Split s_d = e_d + b_d, where e_d counts independent sets contained in at least
one maximum independent set, and b_d counts the remaining sets. Count each set
once, regardless of how many maximum sets contain it.

The recorded results already cover b_1=0 and the low-density condition

\[
3(\alpha-3)b_4\le(\alpha-7)e_4.
\]

The remaining regime has b_1>0, density above that threshold, and adverse
combined correction

\[
D=b_3^2-b_2b_4+2e_3b_3-e_2b_4-e_4b_2<0.
\]

The useful known margin is

\[
e_3^2-e_2e_4\ge
\frac{5\alpha+17}{27(\alpha-3)}e_2e_4.
\]

Primus may prove that this margin absorbs D, or find another route to the
primary target. Failure of that sufficient margin bound is not automatically
a counterexample to the primary target, still less to #993.

## What the ZIP will contain

- A short paste-ready Primus prompt and a self-contained mathematical brief.
- The exact target, definitions, established dependencies, and graded success
  criteria: full proof; exact counterexample; or a genuine partial theorem with
  the remaining gap stated precisely.
- Relevant Lean source with its pinned toolchain and dependency manifest,
  excluding downloaded dependencies and build caches. Historical verification
  and any checks run while packaging will be labelled separately.
- Small exact witness certificates and portable replay code, tested from the
  extracted ZIP. These prevent rediscovery of false componentwise-positivity
  and blocked-log-concavity claims.
- A concise list of failed routes and optional literature leads. Holant is a
  conditional methodological lead, not an applicable theorem; Boros–Moll is a
  proof-design reference, not a missing bound.
- File hashes, provenance, and instructions for returning reproducible results.

## Instructions to Primus

Start with the supplied evidence and frozen target. Do not require every
component or certificate class to have a nonnegative contribution. Do not
assume real-rootedness, log-concavity, mode–mean localization, or an unproved
representation. A larger numerical census is evidence, not a proof.

Prefer one explicit, checkable advance over a broad survey or a speculative
paper. Return readable mathematics, all scripts and certificates used, exact
commands, and an honest evidence grade. If Lean is used, keep the target
unchanged and report its dependencies and axioms; a conditional scalar lemma
must not be presented as the graph theorem.

## Scope and privacy

Exclude confidential referee reports, submission correspondence, credentials,
API/project account metadata, unrelated research, and the repository's full
history. Do not include the sprawling status journal or obsolete prompt archive.

The default is this focused assignment with #993 as context. If you want Primus
to attack the full open problem independently instead, change that scope here
before the ZIP is assembled.

Primus's public announcement describes a Math mode supporting written proofs,
Lean, and computational experiments:
[Primus for Math](https://lab.cloud/news/primus-for-math/).
Its authenticated upload form and file limits could not be inspected here, so
the archive will use ordinary Markdown, Python, JSON, and Lean files rather than
an assumed service-specific schema.
