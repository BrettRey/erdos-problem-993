# `main_v2` Patch Plan (2026-02-21)

## Purpose
Prepare a low-friction patch path for `paper/main_v2.tex` once the remaining symbolic lemmas close, without rewriting the manuscript.

## Anchor Map in `paper/main_v2.tex`
- Intro contribution list:
  - current proof-claim bullet near `paper/main_v2.tex:99`
- Conjecture A block:
  - `\begin{conjecture}[Conjecture A]` near `paper/main_v2.tex:373`
- Mean/tie-fugacity section:
  - `\subsection{The mean bound}` near `paper/main_v2.tex:380`
  - tie-fugacity remark near `paper/main_v2.tex:404`
- Discussion open-problem language:
  - near `paper/main_v2.tex:843` and `paper/main_v2.tex:852`
- Certificates table:
  - Conjecture A and tie-fugacity rows near `paper/main_v2.tex:930` and `paper/main_v2.tex:940`

## Patch Set A (If Full Symbolic Closure Lands)
1. Promote Conjecture A statement:
   - Replace Conjecture at `paper/main_v2.tex:373` with theorem/corollary form.
   - Keep the original conjecture wording in a historical remark, not main statement.

2. Add a short closure theorem block after tie-fugacity setup:
   - New theorem: universal tie-fugacity/bridge inequality chain used in closure.
   - Add concise proof sketch with references to Route 1 + STRONG C2 lemmas.
   - Keep heavy algebra in notes/appendix text; main text gets dependency chain only.

3. Update introduction claims:
   - Keep bullet structure at `paper/main_v2.tex:97`-`paper/main_v2.tex:101`.
   - Replace provisional phrasing with final theorem-level phrasing tied to new labels.

4. Update discussion:
   - Remove “Conjecture A remains open” phrasing near `paper/main_v2.tex:843`.
   - Replace with “ECMS remains open” (or exact remaining open items at closure time).

5. Update certificates table:
   - Move closure rows from “verified computationally” framing to theorem references.
   - Keep computational rows as supporting evidence, not primary proof claims.

## Patch Set B (If Only Partial Symbolic Closure Lands)
1. Keep Conjecture A as conjecture at `paper/main_v2.tex:373`.
2. Add “partial closure” proposition(s) in `\subsection{The mean bound}`:
   - e.g., full-frontier validated inequality with explicit symbolic reduction target.
3. Tighten wording in intro/discussion to avoid overclaim:
   - state “reduced to X symbolic lemma” with explicit equation.
4. Preserve all certificate rows as computational only.

## Minimal Edit Strategy
1. Add new theorem/proposition labels first.
2. Update only local cross-references (no global refactor).
3. Rebuild PDF and check for:
   - broken refs/cites
   - stale “open” language conflicting with final theorem status
4. Keep a single coherent narrative:
   - “proved” only for symbolic results
   - “verified” for computational certificates

## Post-Patch Verification Checklist
1. `rg -n "Conjecture A remains open|remain open|open problem" paper/main_v2.tex`
2. `rg -n "proved|verified|conjecture" paper/main_v2.tex` and manually reconcile phrasing.
3. Build:
   - `cd paper && xelatex main_v2.tex && biber main_v2 && xelatex main_v2.tex && xelatex main_v2.tex`
4. Confirm no theorem/proposition label collisions.
