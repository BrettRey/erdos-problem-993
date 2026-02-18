# Repository Guidelines

## Project Structure & Module Organization
- `search.py`: exhaustive search over non-isomorphic trees; saves counterexamples to `results/`.
- `indpoly.py`: independence polynomial DP and unimodality checks.
- `trees.py`: tree enumeration (nauty `geng` when available, otherwise `networkx`).
- `graph6.py`: graph6 parsing utilities.
- `test_all.py`: `unittest`-based tests (unit + small integration checks).
- `targeted.py`: structured family search (subdivided stars, caterpillars, spiders, brooms).
- `nm_optimizer.py`: evolutionary optimizer for near-miss ratio (multi-arm star search).
- `subdivision_correct.py`: corrected subdivision analysis (A = P_uP_v + xR_uR_v).
- `verify_subdivision_formula.py`: verifies the subdivision-contraction identity.
- `test_forced_bound.py`: Case B bound verification (8.7M trees n<=22).
- `test_spider_extremality.py`: spider maximizes mu among d_leaf<=1 (n<=20).
- `investigate_dleaf1_mean.py`: mean mu vs n/3 for d_leaf<=1 trees.
- `paper/main_v2.tex`: current manuscript (11pp, XeLaTeX + biber).
- `paper/main.tex`: previous manuscript version.
- `notes/`: detailed analysis notes (subdivision_new_findings.md, one_private_status.md, conjecture_A_analysis.md are definitive).
- `README.md`, `STATUS.md`, `CLAUDE.md`: project overview, current status, and agent constraints.

## Build, Test, and Development Commands
- `pip install networkx`: required for the `networkx` tree backend.
- `pip install numpy`: optional speedup for polynomial convolution.
- `brew install nauty`: optional; enables `geng` for fast enumeration.
- `python3 search.py --min-n 1 --max-n 20 --backend auto`: exhaustive check over a range.
- `python3 test_all.py`: run the full test suite with `unittest`.
- Paper build: `cd paper && xelatex main_v2.tex && biber main_v2 && xelatex main_v2.tex && xelatex main_v2.tex`

## Coding Style & Naming Conventions
- Python, 4-space indentation, PEP 8 style.
- Use `snake_case` for functions/variables and lowercase module names.
- Prefer type hints on public functions and keep adjacency lists as `list[list[int]]`.
- No formatter/linter is configured; keep diffs small and readable.

## Testing Guidelines
- Framework: `unittest` in `test_all.py`.
- Run specific tests with: `python3 -m unittest test_all.TestIntegration.test_all_n8_unimodal`.
- Integration checks enumerate trees up to n=10; expect longer runtimes than unit tests.

## Commit & Pull Request Guidelines
- Git history is short; messages are concise and imperative (e.g., "Add MIT license"). Follow that style.
- PRs should include: a summary, commands run, and any performance impact.
- If a counterexample or new search result is added, include the JSON certificate and the exact command/parameters used.
- Avoid committing large generated outputs in `out_erdos993/` or `results/` unless they are essential artifacts.

## Agent-Specific Instructions
- Read `CLAUDE.md` before making claims about results; do not state verification without running computations.
- The current paper is `paper/main_v2.tex`, not `paper/main.tex`.
- See `notes/subdivision_new_findings.md` for the definitive subdivision analysis.
- See `notes/one_private_status.md` for the definitive PNP framework.
- See MEMORY.md for dead ends (do NOT revisit).
