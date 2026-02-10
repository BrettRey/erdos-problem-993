# Repository Guidelines

## Project Structure & Module Organization
- `search.py`: exhaustive search over non-isomorphic trees; saves counterexamples to `results/`.
- `erdos993_hunt.py`: heuristic/parallel hunt over structured families and forests; writes certificates and progress to `out_erdos993/`.
- `indpoly.py`: independence polynomial DP and unimodality checks.
- `trees.py`: tree enumeration (nauty `geng` when available, otherwise `networkx`).
- `graph6.py`: graph6 parsing utilities.
- `test_all.py`: `unittest`-based tests (unit + small integration checks).
- `README.md`, `STATUS.md`, `CLAUDE.md`: project overview, current status, and agent constraints.

## Build, Test, and Development Commands
- `pip install networkx`: required for the `networkx` tree backend.
- `pip install numpy`: optional speedup for polynomial convolution.
- `brew install nauty`: optional; enables `geng` for fast enumeration.
- `python search.py --min-n 1 --max-n 20 --backend auto`: exhaustive check over a range.
- `python erdos993_hunt.py --workers 8 --time-limit 36000`: long-running heuristic search with parallel workers.
- `python test_all.py`: run the full test suite with `unittest`.

## Coding Style & Naming Conventions
- Python, 4-space indentation, PEP 8 style.
- Use `snake_case` for functions/variables and lowercase module names.
- Prefer type hints on public functions and keep adjacency lists as `list[list[int]]`.
- No formatter/linter is configured; keep diffs small and readable.

## Testing Guidelines
- Framework: `unittest` in `test_all.py`.
- Run specific tests with: `python -m unittest test_all.TestIntegration.test_all_n8_unimodal`.
- Integration checks enumerate trees up to n=10; expect longer runtimes than unit tests.

## Commit & Pull Request Guidelines
- Git history is short; messages are concise and imperative (e.g., “Add MIT license”). Follow that style.
- PRs should include: a summary, commands run, and any performance impact.
- If a counterexample or new search result is added, include the JSON certificate and the exact command/parameters used.
- Avoid committing large generated outputs in `out_erdos993/` or `results/` unless they are essential artifacts.

## Agent-Specific Instructions
- Read `CLAUDE.md` before making claims about results; do not state verification without running computations.
