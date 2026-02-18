"""Tree enumeration: auto-selects geng (nauty) or networkx fallback."""

import shutil
import subprocess
from collections.abc import Iterator

from graph6 import parse_graph6


def trees(n: int, backend: str = "auto") -> Iterator[tuple[int, list[list[int]]]]:
    """Yield (n, adj) for each non-isomorphic tree on n vertices.

    Parameters
    ----------
    n : int
        Number of vertices.
    backend : str
        "geng", "networkx", or "auto" (tries geng first).
    """
    if n < 1:
        return
    if n == 1:
        yield (1, [[]])
        return

    if backend == "auto":
        backend = "geng" if shutil.which("geng") else "networkx"

    if backend == "geng":
        yield from trees_geng(n)
    elif backend == "networkx":
        yield from trees_networkx(n)
    else:
        raise ValueError(f"Unknown backend: {backend!r}")


def trees_geng(
    n: int, res: int | None = None, mod: int | None = None
) -> Iterator[tuple[int, list[list[int]]]]:
    """Enumerate trees using nauty's geng.

    Parameters
    ----------
    n : int
        Number of vertices.
    res, mod : int or None
        If both given, enumerate only partition res/mod (0-indexed).
        This enables parallelism: workers get res=0..mod-1.
    """
    edges = n - 1
    cmd = ["geng", "-cq", str(n), f"{edges}:{edges}"]
    if res is not None and mod is not None:
        cmd.extend([f"{res}/{mod}"])

    with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            yield parse_graph6(line)


def trees_geng_raw(
    n: int, res: int | None = None, mod: int | None = None
) -> Iterator[tuple[int, list[list[int]], bytes]]:
    """Like trees_geng but also yields the raw graph6 bytes."""
    edges = n - 1
    cmd = ["geng", "-cq", str(n), f"{edges}:{edges}"]
    if res is not None and mod is not None:
        cmd.extend([f"{res}/{mod}"])

    with subprocess.Popen(cmd, stdout=subprocess.PIPE) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            n_out, adj = parse_graph6(line)
            yield (n_out, adj, line.strip())


def trees_networkx(n: int) -> Iterator[tuple[int, list[list[int]]]]:
    """Enumerate trees using networkx (slower, no external deps)."""
    import networkx as nx

    for T in nx.nonisomorphic_trees(n):
        adj: list[list[int]] = [[] for _ in range(n)]
        for u, v in T.edges():
            adj[u].append(v)
            adj[v].append(u)
        # Sort for deterministic order
        for v in range(n):
            adj[v].sort()
        yield (n, adj)
