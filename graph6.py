"""Parse nauty's graph6 format into adjacency lists."""


def parse_graph6(line: bytes) -> tuple[int, list[list[int]]]:
    """Parse a single graph6-encoded line into (n, adjacency_list).

    Parameters
    ----------
    line : bytes
        A graph6-encoded graph (trailing newline stripped automatically).

    Returns
    -------
    (n, adj) where adj[v] is a sorted list of neighbours of vertex v.
    """
    s = line.strip()
    if s.startswith(b">>graph6<<"):
        s = s[10:]

    pos = 0

    # Decode n
    if s[pos] < 126:
        n = s[pos] - 63
        pos += 1
    elif s[pos] == 126 and s[pos + 1] < 126:
        n = 0
        for i in range(1, 4):
            n = (n << 6) | (s[pos + i] - 63)
        pos += 4
    else:
        n = 0
        for i in range(2, 8):
            n = (n << 6) | (s[pos + i] - 63)
        pos += 8

    # Decode edge bits
    bits = []
    for i in range(pos, len(s)):
        val = s[i] - 63
        for shift in (5, 4, 3, 2, 1, 0):
            bits.append((val >> shift) & 1)

    adj: list[list[int]] = [[] for _ in range(n)]
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                adj[i].append(j)
                adj[j].append(i)
            idx += 1

    return n, adj
