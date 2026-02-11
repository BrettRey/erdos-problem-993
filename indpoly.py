"""Independence polynomial computation and unimodality check for trees."""

try:
    from numpy import convolve as _np_convolve

    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

# int64 safe threshold: each output coeff is a sum of at most
# min(len(a),len(b)) products, each bounded by max_a * max_b.
# Keep total below 2^62 to avoid silent overflow in numpy int64.
_INT64_SAFE = 2**62


def _polymul_python(a: list[int], b: list[int]) -> list[int]:
    """Pure Python convolution with arbitrary-precision integers."""
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def _polymul(a: list[int], b: list[int]) -> list[int]:
    if not a or not b:
        return []
    if _HAS_NUMPY:
        max_a = max(abs(x) for x in a)
        max_b = max(abs(x) for x in b)
        max_terms = min(len(a), len(b))
        if max_a > 0 and max_b > 0 and max_terms * max_a * max_b < _INT64_SAFE:
            return _np_convolve(a, b).astype(int).tolist()
    return _polymul_python(a, b)


def independence_poly(n: int, adj: list[list[int]]) -> list[int]:
    """Compute the independence polynomial of a tree.

    Uses iterative post-order DP to avoid recursion limits.

    Parameters
    ----------
    n : int
        Number of vertices.
    adj : list[list[int]]
        Adjacency list.

    Returns
    -------
    Coefficient list [i_0, i_1, ..., i_alpha] where i_k = number of
    independent sets of size k.
    """
    if n == 0:
        return [1]
    if n == 1:
        return [1, 1]

    # Root at vertex 0, BFS to get parent/children and post-order
    parent = [-1] * n
    children: list[list[int]] = [[] for _ in range(n)]
    visited = [False] * n
    order = []  # will hold post-order
    visited[0] = True
    bfs_queue = [0]
    head = 0
    while head < len(bfs_queue):
        v = bfs_queue[head]
        head += 1
        for u in adj[v]:
            if not visited[u]:
                visited[u] = True
                parent[u] = v
                children[v].append(u)
                bfs_queue.append(u)

    # Post-order from BFS order (reverse of BFS = approximate post-order,
    # but we need true post-order; use iterative DFS instead)
    stack = [(0, False)]
    order = []
    while stack:
        v, processed = stack.pop()
        if processed:
            order.append(v)
            continue
        stack.append((v, True))
        for c in children[v]:
            stack.append((c, False))

    # DP: dp0[v] = poly when v excluded, dp1[v] = poly when v included
    dp0: list[list[int]] = [[] for _ in range(n)]
    dp1: list[list[int]] = [[] for _ in range(n)]

    for v in order:
        if not children[v]:
            # Leaf
            dp0[v] = [1]
            dp1[v] = [0, 1]
        else:
            # dp0[v] = product over children c of (dp0[c] + dp1[c])
            prod = [1]
            for c in children[v]:
                summand = _polyadd(dp0[c], dp1[c])
                prod = _polymul(prod, summand)
            dp0[v] = prod

            # dp1[v] = x * product over children c of dp0[c]
            prod = [1]
            for c in children[v]:
                prod = _polymul(prod, dp0[c])
            dp1[v] = [0] + prod  # multiply by x

    return _polyadd(dp0[0], dp1[0])


def _polyadd(a: list[int], b: list[int]) -> list[int]:
    """Add two polynomials represented as coefficient lists."""
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


def is_unimodal(seq: list[int]) -> bool:
    """Check whether a sequence is unimodal.

    A sequence is unimodal if it never increases after a decrease:
    no index i with seq[i-1] > seq[i] < seq[i+1].

    More precisely: once the sequence starts decreasing, it never
    increases again.
    """
    if len(seq) <= 2:
        return True
    decreasing = False
    for i in range(1, len(seq)):
        if seq[i] > seq[i - 1]:
            if decreasing:
                return False
        elif seq[i] < seq[i - 1]:
            decreasing = True
    return True


def is_log_concave(seq: list[int]) -> bool:
    """Check whether a sequence is log-concave.

    A sequence is log-concave if a_k^2 >= a_{k-1} * a_{k+1} for all
    valid k. Uses integer arithmetic to avoid floating-point issues.
    """
    for k in range(1, len(seq) - 1):
        if seq[k] * seq[k] < seq[k - 1] * seq[k + 1]:
            return False
    return True


def log_concavity_ratio(seq: list[int]) -> tuple[float, int]:
    """Measure worst log-concavity violation.

    Returns (max_ratio, position) where ratio = a_{k-1}*a_{k+1} / a_k^2.
    Ratio > 1 means log-concavity fails at that position.
    Larger ratio = worse violation.
    """
    worst = 0.0
    worst_k = -1
    for k in range(1, len(seq) - 1):
        sq = seq[k] * seq[k]
        if sq == 0:
            continue
        ratio = (seq[k - 1] * seq[k + 1]) / sq
        if ratio > worst:
            worst = ratio
            worst_k = k
    return worst, worst_k


def near_miss_ratio(seq: list[int]) -> tuple[float, int]:
    """Measure how close a unimodal sequence is to violating unimodality.

    After the first strict decrease, tracks the maximum ratio
    a_{j+1}/a_j in the "should be non-increasing" tail. A ratio > 1
    would be a unimodality violation. The closer to 1, the closer to
    a counterexample.

    Returns (max_ratio, position_j) or (0.0, -1) if the sequence
    never decreases.
    """
    if len(seq) <= 2:
        return 0.0, -1

    # Find first strict decrease
    first_descent = -1
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            first_descent = i
            break

    if first_descent == -1:
        return 0.0, -1  # monotone non-decreasing

    # In the descending tail, find max a_{j+1}/a_j
    worst = 0.0
    worst_j = -1
    for j in range(first_descent, len(seq) - 1):
        if seq[j] == 0:
            continue
        ratio = seq[j + 1] / seq[j]
        if ratio > worst:
            worst = ratio
            worst_j = j
    return worst, worst_j
