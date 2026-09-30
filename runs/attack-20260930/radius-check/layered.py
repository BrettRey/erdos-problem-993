"""Exact closed forms (python-flint) for layered trees: root at level 0; every
level-j vertex has b[j] children; level L = len(b) vertices are leaves.
Returns I(T) and, for each level i, J_i = I(T - N[v]) for a level-i vertex v,
plus level sizes and degrees."""
from flint import fmpz_poly
X = fmpz_poly([0, 1]); ONE = fmpz_poly([1])
def layered(b):
    L = len(b)
    E = [None] * (L + 1); P = [None] * (L + 1)
    E[L], P[L] = ONE, X
    for j in range(L - 1, -1, -1):
        E[j] = (E[j + 1] + P[j + 1]) ** b[j]; P[j] = X * E[j + 1] ** b[j]
    I = E[0] + P[0]
    # A0[j], A1[j]: T minus subtree(x), x at level j >= 1, with parent excluded / included
    A0 = [None] * (L + 1); A1 = [None] * (L + 1)
    for j in range(1, L + 1):
        if j == 1:
            A0[1] = (E[1] + P[1]) ** (b[0] - 1); A1[1] = X * E[1] ** (b[0] - 1)
        else:
            A0[j] = (A0[j - 1] + A1[j - 1]) * (E[j] + P[j]) ** (b[j - 1] - 1)
            A1[j] = A0[j - 1] * X * E[j] ** (b[j - 1] - 1)
    J = []
    for i in range(L + 1):
        part = ONE
        if i + 1 <= L:                                   # children removed; grandchildren subtrees remain
            gc = (E[i + 2] + P[i + 2]) ** b[i + 1] if i + 2 <= L else ONE
            part *= gc ** b[i]
        if i >= 1:
            part *= (E[i] + P[i]) ** (b[i - 1] - 1)      # siblings' full subtrees
            if i - 1 >= 1: part *= (A0[i - 1] + A1[i - 1])   # everything above the parent
        J.append(part)
    N = [1]
    for j in range(L): N.append(N[-1] * b[j])
    deg = [b[0]] + [b[j] + 1 for j in range(1, L)] + [1]
    return I, J, N, deg
def coeffs(p): return [int(c) for c in p.coeffs()]
