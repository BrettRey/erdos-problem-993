"""
ERDOS PROBLEM #993 -- an open-problem work packet (single self-contained file).

Std-lib only. All arithmetic exact (Python big ints). No floating point is used
to DECIDE unimodality or log-concavity. Run `python3 erdos993_handoff.py` to
self-test: it reproduces the two published n=26 log-concavity failures
coefficient-for-coefficient and the paper's near-miss champion ratios.

--------------------------------------------------------------------------------
THE PROBLEM
  Tree T on n vertices; i_k(T) = number of independent sets of size k (no two
  adjacent; the empty set counts, i_0 = 1). Independence polynomial
  I(T;x) = sum_k i_k x^k, degree alpha.
  CONJECTURE (Alavi, Malde, Schwenk, Erdos 1987; = Erdos #993): the sequence
  (i_0,...,i_alpha) is UNIMODAL (rises to a peak, then falls; ties allowed) for
  every tree. Open since 1987. General graphs are NOT constrained this way
  (their sequences can be any positive sequence), e.g. a 26-vertex graph
  realizes 1,26,15,20,15,6,1 which dips then rises.

TWO WAYS TO WIN -- refutation is an explicit success mode
  PROVE : unimodality for all trees; or a new nontrivial tree class; or a
          reduction that provably shrinks the open case.
  REFUTE: exhibit one tree whose sequence dips then rises again (a VALLEY).
          A single valid edge list ends a 39-year-old conjecture. This is the
          maximal outcome, not a consolation prize.

GRADED SUB-TARGETS (claim the highest you can actually reach)
  T1  Reproduce the state of the art (run the self-test; know where a
      counterexample must live: n >= 30, outside the known families).
  T2  A new sufficient structural condition for unimodality, WITH PROOF,
      covering trees not already handled by the known families.
  T3  Prove or refute a conditional hinge: the tree "mode <= ceil(mean)"
      localization, or Edge Contraction Mode Stability (see KNOWN RESULTS).
  T4  Beat the near-miss frontier: a tree with near-miss ratio above the best
      known for its vertex count, OR a proof a family cannot exceed 1 - c/n
      for an explicit c.
  T5  Settle it: a full proof, or a verified valley.

BINDING VERIFICATION RULES
  1. Exact integer arithmetic only. NEVER use float root-finding or float
     coefficient arithmetic to decide shape. (A claw-free control once produced
     a fabricated result from numpy.roots rounding noise.) Use exact
     Sturm/sympy or certified Arb isolation if you touch polynomial zeros.
  2. Label every claim THEOREM (proved, proof written out) or EVIDENCE
     (computational, with the exact range checked). "Verified for n <= 20" is
     EVIDENCE, not a theorem.
  3. If you assert a bound, give the constant.
  4. A counterexample MUST ship a certificate: the exact edge list AND the exact
     integer sequence, reproducible by check_counterexample(n, edges). Nothing
     else counts.

SELF-GRADE YOUR REPLY (first line)
  Begin with one of: SOLVED / COUNTEREXAMPLE / PARTIAL / NO-PROGRESS, then grade
  yourself honestly. An honest PARTIAL with one genuinely proved lemma beats a
  grand claim that dissolves on checking. If you only reproduced known facts,
  that is NO-PROGRESS (T1) -- say so.

--------------------------------------------------------------------------------
KNOWN RESULTS (build past these; do not re-derive)
  PROVED unimodal: paths & centipedes (AMSE 1987); regular caterpillars
    (Galvin-Hilyard 2018); Fibonacci trees (Bencs 2018); spiders -> log-concave
    (Li et al. 2025), hence brooms; two families containing the 2025
    log-concavity breakers (2026).
  Tail always safe: i_k strictly decreasing for k >= ceil((2*alpha-1)/3)
    (Levit-Mandrescu 2006) -- any valley sits in the first ~2/3.
  Companion paper (doi:10.5281/zenodo.19100781):
    - THEOREM mean(T) < n/3 for dleaf<=1 trees (decimation + Steiner peeling).
    - THEOREM leaf-attachment: attaching s leaves gives near-miss
      nm(s) = 1 - C/s + O(1/s^2), 4 <= C <= 8 -- so extremal families approach
      the wall but the gap shrinks like C/s and never closes.
  CONDITIONAL hinges (proving/refuting either = T3):
    - mode <= ceil(mean) for trees (OPEN; Darroch 1964 proved the
      Poisson-binomial analogue; verified all dleaf<=1 trees n<=23).
    - Edge Contraction Mode Stability + a tail condition (verified n<=19)
      => a minimal counterexample is homeomorphically irreducible (no degree-2
      vertices). Both conjectural.
  SEARCH FRONTIER (EVIDENCE): all 8,691,747,673 trees on n<=29 checked, zero
    counterexamples -> any counterexample has n>=30. Log-concavity already
    fails (2 trees at n=26, 0 at 27, 19 at 28, plus infinite families) but
    unimodality never has. Best near-miss among all trees through n=28 is
    nm=0.8571 (n=27); tuned multi-arm stars M(s;a1..ak) reach 0.9437 (n=75),
    0.9792 (n=200), 0.9959 (n=1000), all obeying the 1-C/s law.
  ORIENTATION: the small regime is exhausted and extremal families are provably
    bounded away from the wall. A counterexample, if any, is LARGE and OUTSIDE
    these families -- or the conjecture is true. Blind random search over
    n < 30 is wasted; that space is closed.

Full paper: Brett Reynolds, "Mean bounds, structural reductions, and exhaustive
verification for tree independence polynomial unimodality" (2026),
doi:10.5281/zenodo.19100781. Repo: https://github.com/BrettRey/erdos-problem-993
--------------------------------------------------------------------------------
"""

DATA = {
 "problem": "Erdos #993: is (i_0,...,i_alpha) unimodal for every tree?",
 "log_concavity_failures_n26": [
  {
   "id": 1,
   "n": 26,
   "edges": [
    [
     0,
     12
    ],
    [
     1,
     13
    ],
    [
     2,
     14
    ],
    [
     3,
     15
    ],
    [
     4,
     16
    ],
    [
     5,
     17
    ],
    [
     6,
     18
    ],
    [
     7,
     19
    ],
    [
     8,
     20
    ],
    [
     9,
     21
    ],
    [
     10,
     22
    ],
    [
     0,
     23
    ],
    [
     1,
     23
    ],
    [
     2,
     23
    ],
    [
     11,
     23
    ],
    [
     3,
     24
    ],
    [
     4,
     24
    ],
    [
     5,
     24
    ],
    [
     6,
     24
    ],
    [
     11,
     24
    ],
    [
     7,
     25
    ],
    [
     8,
     25
    ],
    [
     9,
     25
    ],
    [
     10,
     25
    ],
    [
     11,
     25
    ]
   ],
   "poly": [
    1,
    26,
    300,
    2040,
    9142,
    28551,
    63933,
    103736,
    121376,
    100144,
    55499,
    18683,
    2979,
    51,
    1
   ],
   "lc_fail_index": 13,
   "lc_ratio": 1.1453287197231834,
   "graph6": "Y???????????_?O?C??_?A??C??C??A???_??C???O?[?_?F`???^???"
  },
  {
   "id": 2,
   "n": 26,
   "edges": [
    [
     0,
     12
    ],
    [
     1,
     13
    ],
    [
     2,
     14
    ],
    [
     3,
     15
    ],
    [
     4,
     16
    ],
    [
     5,
     17
    ],
    [
     6,
     18
    ],
    [
     7,
     19
    ],
    [
     8,
     20
    ],
    [
     9,
     21
    ],
    [
     0,
     22
    ],
    [
     10,
     22
    ],
    [
     1,
     23
    ],
    [
     2,
     23
    ],
    [
     10,
     23
    ],
    [
     11,
     23
    ],
    [
     3,
     24
    ],
    [
     4,
     24
    ],
    [
     5,
     24
    ],
    [
     11,
     24
    ],
    [
     6,
     25
    ],
    [
     7,
     25
    ],
    [
     8,
     25
    ],
    [
     9,
     25
    ],
    [
     11,
     25
    ]
   ],
   "poly": [
    1,
    26,
    300,
    2037,
    9089,
    28147,
    62183,
    98968,
    112870,
    90178,
    48086,
    15498,
    2372,
    48,
    1
   ],
   "lc_fail_index": 13,
   "lc_ratio": 1.0295138888888888,
   "graph6": "Y???????????_?O?C??_?A??C??C??A???_??C?C?O?K@_?F@???|???"
  }
 ],
 "near_miss_champions": [
  {
   "n": 75,
   "s": 66,
   "arms": [
    6,
    2
   ],
   "near_miss": 0.9437,
   "source": "paper Table (multiarm)"
  },
  {
   "n": 100,
   "s": 88,
   "arms": [
    6,
    3,
    2
   ],
   "near_miss": 0.9575,
   "source": "paper Table (multiarm)"
  },
  {
   "n": 200,
   "s": 183,
   "arms": [
    5,
    5,
    4,
    2
   ],
   "near_miss": 0.9792,
   "source": "paper Table (multiarm)"
  },
  {
   "n": 500,
   "s": 483,
   "arms": [
    5,
    5,
    4,
    2
   ],
   "near_miss": 0.9918,
   "source": "paper Table (multiarm)"
  },
  {
   "n": 1000,
   "s": 983,
   "arms": [
    5,
    5,
    4,
    2
   ],
   "near_miss": 0.9959,
   "source": "paper Table (multiarm)"
  }
 ],
 "frontier_facts": {
  "exhaustive_no_counterexample_through_n": 29,
  "total_trees_checked_through_n29": 8691747673,
  "best_near_miss_all_trees_through_n28": {
   "nm": 0.8571425,
   "n": 27,
   "pos_k": 13
  },
  "log_concavity_failure_counts": {
   "26": 2,
   "27": 0,
   "28": 19
  },
  "tail_strictly_decreasing_from_k": "ceil((2*alpha-1)/3)  (Levit-Mandrescu 2006)",
  "counterexample_lower_bound_vertices": 30
 }
}

from fractions import Fraction


def build_adj(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def indpoly(n, edges):
    """Independence polynomial coefficients [i_0,...,i_alpha] as exact ints."""
    adj = build_adj(n, edges)
    parent = [-1] * n
    order = []
    seen = [False] * n
    stack = [0]
    seen[0] = True
    while stack:
        u = stack.pop()
        order.append(u)
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True
                parent[w] = u
                stack.append(w)
    if len(order) != n:
        raise ValueError("not a connected tree on n vertices")

    def mul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    out[i + j] += x * y
        return out

    def add(a, b):
        out = [0] * max(len(a), len(b))
        for i, x in enumerate(a):
            out[i] += x
        for i, x in enumerate(b):
            out[i] += x
        return out

    ex = [None] * n
    inc = [None] * n
    for u in reversed(order):
        e, i = [1], [0, 1]
        for w in adj[u]:
            if parent[w] == u:
                e = mul(e, add(ex[w], inc[w]))
                i = mul(i, ex[w])
        ex[u], inc[u] = e, i
    return add(ex[0], inc[0])


def analyze(poly):
    """All comparisons exact. Returns a dict of shape facts."""
    m = len(poly) - 1
    D = None
    for k in range(1, m + 1):
        if poly[k] < poly[k - 1]:
            D = k
            break
    descended = False
    valley = False
    for k in range(1, m + 1):
        if poly[k] < poly[k - 1]:
            descended = True
        elif descended and poly[k] > poly[k - 1]:
            valley = True
            break
    lc_fail = None
    for k in range(1, m):
        if poly[k] * poly[k] < poly[k - 1] * poly[k + 1]:
            lc_fail = k
            break
    nm = Fraction(0)
    nm_pos = None
    if D is not None:
        for j in range(D, m):
            r = Fraction(poly[j + 1], poly[j])
            if r > nm:
                nm = r
                nm_pos = j
    mode = max(range(m + 1), key=lambda k: poly[k])
    return {
        "alpha": m,
        "mode": mode,
        "first_descent": D,
        "unimodal": not valley,
        "valley": valley,
        "log_concave": lc_fail is None,
        "lc_fail_index": lc_fail,
        "near_miss": float(nm),
        "near_miss_exact": (nm.numerator, nm.denominator),
        "near_miss_pos": nm_pos,
    }


def path(n):
    return n, [[i, i + 1] for i in range(n - 1)]


def star(n):
    return n, [[0, i] for i in range(1, n)]


def multiarm(s, arms):
    """M(s; a1,...,ak): hub 0 with s pendant leaves and paths of given lengths.
    n = 1 + s + sum(arms). The empirically extremal near-miss family."""
    edges, nid = [], 1
    for _ in range(s):
        edges.append([0, nid]); nid += 1
    for a in arms:
        prev = 0
        for _ in range(a):
            edges.append([prev, nid]); prev = nid; nid += 1
    return nid, edges


def broom(handle, bristles):
    return multiarm(bristles, [handle - 1])


def spider(legs, leg_len):
    return multiarm(0, [leg_len] * legs)


def caterpillar(spine, leaves_per_spine):
    edges, nid = [], spine
    for i in range(spine - 1):
        edges.append([i, i + 1])
    for i in range(spine):
        for _ in range(leaves_per_spine):
            edges.append([i, nid]); nid += 1
    return nid, edges


def random_tree(n, seed):
    """Deterministic random labelled tree (attach each new vertex to a uniform earlier one)."""
    state = seed & 0x7FFFFFFF

    def rnd(bound):
        nonlocal state
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        return state % bound

    return n, [[rnd(i), i] for i in range(1, n)]


def graph6_decode(s):
    d = [ord(c) - 63 for c in s]
    n = d[0]
    bits = []
    for b in d[1:]:
        bits += [(b >> i) & 1 for i in range(5, -1, -1)]
    edges, idx = [], 0
    for j in range(1, n):
        for i in range(j):
            if bits[idx]:
                edges.append([i, j])
            idx += 1
    return n, edges


def check_counterexample(n, edges):
    """Gate for a claimed valley. Returns (is_valley, poly, info)."""
    poly = indpoly(n, edges)
    return analyze(poly)["valley"], poly, analyze(poly)


if __name__ == "__main__":
    ok = True
    print("[1] Two published n=26 log-concavity failures (Kadrawi-Levit 2025)")
    for rec in DATA["log_concavity_failures_n26"]:
        poly = indpoly(rec["n"], rec["edges"])
        match = poly == rec["poly"]
        info = analyze(poly)
        print("    tree %d: poly reproduced = %s, log_concave = %s (fails at k=%s), unimodal = %s"
              % (rec["id"], match, info["log_concave"], info["lc_fail_index"], info["unimodal"]))
        ok &= match and (not info["log_concave"]) and info["unimodal"]

    print("[2] Near-miss champions (multi-arm stars; source: paper Table)")
    for c in DATA["near_miss_champions"]:
        n, edges = multiarm(c["s"], c["arms"])
        info = analyze(indpoly(n, edges))
        good = n == c["n"] and abs(info["near_miss"] - c["near_miss"]) < 5e-4
        print("    M(%d; %s) n=%d (expect %d): nm=%.4f (expect %s), gap G=%.2f  %s"
              % (c["s"], c["arms"], n, c["n"], info["near_miss"], c["near_miss"],
                 n * (1 - info["near_miss"]), "OK" if good else "MISMATCH"))
        ok &= good

    print("[3] Sanity: a path is unimodal & log-concave; a valley is detected")
    pinfo = analyze(indpoly(*path(40)))
    fake = analyze([1, 26, 15, 20, 15, 6, 1])
    print("    path(40): unimodal=%s, log_concave=%s" % (pinfo["unimodal"], pinfo["log_concave"]))
    print("    split-graph seq 1,26,15,20,15,6,1: valley=%s (must be True)" % fake["valley"])
    ok &= pinfo["unimodal"] and pinfo["log_concave"] and fake["valley"]

    print("\nSELF-TEST:", "PASS" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)
