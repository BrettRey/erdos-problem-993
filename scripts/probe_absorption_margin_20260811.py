#!/usr/bin/env python3
"""Absorption-margin probe on the Bautista-Ramos pattern-graph families.

verify_partial_sync_obstruction_20260811.py (repo root) builds tree
families F = A + B with positive two-term decompositions of their
independence polynomials and locates where log-concavity (LC) breaks:

  Turan(F)_j = Turan(A)_j + Turan(B)_j + S_j
  Turan(P)_j = p_j^2 - p_{j-1} p_{j+1}
  S_j        = 2 a_j b_j - a_{j-1} b_{j+1} - a_{j+1} b_{j-1}

(out-of-range coefficients read as 0). An LC break at index j means
Turan(F)_j < 0, i.e. the diagonal deficit -S_j exceeds the summand
surplus Turan(A)_j + Turan(B)_j. In the checked corpus, breaks are
confined to reflected depth d = alpha - j <= 3 even though S_j < 0
(the sync-violation precondition) extends to relative depth d/alpha
~ 0.54.

Decision question: at fixed relative depth d/alpha, does the
absorption factor

    rho_j = (Turan(A)_j + Turan(B)_j) / (-S_j)      [defined where S_j < 0]

shrink toward 1 as the family parameter grows (live route to a
counterexample) or stay bounded away from 1 / grow (kill-test
evidence)?

Exact integer/Fraction arithmetic throughout. No floats in any
inequality test; floats appear only as reporting decorations attached
to exact Fraction values.
"""

import json
import os
import sys
import time
from fractions import Fraction

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import verify_partial_sync_obstruction_20260811 as vps  # noqa: E402

pmul = vps.pmul
padd = vps.padd
pshift = vps.pshift
ppow = vps.ppow
RTree = vps.RTree
spider2 = vps.spider2
m_family_instance = vps.m_family_instance
lc_breaks = vps.lc_breaks

P1 = [1, 1]

RESULTS_JSON = os.path.join(REPO_ROOT, "results",
                             "partial_sync_obstruction_20260811.json")
OUT_JSON = os.path.join(REPO_ROOT, "results",
                         "absorption_margin_probe_20260811.json")

BANDS = [
    ("0.05-0.15", Fraction(5, 100), Fraction(15, 100)),
    ("0.15-0.25", Fraction(15, 100), Fraction(25, 100)),
    ("0.25-0.35", Fraction(25, 100), Fraction(35, 100)),
    ("0.35-0.45", Fraction(35, 100), Fraction(45, 100)),
]

M_FAMILIES = [
    ("T12_S24_2", dict(ell=2, nlegs=4, k=2), range(2, 61)),
    ("T13_S24_2", dict(ell=3, nlegs=4, k=2), range(2, 31)),
    ("T17_S25_2", dict(ell=7, nlegs=5, k=2), range(2, 25)),
]

TG_M_RANGE = range(2, 7)     # m = 2..6
TG_T_RANGE = range(4, 13)    # t = 4..12

TG_ANCHORS = {
    "TG_4_6": {"m": 4, "t": 6, "expected": [84, 82, 80, 78]},
    "TG_5_6": {"m": 5, "t": 6, "expected": [105, 103, 101, 99, 97]},
}


# ---------- exact per-index decomposition analysis ----------

def coeff(a, i):
    return a[i] if 0 <= i < len(a) else 0


def turan(P, j):
    return coeff(P, j) * coeff(P, j) - coeff(P, j - 1) * coeff(P, j + 1)


def analyze_decomposition(F, A, B):
    """Per index j with S_j < 0: TA_j, TB_j, S_j, margin_j, rho_j (Fraction).

    Also returns the exact set {j : margin_j < 0} (identity-checked
    against Turan(F)_j directly) and flags if any summand Turan value
    is itself negative (would mean A or B is not log-concave there;
    not expected for these families but checked, not assumed).
    """
    alpha = len(F) - 1
    recs = []
    margin_neg = []
    any_summand_neg = False
    for j in range(len(F)):
        a_j, a_jm1, a_jp1 = coeff(A, j), coeff(A, j - 1), coeff(A, j + 1)
        b_j, b_jm1, b_jp1 = coeff(B, j), coeff(B, j - 1), coeff(B, j + 1)
        TA = a_j * a_j - a_jm1 * a_jp1
        TB = b_j * b_j - b_jm1 * b_jp1
        S = 2 * a_j * b_j - a_jm1 * b_jp1 - a_jp1 * b_jm1
        margin = turan(F, j)
        assert margin == TA + TB + S, (
            f"Turan identity mismatch at j={j}: "
            f"{margin} != {TA} + {TB} + {S}")
        if TA < 0 or TB < 0:
            any_summand_neg = True
        if margin < 0:
            margin_neg.append(j)
        if S < 0:
            d = alpha - j
            rho = Fraction(TA + TB, -S)
            recs.append({
                "j": j, "d": d, "TA": TA, "TB": TB, "S": S,
                "margin": margin, "rho": rho,
            })
    return alpha, recs, margin_neg, any_summand_neg


def band_of(d, alpha):
    rel = Fraction(d, alpha)
    for name, lo, hi in BANDS:
        if lo <= rel <= hi:
            return name
    return None


def band_summary(recs, alpha, absorbed_only=False):
    """min rho per band; None if the band has no qualifying point.

    absorbed_only=True restricts to margin_j >= 0 (S_j < 0 but the
    deficit was absorbed, no actual LC break) -- this is the subset
    the absorption-margin question is really about. It matters for
    families with small alpha, where a relative-depth band can
    overlap the genuine near-top break region (margin_j < 0, rho < 1
    there by construction, since margin = TA+TB+S and rho < 1 iff
    TA+TB < -S iff margin < 0).
    """
    out = {name: None for name, _, _ in BANDS}
    for r in recs:
        if absorbed_only and r["margin"] < 0:
            continue
        name = band_of(r["d"], alpha)
        if name is None:
            continue
        cur = out[name]
        if cur is None or r["rho"] < cur["rho"]:
            out[name] = r
    return out


def frac_str(fr):
    return f"{fr.numerator}/{fr.denominator}"


def band_summary_json(bs):
    out = {}
    for name, r in bs.items():
        if r is None:
            out[name] = None
        else:
            out[name] = {
                "j": r["j"], "d": r["d"],
                "TA": r["TA"], "TB": r["TB"], "S": r["S"],
                "rho_frac": frac_str(r["rho"]),
                "rho_float": float(r["rho"]),
            }
    return out


# ---------- TG_{m,t} family (Bautista-Ramos arXiv:2511.00334) ----------

def T_mt(m, t):
    """Paper's T_{m,t}: root v with m children, each child has t pendant
    P2 legs (Definition 1(iii))."""
    T = RTree()
    v = T.new_vertex()
    T.root = v
    for _ in range(m):
        c = T.new_vertex()
        T.edge(v, c)
        for _ in range(t):
            a = T.new_vertex()
            b = T.new_vertex()
            T.edge(c, a)
            T.edge(a, b)
    return T


def TG_mt(m, t):
    """Paper's TG_{m,t}: m disjoint copies of T_{3,t}, roots joined to a
    new root v0, plus one pendant leaf child of v0 (Definition 1(iv))."""
    T = RTree()
    v0 = T.new_vertex()
    T.root = v0
    for _ in range(m):
        sub = T_mt(3, t)
        base = T.graft(sub)
        T.edge(v0, base + sub.root)
    leaf = T.new_vertex()
    T.edge(v0, leaf)
    return T


def tg_decomposition(m, t):
    """Eq. (4): I(TG_m,t) = I(P1) I(T3,t)^m + x I(S2,t)^{3m}."""
    tree = TG_mt(m, t)
    F = tree.ipoly()
    T3t = T_mt(3, t).ipoly()
    S2t = spider2(t).ipoly()
    A = pmul(P1, ppow(T3t, m))
    B = pshift(ppow(S2t, 3 * m), 1)
    assert padd(A, B) == F, f"Eq.(4) mismatch at (m={m}, t={t})"
    return tree.n, F, A, B


def check_tg_anchors():
    out = {}
    ok = True
    for name, spec in TG_ANCHORS.items():
        m, t, expected = spec["m"], spec["t"], spec["expected"]
        order, F, A, B = tg_decomposition(m, t)
        alpha = len(F) - 1
        computed = sorted(lc_breaks(F), reverse=True)
        match = computed == sorted(expected, reverse=True)
        ok = ok and match
        out[name] = {
            "m": m, "t": t, "alpha": alpha, "order": order,
            "expected": sorted(expected, reverse=True),
            "computed": computed,
            "match": match,
        }
    return ok, out


# ---------- m-direction families ----------

def run_m_family(family, params, m_range, expected_by_m):
    instances = []
    for m in m_range:
        order, F, A, B = m_family_instance(
            params["ell"], params["nlegs"], params["k"], m)
        alpha, recs, margin_neg, any_summand_neg = analyze_decomposition(
            F, A, B)
        expected = expected_by_m.get(m)
        self_check_ok = (expected is not None
                          and sorted(margin_neg) == sorted(expected))
        bs = band_summary(recs, alpha)
        bs_abs = band_summary(recs, alpha, absorbed_only=True)
        deepest = max(recs, key=lambda r: r["d"]) if recs else None
        instances.append({
            "m": m, "order": order, "alpha": alpha,
            "n_Sneg": len(recs),
            "lc_breaks_expected": sorted(expected) if expected is not None
                                   else None,
            "lc_breaks_computed": sorted(margin_neg),
            "self_check_ok": self_check_ok,
            "any_summand_turan_neg": any_summand_neg,
            "deepest_Sneg": (None if deepest is None else {
                "j": deepest["j"], "d": deepest["d"],
                "rel_depth_float": deepest["d"] / alpha,
                "TA": deepest["TA"], "TB": deepest["TB"], "S": deepest["S"],
                "rho_frac": frac_str(deepest["rho"]),
                "rho_float": float(deepest["rho"]),
            }),
            "bands": band_summary_json(bs),
            "bands_absorbed_only": band_summary_json(bs_abs),
        })
    return instances


def load_expected_breaks():
    with open(RESULTS_JSON) as fh:
        d = json.load(fh)
    expected = {}
    for r in d["results"]:
        fam = r["family"]
        if fam in ("T12_S24_2", "T13_S24_2", "T17_S25_2"):
            expected.setdefault(fam, {})[r["m"]] = r["lc_breaks"]
    return expected


def trend_report(instances):
    """For each band, list (m, min_rho_float) across instances that have
    a value there, and classify the trend."""
    report = {}
    for name, _, _ in BANDS:
        pts = [(inst["m"], inst["bands"][name]["rho_float"])
               for inst in instances if inst["bands"][name] is not None]
        if len(pts) < 2:
            trend = "insufficient_data"
        else:
            diffs = [pts[i + 1][1] - pts[i][1] for i in range(len(pts) - 1)]
            if all(x >= -1e-12 for x in diffs) and any(x > 1e-9
                                                         for x in diffs):
                trend = "monotone_up"
            elif all(x <= 1e-12 for x in diffs) and any(x < -1e-9
                                                          for x in diffs):
                trend = "monotone_down"
            elif all(abs(x) < 1e-9 for x in diffs):
                trend = "flat"
            else:
                trend = "mixed"
        report[name] = {"points": pts, "trend": trend}
    return report


def run_tg_family():
    instances = []
    for m in TG_M_RANGE:
        for t in TG_T_RANGE:
            order, F, A, B = tg_decomposition(m, t)
            alpha, recs, margin_neg, any_summand_neg = (
                analyze_decomposition(F, A, B))
            bs = band_summary(recs, alpha)
            bs_abs = band_summary(recs, alpha, absorbed_only=True)
            deepest = max(recs, key=lambda r: r["d"]) if recs else None
            instances.append({
                "m": m, "t": t, "order": order, "alpha": alpha,
                "n_Sneg": len(recs),
                "lc_breaks_computed": sorted(margin_neg),
                "any_summand_turan_neg": any_summand_neg,
                "deepest_Sneg": (None if deepest is None else {
                    "j": deepest["j"], "d": deepest["d"],
                    "rel_depth_float": deepest["d"] / alpha,
                    "TA": deepest["TA"], "TB": deepest["TB"],
                    "S": deepest["S"],
                    "rho_frac": frac_str(deepest["rho"]),
                    "rho_float": float(deepest["rho"]),
                }),
                "bands": band_summary_json(bs),
                "bands_absorbed_only": band_summary_json(bs_abs),
            })
    return instances


def main():
    t0 = time.time()
    expected_breaks = load_expected_breaks()

    families_out = {}
    trends_out = {}
    all_self_checks_ok = True

    for family, params, m_range in M_FAMILIES:
        instances = run_m_family(family, params, m_range,
                                  expected_breaks.get(family, {}))
        bad = [inst["m"] for inst in instances if not inst["self_check_ok"]]
        if bad:
            all_self_checks_ok = False
            print(f"SELF-CHECK FAILURE in {family} at m={bad}",
                  file=sys.stderr)
        families_out[family] = {"params": params, "instances": instances}
        trends_out[family] = trend_report(instances)

    if not all_self_checks_ok:
        print("ABORTING: self-check against "
              "results/partial_sync_obstruction_20260811.json failed.",
              file=sys.stderr)
        sys.exit(1)

    tg_anchors_ok, tg_anchor_detail = check_tg_anchors()
    tg_out = {"anchors_ok": tg_anchors_ok, "anchors": tg_anchor_detail}
    if tg_anchors_ok:
        tg_instances = run_tg_family()
        tg_out["instances"] = tg_instances
        # trend report keyed by t, one series per fixed TG-m value.
        # "bands" mixes in actual near-top breaks when alpha is small
        # (rho < 1 there by construction); "bands_absorbed_only"
        # restricts to margin_j >= 0, the genuine absorption cases.
        tg_trends_by_m = {}
        tg_trends_by_m_absorbed = {}
        for m in TG_M_RANGE:
            sub = [inst for inst in tg_instances if inst["m"] == m]
            renamed = [{"m": inst["t"], "bands": inst["bands"]}
                       for inst in sub]
            tg_trends_by_m[str(m)] = trend_report(renamed)
            renamed_abs = [{"m": inst["t"],
                             "bands": inst["bands_absorbed_only"]}
                            for inst in sub]
            tg_trends_by_m_absorbed[str(m)] = trend_report(renamed_abs)
        tg_out["trends_by_m_fixed_t_varies"] = tg_trends_by_m
        tg_out["trends_by_m_fixed_t_varies_absorbed_only"] = (
            tg_trends_by_m_absorbed)
    else:
        tg_out["note"] = ("TG anchors did not match arXiv:2511.00334's "
                           "stated break indices; TG margin analysis "
                           "skipped.")

    out = {
        "date": "2026-08-11",
        "source_scripts": ["verify_partial_sync_obstruction_20260811.py"],
        "source_results": ["results/partial_sync_obstruction_20260811.json"],
        "source_paper_tg": ("notes/literature/arxiv_2511_00334.txt "
                             "(Bautista-Ramos, arXiv:2511.00334)"),
        "question": ("At fixed relative depth d/alpha, does "
                     "rho_j = (Turan(A)_j+Turan(B)_j)/(-S_j) shrink "
                     "toward 1 as the family parameter grows, or stay "
                     "bounded away from 1 / grow?"),
        "self_checks_ok": all_self_checks_ok,
        "families": families_out,
        "trends": trends_out,
        "tg_family": tg_out,
        "runtime_s": round(time.time() - t0, 2),
    }

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"self_checks_ok={all_self_checks_ok} "
          f"tg_anchors_ok={tg_anchors_ok} "
          f"runtime={out['runtime_s']}s")
    for family in families_out:
        print(f"\n== {family} ==")
        for name, _, _ in BANDS:
            tr = trends_out[family][name]
            print(f"  band {name}: trend={tr['trend']} "
                  f"n_points={len(tr['points'])}")
            if tr["points"]:
                head = tr["points"][:3]
                tail = tr["points"][-3:]
                print(f"    head={head} tail={tail}")


if __name__ == "__main__":
    main()
