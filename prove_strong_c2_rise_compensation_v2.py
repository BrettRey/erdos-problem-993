#!/usr/bin/env python3
"""Focused structural verification for STRONG C2 rise compensation.

Key findings from v1:
1. mode(P) >= m-1 ALWAYS
2. -mismatch = p1*q0 - p0*q1 (cross-product of P,Q)
3. cross term in LC decomposition always >= 0
4. Residual always >= 0
5. lc_surplus >= b0 * delta_b always (0 failures)

This v2 focuses on:
- Proving neg <= rise via the factored identity and mode(P) >= m-1.
- Checking the key structural bound: since mode(P) >= m-1,
  P is still increasing at m-2, so p_{m-1} >= p_{m-2} (P is non-decreasing up to m-1).
  Similarly b_{m-1} >= b_{m-2} (mode(B) >= m-1).
  When mode(P) >= m-1: both P and B are increasing at index m-2 -> m-1.

The identity: rise - neg = p1*(b1-b0) + b1*(q1-q0)
  = p1*delta_b + b1*delta_q

Since delta_b >= 0 and p1 >= 0, the first term is >= 0.
The second term b1*delta_q is >= 0 when q1 >= q0 (Q not dropping at m-1).

Q drops at m-1 in only 2 of 931K trees. So the Q-drop is ultra-rare.

When Q does NOT drop: rise - neg >= 0 trivially (both terms non-negative).
When Q drops: need p1*delta_b >= b1*(-delta_q), i.e., p1*delta_b >= b1*(q0-q1).

Since delta_b = delta_p + delta_q (where delta_p = p1-p0):
  p1*(delta_p + delta_q) >= b1*(q0-q1)
  p1*delta_p + p1*delta_q >= -b1*delta_q
  p1*delta_p >= -(p1+b1)*delta_q = (p1+b1)*(q0-q1)

Since mode(P) >= m-1 implies delta_p >= 0 (P increasing at m-2->m-1):
  p1*delta_p >= 0 >= 0 ... but (p1+b1)*(q0-q1) > 0 in Q-drop.

So this direct approach doesn't close. Need:
  p1*delta_p >= (p1+b1)*(q0-q1) = (p1+b1)|delta_q|

Since p1*delta_p = p1*(p1-p0) = p1^2 - p1*p0
and (p1+b1)|delta_q| = (p1+b1)(q0-q1),

we need p1^2 - p1*p0 >= (p1+b1)(q0-q1).
RHS = p1*q0 - p1*q1 + b1*q0 - b1*q1
LHS = p1^2 - p1*p0

So: p1^2 - p1*p0 >= p1*q0 - p1*q1 + b1*q0 - b1*q1
    p1*(p1 - p0 - q0 + q1) >= b1*(q0-q1)
    p1*(p1 - p0 - q0 + q1) >= b1*(q0-q1)

Note p1 - p0 - q0 + q1 = (p1+q1) - (p0+q0) = b1 - b0 = delta_b.
So: p1*delta_b >= b1*(q0-q1) = b1*|delta_q|

This is exactly the original condition! We've gone in a circle.

BUT HERE'S THE KEY: the problem reduces to p1*delta_b >= b1*|delta_q|
i.e. (p1/b1) >= |delta_q|/delta_b = |delta_q|/(delta_p + delta_q)
    = |delta_q|/(delta_p - |delta_q|)  ... since delta_q < 0

Since delta_b = delta_p + delta_q = delta_p - |delta_q| > 0 (B is increasing),
delta_p > |delta_q|.

So |delta_q|/delta_b = |delta_q|/(delta_p - |delta_q|) = 1/(delta_p/|delta_q| - 1)

And p1/b1 = p1/(p1+q1).

For the bound, we need:
p1/(p1+q1) >= 1/(delta_p/|delta_q| - 1)

Cross-multiplying (both sides positive):
p1*(delta_p/|delta_q| - 1) >= p1 + q1
p1*delta_p/|delta_q| - p1 >= p1 + q1
p1*delta_p/|delta_q| >= 2*p1 + q1
delta_p/|delta_q| >= (2*p1 + q1)/p1 = 2 + q1/p1

So the condition is: (p1-p0)/|q0-q1| >= 2 + q1/p1.

This is saying P grows much faster than Q shrinks, accounting for the Q/P ratio.

Let me check this numerically on the two Q-drop witnesses.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly, is_log_concave


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=23)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--mod", type=int, default=1)
    ap.add_argument("--out", default="results/prove_strong_c2_rise_v2_n23.json")
    args = ap.parse_args()

    stats = {
        "checked": 0,
        "mismatch_neg": 0,
        "combined_neg": 0,
        "rise_fail": 0,
        "q_drop": 0,
        "q_no_drop_mismatch_neg": 0,  # mismatch neg with q1 >= q0 (IMPOSSIBLE per algebra?)

        # When mode(P) >= m-1 and q1 >= q0: mismatch is automatically >= 0.
        # Because: mismatch = p0*b1 - p1*b0 = p0*(p1+q1) - p1*(p0+q0)
        #        = p0*q1 - p1*q0
        # If mode(P) >= m-1: p1 >= p0.
        # If q1 >= q0: then p0*q1 vs p1*q0:
        #   p1 >= p0 and q1 >= q0 does NOT determine the sign.
        #   Example: p0=3,p1=4,q0=2,q1=3 -> mismatch=9-8=1>0.
        #   Example: p0=2,p1=5,q0=3,q1=4 -> mismatch=8-15=-7<0.
        # So mode(P)>=m-1 and q1>=q0 does NOT imply mismatch>=0.

        # Alternative structural approach:
        # Since P = prod I(T_c) and Q = x*prod dp0[c]:
        # The "ratio" approach compares growth rates.
        # P includes contributions from dp1 in each factor; Q does not.
        # So P tends to grow faster than Q (more "mass" at higher k).
        # Specifically: p_k/q_{k+1} = [prod (dp0+dp1)]_k / [prod dp0]_k
        # = product of [(dp0+dp1)/dp0]_k-like ratios.
        # Each factor (dp0+dp1)/dp0 at index k is 1 + dp1_k/dp0_k,
        # and dp1/dp0 is the conditional probability of including the child.
        # As k increases (more vertices in the IS), dp1 contribution grows.

        # KEY APPROACH: LC of P implies p_{m-1}^2 >= p_{m-2}*p_m.
        # Combined with LC of Q: q_{m-1}^2 >= q_{m-2}*q_m (but we're looking at q_{m-2},q_{m-1}).
        # Cauchy-Schwarz on the cross terms.

        # APPROACH A: show combined = lc_P + lc_Q + residual
        # where lc_P = p1^2-pm*p0, lc_Q = q1^2-qm*q0,
        # residual = 2*p1*q1 - pm*q0 - p0*qm - (p1*q0 - p0*q1)
        #          = 2*p1*q1 - pm*q0 - p0*qm - p1*q0 + p0*q1
        #          = p1*(2*q1-q0) + p0*(q1-qm) - pm*q0
        #          = p1*q1 + p1*(q1-q0) + p0*(q1-qm) - pm*q0
        # We proved residual >= 0 in v1 (min = 4).
        # So combined >= lc_P + lc_Q >= 0. DONE!

        # Wait, is that really what v1 showed? Let me re-verify.
        # v1 showed residual >= 4 > 0 for all 77K trees.
        # So combined = lc_P + lc_Q + residual, with all three terms >= 0.
        # That would prove combined >= 0!

        # But wait: residual was defined in v1 as:
        # residual = q1*(p1+p0) + p1*(q1-q0) - p0*qm - pm*q0
        # And total = lc_P + lc_Q + residual = combined
        # And residual >= 4 always (at n<=20).

        # If residual is always >= 0, then combined = lc_P + lc_Q + residual >= 0
        # since lc_P >= 0 (P is LC) and lc_Q >= 0 (Q is LC).
        # THIS WOULD BE THE PROOF!

        # So the question is: CAN WE PROVE residual >= 0 algebraically?
        # residual = p1*(2*q1-q0) + p0*(q1-qm) - pm*q0

        # Let me expand more carefully.
        # residual = p1*q1 + p1*(q1-q0) + p0*q1 - p0*qm - pm*q0

        # Term 1: p1*q1 >= 0 (trivial)
        # Term 2: p1*(q1-q0) = p1*delta_q (could be negative in Q-drop)
        # Term 3: p0*q1 >= 0 (trivial)
        # Term 4: -p0*qm <= 0
        # Term 5: -pm*q0 <= 0

        # Group positive: p1*q1 + p0*q1 = q1*(p1+p0) = q1*(p0+p1)
        # Group negative: p1*delta_q (if Q-drop) - p0*qm - pm*q0
        # Worst case: all negative terms dominate.
        # -p0*qm - pm*q0 is bounded by LC of Q cross with P.

        # Actually, let me think about this differently.
        # residual = 2*p1*q1 - pm*q0 - p0*qm - p1*q0 + p0*q1
        # = (p1*q1 - pm*q0) + (p1*q1 - p0*qm) - p1*q0 + p0*q1
        # = (p1*q1 - pm*q0) + (p1*q1 - p0*qm) - (p1*q0 - p0*q1)

        # The last term -(p1*q0-p0*q1) = mismatch = p0*q1 - p1*q0 (possibly negative).
        # The first two terms: by Cauchy-Schwarz,
        # p1*q1 >= sqrt(pm*p0)*sqrt(qm*q0)  (from LC of P and Q)
        # So p1*q1 - pm*q0 >= sqrt(pm*p0*qm*q0) - pm*q0
        # and p1*q1 - p0*qm >= sqrt(pm*p0*qm*q0) - p0*qm

        # Hmm, let me try a cleaner approach.
        # By AM-GM: pm*q0 + p0*qm >= 2*sqrt(pm*p0*qm*q0) <= 2*p1*q1
        # (The last inequality follows from p1^2>=pm*p0 and q1^2>=qm*q0,
        #  so p1*q1 >= sqrt(pm*p0)*sqrt(qm*q0) = sqrt(pm*p0*qm*q0).)
        # Actually: p1*q1 >= sqrt(pm*p0*qm*q0) by LC.
        # And pm*q0 + p0*qm >= 2*sqrt(pm*p0*qm*q0) by AM-GM.
        # So 2*p1*q1 >= 2*sqrt(pm*p0*qm*q0) <= pm*q0 + p0*qm.
        # This gives: 2*p1*q1 >= pm*q0 + p0*qm? NO! The AM-GM says
        # pm*q0 + p0*qm >= 2*sqrt(pm*p0*qm*q0). Both sides are bounded below
        # by 2*sqrt(...). But that means AM-GM gives pm*q0+p0*qm >= 2*sqrt(...)
        # and LC gives p1*q1 >= sqrt(...). So 2*p1*q1 >= 2*sqrt(...) <= pm*q0+p0*qm.
        # The inequality goes the WRONG WAY: both 2*p1*q1 and pm*q0+p0*qm are >= 2*sqrt(...).
        # We can't conclude either dominates the other.

        # WAIT. Actually, from LC:
        # p1^2 >= pm*p0 and q1^2 >= qm*q0.
        # Multiply: (p1*q1)^2 >= pm*p0*qm*q0.
        # AM-GM on {pm*q0, p0*qm}: pm*q0+p0*qm >= 2*sqrt(pm*q0*p0*qm) = 2*sqrt(pm*p0*qm*q0).
        # So pm*q0+p0*qm >= 2*sqrt(pm*p0*qm*q0).
        # And p1*q1 >= sqrt(pm*p0*qm*q0).
        # So 2*p1*q1 >= 2*sqrt(pm*p0*qm*q0) <= pm*q0+p0*qm.
        # So we have 2*p1*q1 >= pm*q0+p0*qm? Let me check.
        # Let x = pm*p0*qm*q0. Then 2*p1*q1 >= 2*sqrt(x), and pm*q0+p0*qm >= 2*sqrt(x).
        # Both are >= 2*sqrt(x), but we can't compare them to each other!

        # Example: pm=1,p0=4,qm=4,q0=1 -> pm*q0+p0*qm=1+16=17, 2sqrt(16)=8.
        # p1 >= sqrt(4)=2, q1 >= sqrt(4)=2, so 2*p1*q1 >= 8. But 17 > 8, no help.
        # In fact if p1=2,q1=2: 2*p1*q1=8<17. So 2*p1*q1 < pm*q0+p0*qm is possible!

        # So the pure Cauchy-Schwarz approach on the cross term FAILS.
        # This explains why cross_neg = 0 but we can't prove it from LC alone.

        # NEW IDEA: The cross term 2*p1*q1 - pm*q0 - p0*qm is actually the
        # "LC surplus of the mixed sequence at positions (m-2, m-1, m)" where
        # the mixed sequence is sqrt(P)*sqrt(Q) or similar.
        # Actually, 2*p1*q1 - pm*q0 - p0*qm IS the cross term in the expansion
        # (p1+q1)^2 - (pm+qm)*(p0+q0) = lc_P + lc_Q + cross.
        # It equals 2*p1*q1 - pm*q0 - p0*qm.
        # The Cauchy-Schwarz inequality for this specific cross term is
        # related to the FKG inequality or Harris inequality on the tree.
        # Since P and Q come from the SAME tree (same rooted subtrees),
        # they are positively correlated at the coefficient level?
        # This might follow from the FKG lattice structure of independent sets.

        # Alternative: use the Polya frequency / total positivity of the sequence.
        # For PF_2 sequences (LC = totally positive of order 2):
        # The PRODUCT of two PF_2 sequences is PF_2.
        # The SUM is not necessarily PF_2, but the determinant condition...

        # Simplest check: the key quantities at each n.

        "residual_min": None,
        "residual_min_witness": None,
        "cross_min": None,
        "cross_min_witness": None,

        # Is cross >= -mismatch always?
        # cross + mismatch = 2*p1*q1 - pm*q0 - p0*qm + p0*q1 - p1*q0
        # = residual.
        # So cross >= -mismatch iff residual >= 0.
        # This is the same as what we're trying to prove!

        # APPROACH B: DIRECT proof that residual >= 0.
        # residual = p1*(2*q1 - q0) + p0*(q1 - qm) - pm*q0
        #
        # Case 1: q1 >= q0 (no Q-drop). Then 2*q1-q0 >= q1 >= 0.
        #   residual >= p1*q1 + p0*(q1-qm) - pm*q0
        #   If q1 >= qm (mode(Q) >= m-1): residual >= p1*q1 - pm*q0.
        #   Need: p1*q1 >= pm*q0.
        #   From LC: p1^2 >= pm*p0 and q1^2 >= qm*q0 >= q0^2 (if qm >= q0, but not guaranteed).
        #   Actually q1 >= qm doesn't help with q1 vs q0 directly.
        #   But p1*q1 vs pm*q0: by LC of B at m-1, b1^2 >= b2*b0 = (pm+qm)*(p0+q0).
        #   So (p1+q1)^2 >= (pm+qm)*(p0+q0). But we need p1*q1 >= pm*q0, which is
        #   a specific cross term. Not directly implied.

        # NEW APPROACH C: USE Q = x * P' structure.
        # Q = dp[u][1] = x * product dp0[c] =: x*P'.
        # So q_k = P'_{k-1}.
        # q0 = q_{m-2} = P'_{m-3}
        # q1 = q_{m-1} = P'_{m-2}
        # qm = q_m = P'_{m-1}
        #
        # residual = p1*(2*P'_{m-2} - P'_{m-3}) + p0*(P'_{m-2} - P'_{m-1}) - pm*P'_{m-3}
        #
        # P' is also a product of LC sequences (product of dp0[c]).
        # So P' is LC.
        # P'_{m-2}^2 >= P'_{m-3}*P'_{m-1}.
        # If mode(P') >= m-2: P'_{m-2} >= P'_{m-3} (non-decreasing up to mode).
        # Then 2*P'_{m-2} - P'_{m-3} >= P'_{m-2} >= 0.

        # What is mode(P') relative to m-2? Check this.
        "mode_Pprime_ge_m2": 0,
        "mode_Pprime_ge_m1": 0,
        "mode_Pprime_lt_m2": 0,
        "mode_Pprime_eq_m2": 0,

        # If mode(P') >= m-2 and mode(P') = m-2:
        # P'_{m-2} >= P'_{m-1}, so q1-qm = P'_{m-2}-P'_{m-1} >= 0.
        # And P'_{m-2} >= P'_{m-3}, so 2P'_{m-2}-P'_{m-3} >= P'_{m-2}.
        # residual >= p1*P'_{m-2} + p0*0 - pm*P'_{m-3}
        #          = p1*P'_{m-2} - pm*P'_{m-3}
        # By LC of P: p1^2 >= pm*p0.
        # By LC of P': P'_{m-2}^2 >= P'_{m-3}*P'_{m-1} >= P'_{m-3}^2 (if P'_{m-1}>=P'_{m-3}).
        # Hmm, need more.

        # If mode(P') >= m-1:
        # P'_{m-1} >= P'_{m-2} >= P'_{m-3} (non-decreasing up to m-1).
        # q1 - qm = P'_{m-2} - P'_{m-1} <= 0.
        # 2P'_{m-2} - P'_{m-3} >= P'_{m-2} >= 0 (since P'_{m-2} >= P'_{m-3}).
        # residual = p1*(2P'_{m-2}-P'_{m-3}) + p0*(P'_{m-2}-P'_{m-1}) - pm*P'_{m-3}
        # The middle term is <= 0 and the last term is <= 0.
        # Need the first term to dominate.

        # APPROACH D: Schur convexity / majorization.
        # (p_0, p_1, ...) and (P'_0, P'_1, ...) are both LC sequences.
        # The determinant p1*P'_{m-3} - pm*P'_{m-3} ... no, this isn't clean.

        # Let me just collect the data and find what works.

        # RATIO TEST: For Q-drop cases, what is p1*delta_b / (b1*|delta_q|)?
        "qdrop_ratio_min": None,
        "qdrop_ratio_max": None,

        # For mismatch-neg cases without Q-drop:
        # rise - neg = p1*delta_b + b1*delta_q
        # Both terms >= 0, so trivially >= 0.
        "mismatch_neg_no_qdrop": 0,
        "mismatch_neg_qdrop": 0,

        # P' mode relative to m
        "pprime_mode_values": {},  # (mode(P') - (m-2)) -> count

        # Approach E: bound residual via P' LC.
        # residual = p1*(2*q1 - q0) + p0*(q1 - qm) - pm*q0
        # = p1*(2*P'_{m-2} - P'_{m-3}) + p0*(P'_{m-2} - P'_{m-1}) - pm*P'_{m-3}
        #
        # Factor using P' LC: P'_{m-2}^2 >= P'_{m-3}*P'_{m-1}
        # So P'_{m-1} <= P'_{m-2}^2/P'_{m-3} (when P'_{m-3}>0).
        #
        # p0*(P'_{m-2} - P'_{m-1}) >= p0*(P'_{m-2} - P'_{m-2}^2/P'_{m-3})
        #                           = p0*P'_{m-2}*(1 - P'_{m-2}/P'_{m-3})
        #                           = p0*P'_{m-2}*(P'_{m-3}-P'_{m-2})/P'_{m-3}
        # When P' is increasing (P'_{m-2} >= P'_{m-3}): this gives p0*negative, bad.
        # When P' is decreasing: this gives p0*positive, good.
        # So this bound is useful only when mode(P') <= m-2.

        # Let me also track the exact (p,q) values for mismatch-neg trees.
        "mismatch_neg_witnesses": [],
    }

    t0 = time.time()

    for n in range(args.min_n, args.max_n + 1):
        for res in range(args.mod):
            cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
            if args.mod > 1:
                cmd.append(f"{res}/{args.mod}")
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            assert proc.stdout is not None

            for raw in proc.stdout:
                nn, adj = parse_graph6(raw)
                if not is_dleaf_le_1(nn, adj):
                    continue

                poly_t = independence_poly(nn, adj)
                m = mode_index_leftmost(poly_t)
                if m == 0 or m >= len(poly_t) or poly_t[m] == 0:
                    continue

                deg = [len(nb) for nb in adj]
                leaves = [v for v in range(nn) if deg[v] == 1]
                min_parent_deg = min(deg[adj[l][0]] for l in leaves)
                leaf = min(l for l in leaves if deg[adj[l][0]] == min_parent_deg)
                support = adj[leaf][0]

                if deg[support] != 2:
                    continue

                u = [x for x in adj[support] if x != leaf][0]
                b_adj = remove_vertices(adj, {leaf, support})
                b_poly = independence_poly(len(b_adj), b_adj)

                if m - 2 < 0 or m >= len(b_poly) or m - 1 >= len(b_poly) or b_poly[m - 1] == 0:
                    continue

                keep = [v for v in range(nn) if v not in {leaf, support}]
                idx_map = {v: i for i, v in enumerate(keep)}
                u_in_b = idx_map[u]
                p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

                stats["checked"] += 1

                p0 = p_poly[m - 2] if m - 2 < len(p_poly) else 0
                p1 = p_poly[m - 1] if m - 1 < len(p_poly) else 0
                pm = p_poly[m] if m < len(p_poly) else 0
                q0 = q_poly[m - 2] if m - 2 < len(q_poly) else 0
                q1 = q_poly[m - 1] if m - 1 < len(q_poly) else 0
                qm = q_poly[m] if m < len(q_poly) else 0
                b0 = b_poly[m - 2]
                b1 = b_poly[m - 1]
                b2 = b_poly[m]

                mismatch = p0 * b1 - p1 * b0
                lc_surplus = b1 * b1 - b2 * b0
                combined = lc_surplus + mismatch
                neg = -mismatch
                delta_b = b1 - b0
                delta_q = q1 - q0
                rise = b1 * delta_b

                lc_P = p1 * p1 - pm * p0
                lc_Q = q1 * q1 - qm * q0
                cross = 2 * p1 * q1 - pm * q0 - p0 * qm
                residual = cross + mismatch  # = cross - neg_mismatch = cross + p0*q1 - p1*q0

                if combined < 0:
                    stats["combined_neg"] += 1
                if rise - neg < 0:
                    stats["rise_fail"] += 1

                if stats["residual_min"] is None or residual < stats["residual_min"]:
                    stats["residual_min"] = residual
                    stats["residual_min_witness"] = {
                        "n": nn, "m": m,
                        "p0": p0, "p1": p1, "pm": pm,
                        "q0": q0, "q1": q1, "qm": qm,
                        "b0": b0, "b1": b1, "b2": b2,
                        "residual": residual, "lc_P": lc_P, "lc_Q": lc_Q, "cross": cross,
                        "g6": raw.decode("ascii").strip()
                    }

                if stats["cross_min"] is None or cross < stats["cross_min"]:
                    stats["cross_min"] = cross

                # P' mode
                # Q = x*P', so P'_k = q_{k+1}
                # P' coefficients: p'_j = q_{j+1}
                pprime = q_poly[1:] if len(q_poly) > 1 else [0]
                mode_pprime = mode_index_leftmost(pprime) if pprime else 0
                shift_pp = mode_pprime - (m - 2)
                key_pp = str(shift_pp)
                stats["pprime_mode_values"][key_pp] = stats["pprime_mode_values"].get(key_pp, 0) + 1

                if mode_pprime >= m - 2:
                    stats["mode_Pprime_ge_m2"] += 1
                if mode_pprime >= m - 1:
                    stats["mode_Pprime_ge_m1"] += 1
                if mode_pprime < m - 2:
                    stats["mode_Pprime_lt_m2"] += 1
                if mode_pprime == m - 2:
                    stats["mode_Pprime_eq_m2"] += 1

                if mismatch < 0:
                    stats["mismatch_neg"] += 1
                    q_drop = (delta_q < 0)
                    if q_drop:
                        stats["mismatch_neg_qdrop"] += 1
                        stats["q_drop"] += 1
                        if delta_b > 0:
                            ratio = (p1 * delta_b) / (b1 * (-delta_q))
                            if stats["qdrop_ratio_min"] is None or ratio < stats["qdrop_ratio_min"]:
                                stats["qdrop_ratio_min"] = ratio
                            if stats["qdrop_ratio_max"] is None or ratio > stats["qdrop_ratio_max"]:
                                stats["qdrop_ratio_max"] = ratio
                    else:
                        stats["mismatch_neg_no_qdrop"] += 1
                        # Here rise - neg = p1*delta_b + b1*delta_q >= 0 trivially.
                        # But wait, delta_q >= 0 and delta_b >= 0 and p1 >= 0 and b1 >= 0,
                        # so rise - neg >= 0. But mismatch < 0 means neg > 0.
                        # Can neg > rise when delta_q >= 0?
                        # rise - neg = p1*delta_b + b1*delta_q >= 0. YES, trivially!

                    if n <= 23 and len(stats["mismatch_neg_witnesses"]) < 200:
                        w = {
                            "n": nn, "m": m, "g6": raw.decode("ascii").strip(),
                            "p0": p0, "p1": p1, "pm": pm,
                            "q0": q0, "q1": q1, "qm": qm,
                            "b0": b0, "b1": b1, "b2": b2,
                            "mismatch": mismatch, "lc_surplus": lc_surplus,
                            "combined": combined, "neg": neg, "rise": rise,
                            "delta_b": delta_b, "delta_q": delta_q,
                            "lc_P": lc_P, "lc_Q": lc_Q, "cross": cross,
                            "residual": residual,
                            "mode_P": mode_index_leftmost(p_poly),
                            "mode_Pprime": mode_pprime,
                            "q_drop": q_drop,
                        }
                        stats["mismatch_neg_witnesses"].append(w)

            proc.wait()

        print(f"n={n:2d}: checked={stats['checked']:8d} mismatch-={stats['mismatch_neg']:5d} "
              f"qdrop={stats['q_drop']} rise_fail={stats['rise_fail']} "
              f"combined-={stats['combined_neg']} resid_min={stats['residual_min']}",
              flush=True)

    stats["wall_s"] = time.time() - t0

    print("\n" + "=" * 80)
    print("STRONG C2 RISE COMPENSATION V2: STRUCTURAL ANALYSIS")
    print("=" * 80)
    print(f"Trees checked: {stats['checked']:,}")
    print(f"Mismatch negative: {stats['mismatch_neg']}")
    print(f"  - with Q-drop: {stats['mismatch_neg_qdrop']}")
    print(f"  - without Q-drop: {stats['mismatch_neg_no_qdrop']} (trivially handled)")
    print(f"Rise failures: {stats['rise_fail']}")
    print(f"Combined negative: {stats['combined_neg']}")
    print()
    print(f"Residual min: {stats['residual_min']}")
    print(f"Cross min: {stats['cross_min']}")
    print()
    print(f"Mode(P') distribution (shift from m-2):")
    for k in sorted(stats["pprime_mode_values"], key=lambda s: int(s)):
        print(f"  shift {k}: {stats['pprime_mode_values'][k]}")
    print(f"mode(P') >= m-2: {stats['mode_Pprime_ge_m2']}")
    print(f"mode(P') >= m-1: {stats['mode_Pprime_ge_m1']}")
    print(f"mode(P') < m-2: {stats['mode_Pprime_lt_m2']}")
    print()
    if stats["qdrop_ratio_min"] is not None:
        print(f"Q-drop ratio p1*db/(b1*|dq|): min={stats['qdrop_ratio_min']:.6f} max={stats['qdrop_ratio_max']:.6f}")
    print(f"Time: {stats['wall_s']:.1f}s")

    # Print mismatch witnesses for analysis
    print("\n--- Mismatch-negative witnesses ---")
    for w in stats["mismatch_neg_witnesses"]:
        print(f"  n={w['n']} m={w['m']} qdrop={w['q_drop']} "
              f"mismatch={w['mismatch']} lc_surplus={w['lc_surplus']} combined={w['combined']} "
              f"rise={w['rise']} neg={w['neg']} "
              f"p0={w['p0']} p1={w['p1']} q0={w['q0']} q1={w['q1']} "
              f"lc_P={w['lc_P']} lc_Q={w['lc_Q']} cross={w['cross']} resid={w['residual']} "
              f"mode_P={w['mode_P']} mode_P'={w['mode_Pprime']}")

    # Save
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    # Don't save witnesses to JSON (too verbose)
    save_stats = {k: v for k, v in stats.items() if k != "mismatch_neg_witnesses"}
    save_stats["mismatch_neg_witness_count"] = len(stats["mismatch_neg_witnesses"])
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(save_stats, f, indent=2)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
