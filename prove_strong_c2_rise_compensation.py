#!/usr/bin/env python3
"""Algebraic verification for the STRONG C2 rise-compensation inequality.

Target: prove lc_surplus + mismatch >= 0 where
  lc_surplus = b_{m-1}^2 - b_m * b_{m-2}
  mismatch   = p_{m-2} * b_{m-1} - p_{m-1} * b_{m-2}

Equivalently: neg <= rise where
  neg  = p_{m-1} * b_{m-2} - p_{m-2} * b_{m-1}   (when mismatch < 0)
  rise = b_{m-1} * (b_{m-1} - b_{m-2})

Key identity (always true):
  rise - neg = p_{m-1} * (b_{m-1} - b_{m-2}) + b_{m-1} * (q_{m-1} - q_{m-2})

So neg <= rise iff:
  p_{m-1} * (b_{m-1} - b_{m-2}) + b_{m-1} * (q_{m-1} - q_{m-2}) >= 0

This script explores:
1. The ratio p_k/b_k and its relationship to LC of P vs B
2. Whether lc_surplus >= b_{m-2} * (b_{m-1} - b_{m-2}) always holds
   (this would give lc_surplus >= neg immediately via the "rate bound")
3. The "mixed LC" inequality: p_{m-2}/p_{m-1} <= b_{m-2}/b_{m-1}
   (which is exactly mismatch >= 0)
4. The refined approach: showing that when mismatch < 0,
   the LC surplus compensates via factoring through b_{m-1} - b_{m-2}.

Approach 4 (the cleanest path):

Define delta_b = b_{m-1} - b_{m-2} >= 0 (since mode(B) >= m-1 verified).
Define delta_q = q_{m-1} - q_{m-2}.
Define delta_p = p_{m-1} - p_{m-2}.

Note: delta_b = delta_p + delta_q (since b_k = p_k + q_k).

The identity says:
  rise - neg = p_{m-1} * delta_b + b_{m-1} * delta_q
             = p_{m-1} * (delta_p + delta_q) + b_{m-1} * delta_q
             = p_{m-1} * delta_p + (p_{m-1} + b_{m-1}) * delta_q

Also:
  lc_surplus + mismatch = b_{m-1}^2 - b_m * b_{m-2} + p_{m-2} * b_{m-1} - p_{m-1} * b_{m-2}
  = b_{m-1} * (b_{m-1} + p_{m-2}) - b_{m-2} * (b_m + p_{m-1})
  = b_{m-1} * (b_{m-1} + p_{m-2}) - b_{m-2} * (b_m + p_{m-1})

Using b_k = p_k + q_k:
  = (p_{m-1}+q_{m-1}) * (p_{m-1}+q_{m-1}+p_{m-2}) - (p_{m-2}+q_{m-2}) * (p_m+q_m+p_{m-1})

This script verifies candidate structural inequalities on all d_leaf<=1 trees.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from fractions import Fraction
from typing import Any

from conjecture_a_hall_subset_scan import is_dleaf_le_1
from diagnose_bridge_decomposition import compute_hub_polys, mode_index_leftmost, remove_vertices
from graph6 import parse_graph6
from indpoly import independence_poly, is_log_concave


def fresh_stats() -> dict[str, Any]:
    return {
        "checked": 0,
        "mismatch_neg": 0,
        "combined_neg": 0,
        "rise_fail": 0,
        # Approach 1: p_k/b_k ratio monotonicity
        "ratio_p_over_b_increasing_at_m1": 0,  # p_{m-1}/b_{m-1} >= p_{m-2}/b_{m-2}
        "ratio_p_over_b_decreasing_at_m1": 0,
        # Approach 2: lc_surplus >= b_{m-2} * delta_b
        "lc_dominates_quadratic_fail": 0,
        "min_lc_minus_quad": None,
        # Approach 3: factored form
        # lc_surplus = b_{m-1}^2 - b_m*b_{m-2}
        #            = (b_{m-1}-b_{m-2})(b_{m-1}+b_{m-2}) + b_{m-2}(b_{m-1}+b_{m-2}) - b_m*b_{m-2}
        # Actually: lc_surplus = (b_{m-1}-b_{m-2})^2 + b_{m-2}*(2b_{m-1}-b_{m-2}-b_m)
        # Hmm, let's try: lc_surplus = delta_b * b_{m-1} + b_{m-2} * (b_{m-1} - b_m)
        # = delta_b * b_{m-1} - b_{m-2} * (b_m - b_{m-1})
        # when mode(B) >= m-1, b_{m-1} >= b_m, so second term >= 0.
        # So lc_surplus >= delta_b * b_{m-1} >= 0.
        "lc_ge_deltab_times_bm1": 0,  # always?
        "lc_ge_deltab_times_bm1_fail": 0,
        # Approach 4: neg <= delta_b * b_{m-1} (sufficient since lc_surplus >= delta_b * b_{m-1})
        # neg = p_{m-1}*b_{m-2} - p_{m-2}*b_{m-1}
        # neg <= delta_b * b_{m-1}?
        # i.e. p_{m-1}*b_{m-2} - p_{m-2}*b_{m-1} <= (b_{m-1}-b_{m-2})*b_{m-1}
        # i.e. p_{m-1}*b_{m-2} <= b_{m-1}^2 - b_{m-2}*b_{m-1} + p_{m-2}*b_{m-1}
        # i.e. p_{m-1}*b_{m-2} <= b_{m-1}*(b_{m-1} - b_{m-2} + p_{m-2})
        # This is "neg <= rise", which is exactly what we're proving.
        # So approach 4 = approach "rise".

        # Approach 5: factored LC surplus
        # lc_surplus = delta_b * b_{m-1} + b_{m-2} * (b_{m-1} - b_m)
        #            = delta_b * b_{m-1} + b_{m-2} * drop_b
        # where drop_b = b_{m-1} - b_m >= 0 (mode(B)>=m-1 regime: shift=0).
        # In shift=1 regime (mode(B)=m), b_m >= b_{m-1}, so drop_b < 0,
        # but delta_b is larger.
        # mismatch = p_{m-2}*b_{m-1} - p_{m-1}*b_{m-2}
        #          = b_{m-1}*p_{m-2} - b_{m-2}*p_{m-1}
        # For mismatch >= 0 we need p_{m-2}/p_{m-1} >= b_{m-2}/b_{m-1},
        # i.e., P's ratio at (m-2,m-1) >= B's ratio.
        # For LC sequences P and B with mode(P) <= m-2 or mode(P) >= m-1,
        # this relates to where we sit relative to their modes.

        # KEY APPROACH: decompose via modes
        # mode(B) >= m-1 (verified 0 fails)
        # mode(P) -- what is it relative to m-2?
        "mode_P_ge_m2": 0,
        "mode_P_lt_m2": 0,
        "mode_P_eq_m2": 0,
        "mode_P_eq_m1": 0,
        "mode_P_ge_m1": 0,

        # When mode(P) >= m-1: P is increasing at m-2, so p_{m-1} >= p_{m-2}
        # and B is increasing at m-2 (since mode(B) >= m-1), so b_{m-1} >= b_{m-2}
        # Mismatch = p_{m-2}*b_{m-1} - p_{m-1}*b_{m-2}: can it be negative?
        # Only if p_{m-1}/p_{m-2} > b_{m-1}/b_{m-2}.

        # When mode(P) = m-2: P peaks at m-2, so p_{m-2} >= p_{m-1}
        # Then mismatch = p_{m-2}*b_{m-1} - p_{m-1}*b_{m-2} >= 0 since both factors favor it.
        # (p_{m-2} >= p_{m-1} and b_{m-1} >= b_{m-2})
        # This would IMMEDIATELY show mismatch >= 0 when mode(P) <= m-2!
        "mismatch_neg_and_mode_P_le_m2": 0,  # should be 0 by above argument!

        # Approach 6 (MOST PROMISING): When mismatch < 0, mode(P) >= m-1.
        # Then we use LC of both P and B.
        # Both LC at index m-1:
        #   p_{m-1}^2 >= p_{m-2}*p_m      ...(P-LC)
        #   b_{m-1}^2 >= b_{m-2}*b_m      ...(B-LC = lc_surplus >= 0)
        #
        # mismatch < 0 means p_{m-1}*b_{m-2} > p_{m-2}*b_{m-1}
        # i.e. p_{m-1}/b_{m-1} > p_{m-2}/b_{m-2}
        # i.e. the P-fraction is RISING at (m-2, m-1).
        #
        # Factor: lc_surplus + mismatch
        # = b_{m-1}^2 - b_m*b_{m-2} + p_{m-2}*b_{m-1} - p_{m-1}*b_{m-2}
        # = b_{m-1}*(b_{m-1}+p_{m-2}) - b_{m-2}*(b_m+p_{m-1})
        # = b_{m-1}*a_{m-1} - b_{m-2}*a_m          ... where a_k = b_k + p_{k-1}
        # = b_{m-1}*a_{m-1} - b_{m-2}*a_m
        #
        # This is positive iff a_{m-1}/a_m >= b_{m-2}/b_{m-1}
        # i.e. lambda_m(A) >= lambda_{m-1}(B)  (the STRONG C2 condition!)
        #
        # Approach 7: use that a_k = b_k + p_{k-1}, and A = I(B)(1+x) + xP,
        # so A inherits LC from B and P.

        # Check: is A always LC at index m?
        "A_lc_at_m": 0,
        "A_not_lc_at_m": 0,

        # Check: a_{m-1}^2 >= a_{m-2}*a_m?
        "A_lc_at_m1": 0,
        "A_not_lc_at_m1": 0,

        # The Cauchy-Schwarz approach:
        # If B is LC and P is LC, is (B + xP) LC?
        # (B + xP)_k = b_k + p_{k-1}
        # LC at k means (b_k+p_{k-1})^2 >= (b_{k-1}+p_{k-2})*(b_{k+1}+p_k)
        # Expanding: b_k^2 + 2*b_k*p_{k-1} + p_{k-1}^2
        #   >= b_{k-1}*b_{k+1} + b_{k-1}*p_k + p_{k-2}*b_{k+1} + p_{k-2}*p_k
        # Using B-LC: b_k^2 >= b_{k-1}*b_{k+1}, so we need:
        #   2*b_k*p_{k-1} + p_{k-1}^2 >= b_{k-1}*p_k + p_{k-2}*b_{k+1} + p_{k-2}*p_k
        # Using P-LC: p_{k-1}^2 >= p_{k-2}*p_k, so we need:
        #   2*b_k*p_{k-1} >= b_{k-1}*p_k + p_{k-2}*b_{k+1}
        # This is 2*b_k*p_{k-1} >= b_{k-1}*p_k + p_{k-2}*b_{k+1}
        # = p_k*b_{k-1} + p_{k-2}*b_{k+1}
        # By AM-GM on the cross terms? Not obvious.
        # But this is related to the "interleaving" of B and P growth rates.

        # Approach 8 (DIRECT CAUCHY-SCHWARZ):
        # lc_surplus + mismatch
        # = b_{m-1}^2 - b_m*b_{m-2} + p_{m-2}*b_{m-1} - p_{m-1}*b_{m-2}
        # Factor differently:
        # = b_{m-1}*(b_{m-1} - b_{m-2}) + b_{m-2}*(b_{m-1} - b_m) + p_{m-2}*b_{m-1} - p_{m-1}*b_{m-2}
        # = delta_b * b_{m-1} + b_{m-2}*drop + p_{m-2}*b_{m-1} - p_{m-1}*b_{m-2}
        # = b_{m-1}*(delta_b + p_{m-2}) - b_{m-2}*(p_{m-1} - drop)
        # where drop = b_{m-1} - b_m = -(b_m - b_{m-1}).
        # When mode(B) >= m-1 and mode(B) = m-1 (shift=0): drop >= 0.
        # Then: = b_{m-1}*(delta_b + p_{m-2}) - b_{m-2}*(p_{m-1} - drop)
        #       = b_{m-1}*(delta_b + p_{m-2}) - b_{m-2}*(p_{m-1} - b_{m-1} + b_m)

        # Let's try yet another factoring.
        # lc_surplus + mismatch = b_{m-1}*a_{m-1} - b_{m-2}*a_m
        # where a_k = b_k + p_{k-1}.
        # = b_{m-1}*(b_{m-1}+p_{m-2}) - b_{m-2}*(b_m+p_{m-1})
        #
        # Both A and B are (conjecturally) LC. So:
        # a_{m-1}^2 >= a_{m-2}*a_m  ... if A is LC at m-1
        # gives: a_{m-1}/a_m >= a_{m-2}/a_{m-1}
        # We want: a_{m-1}/a_m >= b_{m-2}/b_{m-1}
        # i.e. a_{m-1}*b_{m-1} >= a_m*b_{m-2}
        #
        # If a_{m-2}/a_{m-1} >= b_{m-2}/b_{m-1}, then LC of A gives:
        # a_{m-1}/a_m >= a_{m-2}/a_{m-1} >= b_{m-2}/b_{m-1} and we're done.
        #
        # So check: a_{m-2}/a_{m-1} >= b_{m-2}/b_{m-1}?
        # i.e. a_{m-2}*b_{m-1} >= a_{m-1}*b_{m-2}?
        # a_{m-2} = b_{m-2} + p_{m-3}
        # a_{m-1} = b_{m-1} + p_{m-2}
        # (b_{m-2}+p_{m-3})*b_{m-1} >= (b_{m-1}+p_{m-2})*b_{m-2}
        # b_{m-2}*b_{m-1} + p_{m-3}*b_{m-1} >= b_{m-1}*b_{m-2} + p_{m-2}*b_{m-2}
        # p_{m-3}*b_{m-1} >= p_{m-2}*b_{m-2}
        # This is "mismatch at (m-3, m-2)" for B.
        "lower_mismatch_neg": 0,  # p_{m-3}*b_{m-1} < p_{m-2}*b_{m-2}
        "lower_mismatch_pos": 0,

        # Verify the identity: rise - neg = p1*delta_b + b1*delta_q
        "identity_fail": 0,

        # Approach 9: LC-surplus decomposition via P and Q LC surpluses
        # b_{m-1}^2 - b_m*b_{m-2} = (p_{m-1}+q_{m-1})^2 - (p_m+q_m)*(p_{m-2}+q_{m-2})
        # = [p_{m-1}^2 - p_m*p_{m-2}] + [q_{m-1}^2 - q_m*q_{m-2}]
        #   + [2*p_{m-1}*q_{m-1} - p_m*q_{m-2} - p_{m-2}*q_m]
        # = lc_P + lc_Q + cross
        # where cross = 2*p_{m-1}*q_{m-1} - p_m*q_{m-2} - p_{m-2}*q_m
        # By Cauchy-Schwarz / AM-GM on (p_{m-1},q_{m-1}) vs (p_m,q_{m-2}) and (p_{m-2},q_m):
        # cross = p_{m-1}*q_{m-1} - p_m*q_{m-2} + p_{m-1}*q_{m-1} - p_{m-2}*q_m
        # Each of these: p_{m-1}*q_{m-1} >= sqrt(p_m*p_{m-2}) * sqrt(q_m*q_{m-2})
        # by LC... >= sqrt(p_m*q_{m-2} * p_{m-2}*q_m)  by rearrangement.
        # Hmm, need: p_{m-1}*q_{m-1} >= sqrt(p_m*q_{m-2})*sqrt(p_{m-2}*q_m)?
        # = sqrt(p_m*p_{m-2}*q_m*q_{m-2})
        # By LC: p_{m-1}^2 >= p_m*p_{m-2} and q_{m-1}^2 >= q_m*q_{m-2}
        # So p_{m-1}*q_{m-1} >= sqrt(p_0*p_{m-2})*sqrt(q_m*q_{m-2})
        # Yes! By AM-GM: X*Y >= sqrt(A*B)*sqrt(C*D) when X^2 >= A*B and Y^2 >= C*D.
        # Actually that gives X*Y >= sqrt(AB)*sqrt(CD) = sqrt(ABCD).
        # And we need p_{m-1}*q_{m-1} - p_m*q_{m-2} >= 0 and same for the other.
        # By CBS: (p_{m-1}*q_{m-1})^2 >= p_m*p_{m-2} * q_m*q_{m-2}
        # But p_m*q_{m-2} vs p_{m-1}*q_{m-1} is not directly handled.
        # Let r = p_m/p_{m-1}, s = q_{m-2}/q_{m-1} (both <= 1 when in rising part).
        # Then p_m*q_{m-2} / (p_{m-1}*q_{m-1}) = r*s.
        # cross >= 0 iff r*s + (1/r)*(1/s) <= 2 * ... no, that's not right.

        # Let me just check cross >= 0 empirically.
        "cross_neg": 0,
        "cross_pos": 0,

        # Approach 10: DIRECT bound on -mismatch
        # -mismatch = p_{m-1}*b_{m-2} - p_{m-2}*b_{m-1}
        #           = p_{m-1}*(p_{m-2}+q_{m-2}) - p_{m-2}*(p_{m-1}+q_{m-1})
        #           = p_{m-1}*q_{m-2} - p_{m-2}*q_{m-1}
        # So -mismatch = p_{m-1}*q_{m-2} - p_{m-2}*q_{m-1}
        # And lc_surplus = lc_P + lc_Q + cross (from approach 9)
        # So lc_surplus + mismatch = lc_P + lc_Q + cross - (p_{m-1}*q_{m-2} - p_{m-2}*q_{m-1})
        # = lc_P + lc_Q + (2*p_{m-1}*q_{m-1} - p_m*q_{m-2} - p_{m-2}*q_m)
        #   - (p_{m-1}*q_{m-2} - p_{m-2}*q_{m-1})
        # = lc_P + lc_Q + 2*p_{m-1}*q_{m-1} - p_m*q_{m-2} - p_{m-2}*q_m
        #   - p_{m-1}*q_{m-2} + p_{m-2}*q_{m-1}
        # = lc_P + lc_Q + p_{m-1}*q_{m-1} + p_{m-2}*q_{m-1}
        #   + p_{m-1}*q_{m-1} - p_m*q_{m-2} - p_{m-2}*q_m - p_{m-1}*q_{m-2}
        # Hmm let me redo this more carefully in the code.

        # THE KEY SIMPLIFICATION:
        # -mismatch = p_{m-1}*q_{m-2} - p_{m-2}*q_{m-1}
        # This is the "cross-LC" of (P,Q) at indices (m-2, m-1).
        # It measures whether P grows faster than Q near m-1.
        "cross_pq_neg": 0,  # p_{m-1}*q_{m-2} - p_{m-2}*q_{m-1} < 0
        "cross_pq_pos": 0,  # = -mismatch > 0 means P grows faster than Q

        # Verify: -mismatch = p_{m-1}*q_{m-2} - p_{m-2}*q_{m-1}
        "neg_mismatch_identity_fail": 0,

        # Now: lc_surplus = lc_P + lc_Q + cross  (approach 9)
        # And we need lc_P + lc_Q + cross >= -mismatch = p_{m-1}*q_{m-2} - p_{m-2}*q_{m-1}
        # i.e. lc_P + lc_Q + cross + p_{m-2}*q_{m-1} - p_{m-1}*q_{m-2} >= 0
        # = lc_P + lc_Q + (2*p1*q1 - pm*q0 - p0*qm) + p0*q1 - p1*q0
        # where p0=p_{m-2}, p1=p_{m-1}, pm=p_m, q0=q_{m-2}, q1=q_{m-1}, qm=q_m
        # = lc_P + lc_Q + 2*p1*q1 - pm*q0 - p0*qm + p0*q1 - p1*q0
        # = lc_P + lc_Q + p1*q1 + p0*q1 + p1*q1 - pm*q0 - p0*qm - p1*q0
        # = lc_P + lc_Q + q1*(p1+p0) + p1*q1 - q0*(pm+p1) - p0*qm
        # Hmm, this is getting complex. Let me just check things numerically.

        # NEW APPROACH 11:
        # -mismatch = p1*q0 - p0*q1 (using shorthand p0=p_{m-2}, etc.)
        # lc_surplus = (p1+q1)^2 - (pm+qm)*(p0+q0)
        #
        # We want: (p1+q1)^2 - (pm+qm)*(p0+q0) >= p1*q0 - p0*q1
        #
        # LHS = p1^2 + q1^2 + 2*p1*q1 - pm*p0 - pm*q0 - qm*p0 - qm*q0
        # RHS = p1*q0 - p0*q1
        #
        # LHS - RHS = p1^2 + q1^2 + 2*p1*q1 - pm*p0 - pm*q0 - qm*p0 - qm*q0 - p1*q0 + p0*q1
        # = (p1^2 - pm*p0) + (q1^2 - qm*q0) + (2*p1*q1 - pm*q0 - qm*p0 - p1*q0 + p0*q1)
        # = lc_P + lc_Q + (2*p1*q1 - pm*q0 - qm*p0 - p1*q0 + p0*q1)
        # = lc_P + lc_Q + [p1*(2*q1 - q0) + p0*(q1 - qm) - pm*q0]
        # = lc_P + lc_Q + [p1*(q1 + q1 - q0) + p0*q1 - p0*qm - pm*q0]
        # = lc_P + lc_Q + [p1*q1 + p1*(q1-q0) + p0*q1 - p0*qm - pm*q0]
        # = lc_P + lc_Q + [q1*(p1+p0) + p1*(q1-q0) - p0*qm - pm*q0]
        #
        # Define: residual = q1*(p1+p0) + p1*(q1-q0) - p0*qm - pm*q0
        # When q1 >= q0: residual >= q1*(p1+p0) - p0*qm - pm*q0
        # When q1 < q0 (Q-drop): p1*(q1-q0) < 0 but bounded.
        #
        # Since lc_P and lc_Q are >= 0, we need: lc_P + lc_Q + residual >= 0.
        "min_residual": None,
        "residual_neg": 0,
        "lc_P_plus_lc_Q_plus_residual_neg": 0,

        # Check specifically: lc_P alone >= |residual| when residual < 0?
        "lc_P_ge_neg_residual": 0,
        "lc_P_lt_neg_residual": 0,

        # CLEANEST APPROACH (Approach 12):
        # Since -mismatch = p1*q0 - p0*q1, and using q_k = b_k - p_k:
        # -mismatch = p1*(b0-p0) - p0*(b1-p1) = p1*b0 - p1*p0 - p0*b1 + p0*p1
        # = p1*b0 - p0*b1   (the p0*p1 terms cancel!)
        # This confirms -mismatch = p1*b0 - p0*b1 (as expected).
        # But also -mismatch = p1*q0 - p0*q1 -- the cross-product of P and Q sequences!
        #
        # Now P is a product of subtree IS polynomials, Q = x * product_of_dp0s.
        # Q = x * P', where P' = product over children c of dp0[c].
        # Actually Q = dp[u][1] = x * product_{c children of u} dp0[c].
        # And P = dp[u][0] = product_{c children of u} (dp0[c] + dp1[c])
        #       = product_{c} I(subtree_c).
        # So P' (the product of dp0s) is different from P.
        #
        # Let P' = product_{c} dp0[c]. Then Q = xP'.
        # q_k = P'_{k-1} (Q is P' shifted by 1).
        # So q_{m-2} = P'_{m-3}, q_{m-1} = P'_{m-2}.
        # -mismatch = p_{m-1}*P'_{m-3} - p_{m-2}*P'_{m-2}
        #
        # When P' is LC (it's a product of LC sequences):
        # P'_{m-2}^2 >= P'_{m-3}*P'_{m-1}
        # This gives P'_{m-2}/P'_{m-3} >= P'_{m-1}/P'_{m-2}.
        # For -mismatch >= 0: p_{m-1}/p_{m-2} >= P'_{m-2}/P'_{m-3}.
        # Hmm, need to compare P's growth rate to P''s growth rate.
        # P = product of (dp0+dp1) = product of I(T_c)
        # P' = product of dp0
        # The ratio p_k/P'_k measures "how much the dp1 terms contribute at level k".
        # As k increases, dp1 contributes more (it's x * product dp0), so P grows faster.
        # So p_{m-1}/P'_{m-2} >= p_{m-2}/P'_{m-3}? This would give -mismatch >= 0.
        # But we know -mismatch can be negative (129 cases), so this is wrong.
        #
        # Let me just track P' and verify.
        "pprime_computed": 0,

        # APPROACH 13 (simplest sufficient):
        # Recall combined = b1*a1 - b0*a2  where a_k = b_k + p_{k-1}, b_k = IS(B) coeffs.
        # = b1*(b1+p0) - b0*(b2+p1)   [where b2=b_m, p0=p_{m-2}, p1=p_{m-1}]
        # Since A is LC at index m (verified computationally), a1^2 >= a0*a2 where a0=a_{m-2}, a1=a_{m-1}, a2=a_m.
        # So a2 <= a1^2/a0.
        # combined >= b1*a1 - b0*a1^2/a0 = a1*(b1 - b0*a1/a0) = a1*b1*(1 - (b0/b1)*(a1/a0))
        # For this to be >= 0: b0/b1 * a1/a0 <= 1, i.e. a1/a0 <= b1/b0.
        # a1/a0 = (b1+p0)/(b0+p_{m-3})
        # b1/b0 = b_{m-1}/b_{m-2}
        # So need: (b1+p0)/(b0+p_{m-3}) <= b1/b0
        # i.e. b0*(b1+p0) <= b1*(b0+p_{m-3})
        # i.e. b0*p0 <= b1*p_{m-3}
        # i.e. p0/p_{m-3} <= b1/b0 (= growth rate of B at m-1)
        # This compares P's growth rate at (m-3,m-2) to B's growth rate at (m-2,m-1).
        # Need to check...
        "growth_comparison_fail": 0,

        "wall_s": 0.0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=20)
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--out", default="results/prove_strong_c2_rise_n20.json")
    args = ap.parse_args()

    stats = fresh_stats()
    t0 = time.time()

    for n in range(args.min_n, args.max_n + 1):
        cmd = [args.geng, "-q", str(n), f"{n-1}:{n-1}", "-c"]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert proc.stdout is not None

        n_checked_this = 0
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
            idx = {v: i for i, v in enumerate(keep)}
            u_in_b = idx[u]
            p_poly, q_poly = compute_hub_polys(b_adj, u_in_b)

            stats["checked"] += 1
            n_checked_this += 1

            # Shorthand
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
            rise = b1 * delta_b
            delta_q = q1 - q0

            if combined < 0:
                stats["combined_neg"] += 1

            # Verify identity: rise - neg = p1*delta_b + b1*delta_q
            rise_minus_neg = rise - neg
            identity_rhs = p1 * delta_b + b1 * delta_q
            if rise_minus_neg != identity_rhs:
                stats["identity_fail"] += 1

            # Verify: -mismatch = p1*q0 - p0*q1
            neg_check = p1 * q0 - p0 * q1
            if neg != neg_check:
                stats["neg_mismatch_identity_fail"] += 1

            if mismatch < 0:
                stats["mismatch_neg"] += 1
                if rise - neg < 0:
                    stats["rise_fail"] += 1

            # Mode of P
            mode_P = mode_index_leftmost(p_poly) if p_poly else 0
            if mode_P >= m - 1:
                stats["mode_P_ge_m1"] += 1
            if mode_P == m - 1:
                stats["mode_P_eq_m1"] += 1
            if mode_P == m - 2:
                stats["mode_P_eq_m2"] += 1
            if mode_P >= m - 2:
                stats["mode_P_ge_m2"] += 1
            if mode_P < m - 2:
                stats["mode_P_lt_m2"] += 1

            # KEY CHECK: mismatch negative AND mode(P) <= m-2?
            if mismatch < 0 and mode_P <= m - 2:
                stats["mismatch_neg_and_mode_P_le_m2"] += 1

            # P/B ratio monotonicity
            if b0 > 0 and b1 > 0:
                # p_{m-1}/b_{m-1} vs p_{m-2}/b_{m-2}: cross-multiply
                if p1 * b0 >= p0 * b1:
                    stats["ratio_p_over_b_increasing_at_m1"] += 1
                else:
                    stats["ratio_p_over_b_decreasing_at_m1"] += 1

            # LC surplus >= delta_b * b_{m-1}?
            # lc_surplus = b1^2 - b2*b0 = b1*(b1-b0) + b0*(b1-b2)
            #            = rise + b0*(b1-b2) = rise + b0*drop
            quad = delta_b * b1
            if lc_surplus >= quad:
                stats["lc_ge_deltab_times_bm1"] += 1
            else:
                stats["lc_ge_deltab_times_bm1_fail"] += 1

            # LC surplus >= b0 * delta_b (quadratic bound)?
            quad2 = b0 * delta_b
            if lc_surplus < quad2:
                stats["lc_dominates_quadratic_fail"] += 1
                if stats["min_lc_minus_quad"] is None or (lc_surplus - quad2) < stats["min_lc_minus_quad"]:
                    stats["min_lc_minus_quad"] = lc_surplus - quad2

            # A polynomial and its LC
            a_poly_m2 = b_poly[m - 2] + (p_poly[m - 3] if m - 3 >= 0 and m - 3 < len(p_poly) else 0)
            a_poly_m1 = b1 + p0
            a_poly_m = b2 + p1
            a_lc_at_m = a_poly_m1 * a_poly_m1 >= a_poly_m2 * a_poly_m
            if a_lc_at_m:
                stats["A_lc_at_m"] += 1
            else:
                stats["A_not_lc_at_m"] += 1

            # A LC at m-1
            if m - 3 >= 0:
                a_m3 = (b_poly[m - 3] if m - 3 < len(b_poly) else 0) + (p_poly[m - 4] if m - 4 >= 0 and m - 4 < len(p_poly) else 0)
                a_lc_at_m1 = a_poly_m2 * a_poly_m2 >= a_m3 * a_poly_m1
                if a_lc_at_m1:
                    stats["A_lc_at_m1"] += 1
                else:
                    stats["A_not_lc_at_m1"] += 1

            # LC surplus decomposition (Approach 9)
            lc_P = p1 * p1 - pm * p0
            lc_Q = q1 * q1 - qm * q0
            cross = 2 * p1 * q1 - pm * q0 - p0 * qm
            if cross < 0:
                stats["cross_neg"] += 1
            else:
                stats["cross_pos"] += 1

            # Verify decomposition
            assert lc_P + lc_Q + cross == lc_surplus, f"Decomposition mismatch: {lc_P}+{lc_Q}+{cross} != {lc_surplus}"

            # Residual (Approach 11)
            residual = q1 * (p1 + p0) + p1 * (q1 - q0) - p0 * qm - pm * q0
            total = lc_P + lc_Q + residual
            assert total == combined, f"Residual decomposition: {lc_P}+{lc_Q}+{residual} != {combined}"

            if residual < 0:
                stats["residual_neg"] += 1
            if stats["min_residual"] is None or residual < stats["min_residual"]:
                stats["min_residual"] = residual
            if lc_P + lc_Q + residual < 0:
                stats["lc_P_plus_lc_Q_plus_residual_neg"] += 1
            if residual < 0:
                if lc_P >= -residual:
                    stats["lc_P_ge_neg_residual"] += 1
                else:
                    stats["lc_P_lt_neg_residual"] += 1

            # Cross-product sign
            cross_pq = p1 * q0 - p0 * q1
            if cross_pq < 0:
                stats["cross_pq_neg"] += 1
            elif cross_pq > 0:
                stats["cross_pq_pos"] += 1

            # Lower mismatch
            if m - 3 >= 0:
                p_m3 = p_poly[m - 3] if m - 3 < len(p_poly) else 0
                lower_mismatch = p_m3 * b1 - p0 * b0
                if lower_mismatch < 0:
                    stats["lower_mismatch_neg"] += 1
                else:
                    stats["lower_mismatch_pos"] += 1

            # Growth comparison for Approach 13
            if m - 3 >= 0:
                p_m3 = p_poly[m - 3] if m - 3 < len(p_poly) else 0
                # Need: p0*b0 <= p_m3*b1 for the approach to work
                if b0 > 0 and p_m3 > 0:
                    if p0 * b0 > p_m3 * b1:
                        stats["growth_comparison_fail"] += 1

        proc.wait()
        print(
            f"n={n:2d}: checked={n_checked_this:8d}",
            flush=True,
        )

    stats["wall_s"] = time.time() - t0

    # Print summary
    print("\n" + "=" * 80)
    print("STRONG C2 RISE COMPENSATION: STRUCTURAL ANALYSIS")
    print("=" * 80)
    print(f"Trees checked: {stats['checked']:,}")
    print(f"Identity verified: {stats['identity_fail']} failures")
    print(f"Neg mismatch identity: {stats['neg_mismatch_identity_fail']} failures")
    print()
    print(f"Mismatch negative: {stats['mismatch_neg']}")
    print(f"Rise failures: {stats['rise_fail']}")
    print(f"Combined negative: {stats['combined_neg']}")
    print()
    print("--- Mode(P) distribution ---")
    print(f"mode(P) >= m-1: {stats['mode_P_ge_m1']}")
    print(f"mode(P) = m-1: {stats['mode_P_eq_m1']}")
    print(f"mode(P) = m-2: {stats['mode_P_eq_m2']}")
    print(f"mode(P) >= m-2: {stats['mode_P_ge_m2']}")
    print(f"mode(P) < m-2: {stats['mode_P_lt_m2']}")
    print(f"** mismatch neg AND mode(P) <= m-2: {stats['mismatch_neg_and_mode_P_le_m2']} **")
    print()
    print("--- Ratio p/b ---")
    print(f"p/b increasing at m-1: {stats['ratio_p_over_b_increasing_at_m1']}")
    print(f"p/b decreasing at m-1: {stats['ratio_p_over_b_decreasing_at_m1']} (= mismatch neg)")
    print()
    print("--- LC surplus bounds ---")
    print(f"lc_surplus >= delta_b * b1: {stats['lc_ge_deltab_times_bm1']} / fail: {stats['lc_ge_deltab_times_bm1_fail']}")
    print(f"lc_surplus >= b0 * delta_b (quadratic) fail: {stats['lc_dominates_quadratic_fail']}")
    print()
    print("--- A log-concavity ---")
    print(f"A LC at m: {stats['A_lc_at_m']} / not: {stats['A_not_lc_at_m']}")
    print(f"A LC at m-1: {stats['A_lc_at_m1']} / not: {stats['A_not_lc_at_m1']}")
    print()
    print("--- LC decomposition ---")
    print(f"Cross term negative: {stats['cross_neg']}")
    print(f"Cross PQ sign: neg={stats['cross_pq_neg']} pos={stats['cross_pq_pos']}")
    print(f"Residual negative: {stats['residual_neg']}")
    print(f"Min residual: {stats['min_residual']}")
    print(f"lc_P + lc_Q + residual < 0: {stats['lc_P_plus_lc_Q_plus_residual_neg']}")
    print(f"lc_P >= |neg residual|: {stats['lc_P_ge_neg_residual']} / <: {stats['lc_P_lt_neg_residual']}")
    print()
    print(f"Lower mismatch neg: {stats['lower_mismatch_neg']} pos: {stats['lower_mismatch_pos']}")
    print(f"Growth comparison fail (approach 13): {stats['growth_comparison_fail']}")
    print(f"Time: {stats['wall_s']:.1f}s")

    # Write output
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
