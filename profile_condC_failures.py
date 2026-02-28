"""Profile Strong Condition C product-closure failures when J <= E is absent.

Generates random synthetic pairs (I, E) with:
  - I = E + x*J, so I >= E coefficientwise and I[0] = E[0]
  - E[0] = 1
  - E is log-concave
  - J has nonneg coefficients
  - Strong Condition C holds for (I, E)

For pairs where J > E at some coefficient (the "relaxed" case), multiplies
two such pairs and checks if Condition C holds for the product.

Reports detailed statistics on failures vs passes, thresholds, and correlations.
"""

import random
import math
from collections import defaultdict


# ── Polynomial helpers ──────────────────────────────────────────────────

def polymul(a, b):
    """Pure-Python polynomial multiplication."""
    if not a or not b:
        return []
    la, lb = len(a), len(b)
    out = [0] * (la + lb - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[i + j] += ai * bj
    return out


def polyadd(a, b):
    la, lb = len(a), len(b)
    out = [0] * max(la, lb)
    for i in range(la):
        out[i] += a[i]
    for i in range(lb):
        out[i] += b[i]
    return out


def coeff(poly, k):
    return poly[k] if 0 <= k < len(poly) else 0


def is_log_concave(seq):
    for k in range(1, len(seq) - 1):
        if seq[k] * seq[k] < seq[k - 1] * seq[k + 1]:
            return False
    return True


# ── Condition C check (returns value and details) ──────────────────────

def check_condC(I_poly, E_poly):
    """Check Condition C for (I, E). Returns True if passes."""
    d_seq = []
    L = max(len(I_poly), len(E_poly)) + 1
    for k in range(L):
        dk = coeff(I_poly, k + 1) * coeff(E_poly, k) - coeff(I_poly, k) * coeff(E_poly, k + 1)
        d_seq.append(dk)
    for k in range(1, L):
        ekm1 = coeff(E_poly, k - 1)
        if ekm1 == 0:
            continue
        ek = coeff(E_poly, k)
        ekp1 = coeff(E_poly, k + 1)
        ikm1 = coeff(I_poly, k - 1)
        ck = ek * ek - ekm1 * ekp1
        val = ekm1 * d_seq[k] + ek * d_seq[k - 1] + ikm1 * ck
        if val < 0:
            return False
    return True


def condC_details(I_poly, E_poly):
    """Return (passes, worst_val, worst_k, min_margin_val, min_margin_k).

    worst_val/k: most negative Condition C value (only meaningful if fails).
    min_margin_val/k: smallest nonneg value (closest to failure).
    """
    d_seq = []
    L = max(len(I_poly), len(E_poly)) + 1
    for k in range(L):
        dk = coeff(I_poly, k + 1) * coeff(E_poly, k) - coeff(I_poly, k) * coeff(E_poly, k + 1)
        d_seq.append(dk)

    worst_val = float('inf')
    worst_k = -1
    min_margin_val = float('inf')
    min_margin_k = -1
    passes = True

    for k in range(1, L):
        ekm1 = coeff(E_poly, k - 1)
        if ekm1 == 0:
            continue
        ek = coeff(E_poly, k)
        ekp1 = coeff(E_poly, k + 1)
        ikm1 = coeff(I_poly, k - 1)
        ck = ek * ek - ekm1 * ekp1
        val = ekm1 * d_seq[k] + ek * d_seq[k - 1] + ikm1 * ck
        if val < worst_val:
            worst_val = val
            worst_k = k
        if val >= 0 and val < min_margin_val:
            min_margin_val = val
            min_margin_k = k

    if worst_val < 0:
        passes = False

    return passes, worst_val, worst_k, min_margin_val, min_margin_k


# ── Random LC polynomial generator ────────────────────────────────────

def random_lc_poly(deg, max_coeff=10):
    """Generate a random UNIMODAL log-concave polynomial with nonneg integer coeffs.

    Strategy: start with coeff[0] = 1, build greedily ensuring LC at each step.
    Once a coefficient is 0, all subsequent must be 0 (LC + nonneg).
    """
    if deg == 0:
        return [1]

    coeffs = [1]
    hit_zero = False
    for k in range(1, deg + 1):
        if hit_zero:
            coeffs.append(0)
            continue

        # LC constraint from previous: coeffs[k-1]^2 >= coeffs[k-2]*coeffs[k]
        # => coeffs[k] <= coeffs[k-1]^2 / coeffs[k-2] if k >= 2 and coeffs[k-2] > 0
        if k >= 2 and coeffs[k - 2] > 0:
            upper = coeffs[k - 1] ** 2 // coeffs[k - 2]
        elif k >= 2 and coeffs[k - 2] == 0:
            # If coeffs[k-2] = 0, then coeffs[k-1] must also be 0 for LC
            # (since coeffs[k-1]^2 >= coeffs[k-2]*coeffs[k] = 0 is trivially true,
            #  but coeffs[k-2]^2 >= coeffs[k-3]*coeffs[k-1] requires coeffs[k-1]=0 if coeffs[k-2]=0)
            upper = 0
        else:
            upper = max_coeff

        upper = min(upper, max_coeff)
        if upper <= 0:
            coeffs.append(0)
            hit_zero = True
        else:
            val = random.randint(0, upper)
            coeffs.append(val)
            if val == 0:
                hit_zero = True

    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()

    return coeffs


def random_nonneg_poly(deg, max_coeff=10):
    """Random polynomial with nonneg integer coefficients, coeff[0] >= 0."""
    coeffs = [random.randint(0, max_coeff) for _ in range(deg + 1)]
    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


# ── J > E metrics ──────────────────────────────────────────────────────

def j_exceeds_e_metrics(J, E):
    """Compute metrics for how much J exceeds E.

    Returns dict with:
      has_excess: bool (True if J_k > E_k for some k)
      max_excess: max(J_k - E_k) over k where J_k > E_k
      max_E: max(E_k) over all k
      rel_excess: max_excess / max_E if max_E > 0
      num_excess_positions: how many k have J_k > E_k
      excess_positions: list of k where J_k > E_k
    """
    max_excess = 0
    excess_positions = []
    max_E = max(E) if E else 0

    L = max(len(J), len(E))
    for k in range(L):
        jk = coeff(J, k)
        ek = coeff(E, k)
        if jk > ek:
            diff = jk - ek
            if diff > max_excess:
                max_excess = diff
            excess_positions.append(k)

    has_excess = len(excess_positions) > 0
    rel_excess = max_excess / max_E if max_E > 0 else 0.0

    return {
        'has_excess': has_excess,
        'max_excess': max_excess,
        'max_E': max_E,
        'rel_excess': rel_excess,
        'num_excess_positions': len(excess_positions),
        'excess_positions': excess_positions,
    }


# ── Main profiling ─────────────────────────────────────────────────────

def gen_pair_constrained(deg_E, deg_J, max_coeff, force_je=False):
    """Generate a single valid (I, E, J, metrics) tuple or None.

    If force_je=True, generate J with J_k <= E_k for all k.
    """
    E = random_lc_poly(deg_E, max_coeff)
    if not is_log_concave(E):
        return None, 'lc'

    if force_je:
        # Generate J with J_k <= E_{k+1} for each k (since I[k+1] = E[k+1] + J[k])
        # Actually J_k should be <= E_k for the J<=E test.
        J = []
        for k in range(deg_J + 1):
            ek = coeff(E, k)
            if ek <= 0:
                J.append(0)
            else:
                J.append(random.randint(0, ek))
    else:
        J = random_nonneg_poly(deg_J, max_coeff)

    # I = E + x*J
    xJ = [0] + list(J)
    I = polyadd(list(E), xJ)
    while len(I) > 1 and I[-1] == 0:
        I.pop()

    if not check_condC(I, E):
        return None, 'condC'

    # Check I is unimodal (tree IS polynomials should be unimodal)
    decreasing = False
    for idx in range(1, len(I)):
        if I[idx] > I[idx - 1]:
            if decreasing:
                return None, 'not_unimodal'
        elif I[idx] < I[idx - 1]:
            decreasing = True

    metrics = j_exceeds_e_metrics(J, E)

    if len(E) <= 1 or len(I) <= 1:
        return None, 'trivial'

    # If we forced J<=E but it still has excess (shouldn't happen), reject
    if force_je and metrics['has_excess']:
        return None, 'force_fail'

    return (I, E, J, metrics), 'ok'


def main():
    random.seed(42)

    N_JE_PAIRS = 30_000       # J <= E pairs (harder to generate)
    N_JGT_PAIRS = 170_000     # J > E pairs (easy to generate)
    DEG_RANGE = (3, 8)
    MAX_COEFF = 10

    # Phase 1: Generate valid (I, E) pairs satisfying all constraints
    print("Phase 1: Generating valid (I, E) pairs...")
    print(f"  Target: {N_JE_PAIRS:,d} J<=E pairs + {N_JGT_PAIRS:,d} J>E pairs")
    print(f"  Degree {DEG_RANGE[0]}-{DEG_RANGE[1]}, coeffs 0-{MAX_COEFF}")

    pairs_je = []    # J <= E
    pairs_jgt = []   # J > E
    attempts = 0
    rejected = {'lc': 0, 'condC': 0, 'trivial': 0, 'force_fail': 0, 'not_unimodal': 0}

    # Generate J <= E pairs
    while len(pairs_je) < N_JE_PAIRS:
        attempts += 1
        deg_E = random.randint(DEG_RANGE[0], DEG_RANGE[1])
        deg_J = random.randint(DEG_RANGE[0], min(deg_E, DEG_RANGE[1]))
        result, reason = gen_pair_constrained(deg_E, deg_J, MAX_COEFF, force_je=True)
        if result is None:
            rejected[reason] += 1
            continue
        pairs_je.append(result)
        if len(pairs_je) % 10000 == 0:
            print(f"    J<=E: {len(pairs_je):,d} collected ({attempts:,d} attempts)")

    print(f"  J<=E pairs: {len(pairs_je):,d} in {attempts:,d} attempts")

    # Generate J > E pairs (unconstrained J, filter to those with excess)
    attempts2 = 0
    while len(pairs_jgt) < N_JGT_PAIRS:
        attempts2 += 1
        deg_E = random.randint(DEG_RANGE[0], DEG_RANGE[1])
        deg_J = random.randint(DEG_RANGE[0], DEG_RANGE[1])
        result, reason = gen_pair_constrained(deg_E, deg_J, MAX_COEFF, force_je=False)
        if result is None:
            rejected[reason] += 1
            continue
        # Keep only if J > E at some coeff
        if not result[3]['has_excess']:
            # It's a J<=E pair, add to je pool if still needed
            continue
        pairs_jgt.append(result)
        if len(pairs_jgt) % 50000 == 0:
            print(f"    J>E:  {len(pairs_jgt):,d} collected ({attempts2:,d} attempts)")

    print(f"  J>E pairs:  {len(pairs_jgt):,d} in {attempts2:,d} attempts")

    pairs = pairs_je + pairs_jgt
    total_attempts = attempts + attempts2
    print(f"\n  Total pairs: {len(pairs):,d} in {total_attempts:,d} attempts")
    print(f"  Rejected: {rejected}")
    print(f"  J<=E: {len(pairs_je):,d}, J>E: {len(pairs_jgt):,d}")

    # ── Phase 2: Product closure test ──────────────────────────────────

    # Ensure balanced sampling: 50K both_je, 75K one_jgt, 75K both_jgt
    N_BOTH_JE = 50_000
    N_ONE_JGT = 75_000
    N_BOTH_JGT = 75_000
    N_PRODUCTS = N_BOTH_JE + N_ONE_JGT + N_BOTH_JGT
    print(f"\nPhase 2: Testing {N_PRODUCTS:,d} products for Condition C closure...")
    print(f"  both_je: {N_BOTH_JE:,d}, one_jgt: {N_ONE_JGT:,d}, both_jgt: {N_BOTH_JGT:,d}")

    # Categories:
    #   both_je: both factors have J <= E
    #   one_jgt: exactly one factor has J > E
    #   both_jgt: both factors have J > E

    results = {
        'both_je': {'total': 0, 'fails': 0, 'fail_details': [],
                     'min_margin': float('inf'),
                     'min_margin_k': -1, 'min_margin_pair': None},
        'one_jgt': {'total': 0, 'fails': 0, 'fail_details': []},
        'both_jgt': {'total': 0, 'fails': 0, 'fail_details': []},
    }

    # Track failure details
    fail_by_deg = defaultdict(int)
    fail_by_k = defaultdict(int)
    fail_magnitudes = []
    fail_rel_excesses = []
    pass_jgt_min_margins = []
    pass_jgt_rel_excesses = []

    # Build trial list with controlled categories
    trial_plan = []
    trial_plan += [('both_je', False, False)] * N_BOTH_JE
    trial_plan += [('one_jgt', True, False)] * N_ONE_JGT
    trial_plan += [('both_jgt', True, True)] * N_BOTH_JGT
    random.shuffle(trial_plan)

    for trial_idx, (cat, f1_jgt, f2_jgt) in enumerate(trial_plan):
        # Pick factors from appropriate pools
        if f1_jgt:
            I1, E1, J1, m1 = pairs_jgt[random.randint(0, len(pairs_jgt) - 1)]
        else:
            I1, E1, J1, m1 = pairs_je[random.randint(0, len(pairs_je) - 1)]
        if f2_jgt:
            I2, E2, J2, m2 = pairs_jgt[random.randint(0, len(pairs_jgt) - 1)]
        else:
            I2, E2, J2, m2 = pairs_je[random.randint(0, len(pairs_je) - 1)]

        cat1_jgt = m1['has_excess']
        cat2_jgt = m2['has_excess']

        # Form product
        A = polymul(I1, I2)
        B = polymul(E1, E2)

        passes, worst_val, worst_k, min_margin_val, min_margin_k = condC_details(A, B)

        results[cat]['total'] += 1

        if not passes:
            results[cat]['fails'] += 1

            prod_deg = len(A) - 1
            fail_by_deg[prod_deg] += 1
            fail_by_k[worst_k] += 1
            fail_magnitudes.append(worst_val)

            # Record the max J > E relative excess for the factors
            re1 = m1['rel_excess'] if cat1_jgt else 0.0
            re2 = m2['rel_excess'] if cat2_jgt else 0.0
            max_re = max(re1, re2)
            fail_rel_excesses.append(max_re)

            if len(results[cat]['fail_details']) < 10:
                results[cat]['fail_details'].append({
                    'I1': I1, 'E1': E1, 'J1': J1,
                    'I2': I2, 'E2': E2, 'J2': J2,
                    'm1': m1, 'm2': m2,
                    'worst_val': worst_val, 'worst_k': worst_k,
                    'A': A, 'B': B,
                })
        else:
            # Pass
            if cat != 'both_je':
                # Track margin for J > E passes
                if min_margin_val < float('inf'):
                    pass_jgt_min_margins.append(min_margin_val)
                    re1 = m1['rel_excess'] if cat1_jgt else 0.0
                    re2 = m2['rel_excess'] if cat2_jgt else 0.0
                    pass_jgt_rel_excesses.append(max(re1, re2))
            else:
                if min_margin_val < results['both_je']['min_margin']:
                    results['both_je']['min_margin'] = min_margin_val
                    results['both_je']['min_margin_k'] = min_margin_k
                    results['both_je']['min_margin_pair'] = (I1, E1, I2, E2)

        if (trial_idx + 1) % 50000 == 0:
            total_fails = sum(r['fails'] for r in results.values())
            print(f"    {trial_idx+1:,d} products tested, {total_fails} failures so far")

    # ── Phase 3: Report ────────────────────────────────────────────────

    print("\n" + "=" * 72)
    print("RESULTS")
    print("=" * 72)

    total_tested = sum(r['total'] for r in results.values())
    total_fails = sum(r['fails'] for r in results.values())
    print(f"\nTotal products tested: {total_tested:,d}")
    print(f"Total failures:       {total_fails}")

    for cat, label in [('both_je', 'Both J <= E'),
                        ('one_jgt', 'One J > E'),
                        ('both_jgt', 'Both J > E')]:
        r = results[cat]
        rate = r['fails'] / r['total'] if r['total'] > 0 else 0
        print(f"\n  {label}:")
        print(f"    Tested: {r['total']:,d}")
        print(f"    Fails:  {r['fails']} ({rate:.4%})")

    # ── J <= E closest call ────────────────────────────────────────────

    print("\n" + "-" * 72)
    print("CLOSEST CALL (J <= E pairs, smallest Condition C margin at product)")
    r = results['both_je']
    if r['min_margin'] < float('inf'):
        print(f"  Min margin value: {r['min_margin']}")
        print(f"  At coefficient k: {r['min_margin_k']}")
        if r['min_margin_pair'] is not None:
            I1, E1, I2, E2 = r['min_margin_pair']
            print(f"  Factor 1: I={I1}, E={E1}")
            print(f"  Factor 2: I={I2}, E={E2}")
    else:
        print("  (No valid margin recorded)")

    # ── Failure analysis ───────────────────────────────────────────────

    if total_fails > 0:
        print("\n" + "-" * 72)
        print("FAILURE ANALYSIS")

        print(f"\n  Failures by product degree:")
        for deg in sorted(fail_by_deg):
            print(f"    deg {deg}: {fail_by_deg[deg]}")

        print(f"\n  Failures by worst-k (Condition C violated at coefficient k):")
        for k in sorted(fail_by_k):
            print(f"    k={k}: {fail_by_k[k]}")

        print(f"\n  Violation magnitudes (how negative):")
        fail_magnitudes.sort()
        print(f"    Min (worst): {fail_magnitudes[0]}")
        print(f"    Max (mildest): {fail_magnitudes[-1]}")
        print(f"    Median: {fail_magnitudes[len(fail_magnitudes)//2]}")
        if len(fail_magnitudes) >= 10:
            print(f"    10th percentile: {fail_magnitudes[len(fail_magnitudes)//10]}")
            print(f"    90th percentile: {fail_magnitudes[9*len(fail_magnitudes)//10]}")

        print(f"\n  Relative J-excess of failing pairs (max(J_k-E_k)/max(E_k)):")
        fail_rel_excesses.sort()
        print(f"    Min: {fail_rel_excesses[0]:.4f}")
        print(f"    Max: {fail_rel_excesses[-1]:.4f}")
        print(f"    Median: {fail_rel_excesses[len(fail_rel_excesses)//2]:.4f}")

        # Threshold analysis: bin by relative excess
        print(f"\n  Failure rate by relative excess bin:")
        # Collect ALL products (not just fails) with their rel_excess and outcome
        # We'll need to re-do this more carefully. Instead, let's use what we have.
        # We can compute from fail_rel_excesses and pass_jgt_rel_excesses.
        all_jgt_excesses = [(re, 'fail') for re in fail_rel_excesses] + \
                           [(re, 'pass') for re in pass_jgt_rel_excesses]
        all_jgt_excesses.sort()

        bins = [0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 5.0, float('inf')]
        for lo, hi in zip(bins[:-1], bins[1:]):
            in_bin = [x for x in all_jgt_excesses if lo <= x[0] < hi]
            n_fail = sum(1 for x in in_bin if x[1] == 'fail')
            n_total = len(in_bin)
            if n_total > 0:
                rate = n_fail / n_total
                hi_str = f"{hi:.1f}" if hi < float('inf') else "inf"
                print(f"    [{lo:.1f}, {hi_str}): {n_fail}/{n_total} = {rate:.4%}")

        # ── Sample failures ────────────────────────────────────────────

        print(f"\n  Sample failures (up to 5):")
        for cat in ['both_je', 'one_jgt', 'both_jgt']:
            details = results[cat]['fail_details'][:5]
            if not details:
                continue
            print(f"\n    Category: {cat}")
            for idx, d in enumerate(details):
                print(f"\n    Failure #{idx+1}:")
                print(f"      I1 = {d['I1']}")
                print(f"      E1 = {d['E1']}")
                print(f"      J1 = {d['J1']}")
                print(f"      J1>E1 excess: max_excess={d['m1']['max_excess']}, "
                      f"rel={d['m1']['rel_excess']:.4f}, "
                      f"positions={d['m1']['excess_positions']}")
                print(f"      I2 = {d['I2']}")
                print(f"      E2 = {d['E2']}")
                print(f"      J2 = {d['J2']}")
                print(f"      J2>E2 excess: max_excess={d['m2']['max_excess']}, "
                      f"rel={d['m2']['rel_excess']:.4f}, "
                      f"positions={d['m2']['excess_positions']}")
                print(f"      Product CondC worst value: {d['worst_val']} at k={d['worst_k']}")
                print(f"      Product A = {d['A']}")
                print(f"      Product B = {d['B']}")

    # ── Pass analysis for J > E ────────────────────────────────────────

    if pass_jgt_min_margins:
        print("\n" + "-" * 72)
        print("PASS ANALYSIS (J > E pairs that PASSED)")
        pass_jgt_min_margins.sort()
        print(f"  Count: {len(pass_jgt_min_margins):,d}")
        print(f"  Smallest margin: {pass_jgt_min_margins[0]}")
        print(f"  Median margin:   {pass_jgt_min_margins[len(pass_jgt_min_margins)//2]}")
        if len(pass_jgt_min_margins) >= 10:
            print(f"  10th percentile: {pass_jgt_min_margins[len(pass_jgt_min_margins)//10]}")

    # ── Correlation check ──────────────────────────────────────────────

    if total_fails > 0 and fail_rel_excesses and pass_jgt_rel_excesses:
        print("\n" + "-" * 72)
        print("CORRELATION: Relative excess vs failure")
        fail_mean = sum(fail_rel_excesses) / len(fail_rel_excesses)
        pass_mean = sum(pass_jgt_rel_excesses) / len(pass_jgt_rel_excesses) if pass_jgt_rel_excesses else 0
        print(f"  Mean relative excess (failures):  {fail_mean:.4f}")
        print(f"  Mean relative excess (passes):    {pass_mean:.4f}")
        print(f"  Ratio fail/pass:                  {fail_mean/pass_mean:.4f}" if pass_mean > 0 else "  (no passes)")

    # ── One-factor vs both-factor analysis ─────────────────────────────

    if total_fails > 0:
        one_fails = results['one_jgt']['fails']
        both_fails = results['both_jgt']['fails']
        one_total = results['one_jgt']['total']
        both_total = results['both_jgt']['total']
        print("\n" + "-" * 72)
        print("ONE-FACTOR vs BOTH-FACTOR")
        if one_total > 0:
            print(f"  One J>E:  {one_fails}/{one_total} = {one_fails/one_total:.4%}")
        if both_total > 0:
            print(f"  Both J>E: {both_fails}/{both_total} = {both_fails/both_total:.4%}")

    # ── Degree analysis ────────────────────────────────────────────────

    print("\n" + "-" * 72)
    print("DEGREE DISTRIBUTION OF GENERATED PAIRS")
    deg_dist = defaultdict(int)
    for I, E, J, m in pairs:
        d = len(I) - 1
        deg_dist[d] += 1
    for d in sorted(deg_dist):
        print(f"  deg {d}: {deg_dist[d]:,d}")

    je_excess_dist = defaultdict(int)
    for _, _, _, m in pairs:
        if m['has_excess']:
            # Bin relative excess
            re = m['rel_excess']
            if re < 0.1:
                je_excess_dist['<0.1'] += 1
            elif re < 0.5:
                je_excess_dist['0.1-0.5'] += 1
            elif re < 1.0:
                je_excess_dist['0.5-1.0'] += 1
            elif re < 2.0:
                je_excess_dist['1.0-2.0'] += 1
            else:
                je_excess_dist['>=2.0'] += 1
    if je_excess_dist:
        print(f"\n  Relative excess distribution (J > E pairs only):")
        for bucket in ['<0.1', '0.1-0.5', '0.5-1.0', '1.0-2.0', '>=2.0']:
            if bucket in je_excess_dist:
                print(f"    {bucket}: {je_excess_dist[bucket]:,d}")

    print("\n" + "=" * 72)
    print("DONE")


if __name__ == '__main__':
    main()
