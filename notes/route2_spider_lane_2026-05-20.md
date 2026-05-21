# Route-2 Spider Lane Scan (2026-05-20)

## Purpose

The all-degree-2-support endpoint scan found its minimum endpoint margin at the
spider `S(1,2^7,4)`. This suggested that the old mixed-spider lane
`S(2^k,1^j)` was too narrow. I added a structured scanner for
`S(1^j,2^a,r)`.

Script:

```bash
python3 gpt_attack/route2_spider_lane_scan.py
```

It uses the spider product formula

```text
I_T(x) = prod_L P_L(x) + x prod_L P_{L-1}(x),
```

where `P_L` is the independence polynomial of a path on `L` vertices.

## Scans

Smoke/exact-compatible scan:

```bash
python3 gpt_attack/route2_spider_lane_scan.py \
  --a-max 20 --r-max 12 --j-max 1 \
  --out /tmp/spider_lane_smoke.json
```

Output:

- checked `819` records,
- minimum endpoint margin: `0.164846737171` at `S(2^12,8)`, removing a length-2 arm,
- minimum route-2 slack: `0.180874563128` at `S(2^18,8)`, removing a length-2 arm.

Wider float scan:

```bash
python3 gpt_attack/route2_spider_lane_scan.py \
  --float --a-max 200 --r-max 80 --j-max 1 \
  --out results/route2_spider_lane_j0_1_a200_r80_float.json
```

Output:

- checked `62,643` records,
- minimum endpoint margin: `0.164846737171` at `S(2^12,8)`, removing a length-2 arm,
- minimum route-2 slack: `0.168459434752` at `S(2^198,8)`, removing a length-2 arm,
- minimum stronger-threshold slack: `0.168127340824` at `S(2^199)`, removing a length-2 arm,
- maximum `tau`-deficit: `0.163368617465` at `S(2^198,8)`, removing a length-2 arm.

Wider unsaved float scan:

```bash
python3 gpt_attack/route2_spider_lane_scan.py \
  --float --a-max 500 --r-max 80 --j-max 1
```

Output:

- checked `156,843` records,
- minimum endpoint margin still `0.164846737171` at `S(2^12,8)`,
- minimum route-2 slack: `0.167392925993` at `S(2^498,8)`,
- minimum stronger-threshold slack: `0.167250374813` at `S(2^499)`,
- maximum `tau`-deficit: `0.165338645418` at `S(2^500,4)`.

Fixed-`r` asymptotic probe, after adding log-scaled mean evaluation:

```bash
python3 - <<'PY'
from gpt_attack.route2_spider_lane_scan import path_polys, route2_float_for_removed_arm
paths=path_polys(100)
for r in [2,3,4,5,6,7,8,9,10,12,16,20,40,80]:
    best=None
    for a in range(3,1001):
        counts={2:a}
        counts[r]=counts.get(r,0)+1
        rec=route2_float_for_removed_arm(counts,2,paths)
        if rec and rec['lam_margin'] is not None:
            if best is None or rec['lam_margin'] < best[0]:
                best=(rec['lam_margin'],a,rec)
    print(r, best)
PY
```

Summary:

- For all sampled fixed `r`, route-2 slack remains positive and approaches the
  same `1/6` asymptote from above.
- The smallest endpoint margin in this fixed-`r` probe remains the finite
  witness `S(2^12,8)`, with margin `0.164846737171`.
- Representative large-`a` values:
  - `r=2`, `a=1000`: endpoint margin `0.168329343973`, route-2 slack
    `0.168662674650`;
  - `r=4`, `a=1000`: endpoint margin `0.167613926779`, route-2 slack
    `0.167946462863`;
  - `r=8`, `a=1000`: endpoint margin `0.168024666079`, route-2 slack
    `0.168356662360`;
  - `r=80`, `a=1000`: endpoint margin `0.169662526627`, route-2 slack
    `0.169984214548`.

The `r=2` line is off by one relative to the pure-spider notation in
`notes/route2_pure_spider_asymptotic_2026-05-20.md`, because this probe
constructs `S(2^a,r)` and then adds the `r=2` arm to the same count.

## Exact Candidate Verification

Exact rational verification was run with:

```bash
python3 - <<'PY'
import sys
sys.set_int_max_str_digits(1000000)
from gpt_attack.route2_spider_lane_scan import path_polys, route2_for_removed_arm, fmt_counts
candidates=[({2:12,8:1},2),({2:198,8:1},2),({2:199},2),({1:1,2:200,80:1},2)]
paths=path_polys(100)
for counts, rem in candidates:
    rec=route2_for_removed_arm(counts,rem,paths)
    print(fmt_counts(counts), rem, float(rec["lam_margin"] or 0), float(rec["route2_slack"]))
PY
```

Key exact values:

- `S(2^12,8)`, remove length-2 arm:
  - `lambda = 5978621/5989843`,
  - `tau = 175705/193406`,
  - endpoint margin `M_lam = 0.164846737171075`,
  - route-2 slack `0.182972176581926`.
- `S(2^198,8)`, remove length-2 arm:
  - endpoint margin `M_lam = 0.166825739544753`,
  - route-2 slack `0.168459434752071`.
- `S(2^199)`, remove length-2 arm:
  - endpoint margin `M_lam = 0.168333417085427`,
  - route-2 slack `0.17`,
  - stronger-threshold slack `0.16812734082397`.

## Interpretation

The endpoint margin has a finite dip below `1/6` at `S(2^12,8)` and then
appears to climb back toward `1/6` along the long `S(2^a,8)` lane.

The route-2 and stronger-threshold slacks decrease toward `1/6` from above in
the scanned lanes. The maximum `tau`-deficit also appears to approach `1/6`,
with `S(2^a,4)` becoming the large-`a` extremal for that diagnostic.

## Next Symbolic Target

The next useful proof target is no longer just mixed spiders `S(2^k,1^j)`.
It should cover at least:

```text
S(2^a,r), removing a length-2 arm, especially r = 8;
S(2^a,4), removing the length-4 arm, for maximum tau-deficit;
S(1,2^a,r), as a check that the extra unit arm does not reduce the endpoint margin.
```

A plausible finite/asymptotic strategy:

1. Derive closed forms for the Route-2 endpoint margin on `S(2^a,r)` for fixed
   `r`.
2. Prove positivity for all sufficiently large `a`, with limiting margin
   `>= 1/6`.
3. Verify the finite window exactly for small `a` and moderate `r`.
