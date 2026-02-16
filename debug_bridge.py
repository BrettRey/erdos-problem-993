
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scripts.bridge_stitch_scan import get_rooted_signatures, solve_bridge_poly
from indpoly import log_concavity_ratio

print("Checking N=1 signatures...")
s1 = get_rooted_signatures(1)
print(f"Count: {len(s1)}. Example: {s1[0]}")

print("Checking N=3 signatures...")
s3 = get_rooted_signatures(3)
print(f"Count: {len(s3)}. Example: {s3[0]}")

print("Checking N=15 signatures...")
# Just check first few to see if they are valid
try:
    s15 = get_rooted_signatures(15)
    print(f"Count: {len(s15)}")
    
    small_sig_count = 0
    for i, sig in enumerate(s15):
        if len(sig[0]) < 5 or len(sig[1]) < 5:
            small_sig_count += 1
            if i < 5: print(f"Found small signature at {i}: {sig}")
            
    print(f"Small signatures in N=15: {small_sig_count}")
    
except Exception as e:
    print(f"Error generating N=15: {e}")

# Check specific poly 
bad_poly = [1, 3, 11, 36, 59, 46, 18, 3]
lc, k = log_concavity_ratio(bad_poly)
print(f"Bad poly LC ratio: {lc} at k={k}")

# Try to reproduce 1.22
# If it came from N1=1, N2=15?
if 's15' in locals() and s15:
    print("Testing pairs N1=1 x N2=15...")
    max_lc = 0
    for sig1 in s1:
        for sig2 in s15:
            p = solve_bridge_poly(sig1, sig2)
            lc, _ = log_concavity_ratio(p)
            if lc > 1.2:
                print(f"Found High LC: {lc}")
                print(f"Poly: {p}")
                print(f"Sig1: {sig1}")
                print(f"Sig2: {sig2}")
                break
        if max_lc > 1.2: break
