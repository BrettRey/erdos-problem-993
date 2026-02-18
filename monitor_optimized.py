import time
import re
import sys
import os

log_file = "out_radio_optimized/search.log"
best_nm = 0.0
start_time = time.time()
duration = 120  # Monitor for 2 minutes

print(f"Monitoring {log_file} for {duration} seconds...")
print("-" * 60)

# Track file position
try:
    with open(log_file, "r") as f:
        f.seek(0, os.SEEK_END)
        pos = f.tell()
except FileNotFoundError:
    pos = 0

while time.time() - start_time < duration:
    if os.path.exists(log_file):
        with open(log_file, "r") as f:
            f.seek(pos)
            lines = f.readlines()
            pos = f.tell()
            
            for line in lines:
                line = line.strip()
                if "COUNTEREXAMPLE FOUND" in line:
                    print("\n" + "!"*60)
                    print(f"🚨 ALERT: {line}")
                    print("!"*60 + "\n")
                    sys.exit(0)
                
                # [Iter 3] Best: ...
                if "[Iter" in line and "Best:" in line:
                    print(line)
                    
                # Score=18.82 (NM=0.9704, Re=0.4107, Im=5.2780)
                if "Score=" in line:
                    try:
                        # Extract metrics
                        nm_match = re.search(r"NM=([\d\.]+)", line)
                        if nm_match:
                            nm = float(nm_match.group(1))
                            if nm > best_nm:
                                best_nm = nm
                                print(f"  >>> NEW BEST NEAR-MISS RATIO: {nm:.5f}")
                            if nm > 0.99:
                                print(f"  🔥 CRITICAL: Approaching counterexample threshold! (NM={nm})")
                    except Exception:
                        pass
                    print(f"  {line}")

    time.sleep(2)

print("-" * 60)
print(f"Monitoring complete. Highest NM seen this session: {best_nm:.5f}")
