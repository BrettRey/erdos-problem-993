"""K6: attack R1 with a centre adjacent to h hub-stars (each hub: m cherries),
so the centre's closed neighbourhood collects h positive hub defects.
Exact checks via k5_repair_families.check. Output JSON lines."""
import json
import sys
from k3_families import B
from k5_repair_families import check


def star_of_hubs(h, m, t=2, centre_leaves=0):
    b = B()
    for _ in range(h):
        hub = b.add(0)
        for _ in range(m):
            c = b.add(hub)
            for _ in range(t):
                b.add(c)
    for _ in range(centre_leaves):
        b.add(0)
    return b.adj()


if __name__ == "__main__":
    for h in [2, 3, 4, 5]:
        for m in [9, 11, 13, 15]:
            if 1 + h * (1 + 3 * m) > 150:
                continue
            r = check(star_of_hubs(h, m))
            r["family"] = f"starofhubs_h{h}_m{m}"
            print(json.dumps(r), flush=True)
