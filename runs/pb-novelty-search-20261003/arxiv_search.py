"""arXiv novelty and delta sweep for the Poisson-binomial first-descent paper (2026-10-03).
General queries over all dates, then a delta filter for submissions since 2026-07-01."""
import json, re, time, urllib.parse, urllib.request
GENERAL = [
    'all:"Poisson binomial" AND all:"log-concave"',
    'all:"Poisson binomial" AND all:Turan',
    'all:"Poisson-binomial" AND all:"log-concavity"',
    'all:"sums of independent Bernoulli" AND all:"log-concave"',
    'all:"Bernoulli sums" AND all:"log-concave"',
    'all:"real-rooted" AND all:Turan AND all:variance',
    'all:"Newton inequalities" AND all:variance',
    'all:"Newton\'s inequalities" AND all:improvement',
    'all:"ultra log-concave" AND all:variance AND all:mode',
    'all:"strongly log-concave" AND all:Bernoulli',
    'all:"discrete curvature" AND all:"log-concave"',
    'all:"c-log-concave"',
    'all:"Shepp-Olkin"',
]
DELTA = [
    'all:"Poisson binomial"',
    'all:"Poisson-binomial"',
    'all:"log-concave" AND all:Bernoulli',
    'all:"ultra-log-concave"',
    'all:"ultra log-concave"',
    'all:"real-rooted" AND all:probability',
    'all:Turan AND all:"log-concave"',
]
def run(q, n=100):
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode(
        {"search_query": q, "start": 0, "max_results": n, "sortBy": "submittedDate", "sortOrder": "descending"})
    for attempt in range(4):
        try:
            x = urllib.request.urlopen(url, timeout=90).read().decode()
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(10)
    out = []
    for e in re.findall(r"<entry>(.*?)</entry>", x, re.S):
        g = lambda tag: " ".join(re.search(rf"<{tag}>(.*?)</{tag}>", e, re.S).group(1).split())
        out.append({"id": g("id"), "published": g("published")[:10], "title": g("title"),
                    "summary": g("summary"), "authors": re.findall(r"<name>(.*?)</name>", e)})
    return out
res = {"general": {}, "delta": {}}
for q in GENERAL:
    r = run(q); res["general"][q] = r; print(f"[{len(r):3d}] {q}"); time.sleep(3)
for q in DELTA:
    r = [e for e in run(q) if e["published"] >= "2026-07-01"]
    res["delta"][q] = r; print(f"[{len(r):3d} since 2026-07-01] {q}"); time.sleep(3)
json.dump(res, open("arxiv_results.json", "w"), indent=1)
