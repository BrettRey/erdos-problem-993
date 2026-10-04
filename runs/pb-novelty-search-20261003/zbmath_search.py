"""zbMATH Open novelty sweep for the Poisson-binomial first-descent paper (2026-10-03).
Runs each query against api.zbmath.org, saves raw JSON, prints deduplicated hits."""
import json, time, urllib.error, urllib.parse, urllib.request
QUERIES = [
    "Poisson binomial log-concave",
    "Poisson binomial Turan",
    "Poisson binomial variance mode",
    "Poisson binomial ultra log-concave",
    "Bernoulli sums log-concave",
    "sums of independent Bernoulli log-concavity",
    "independent Bernoulli ultra log-concave",
    "strongly log-concave discrete Bernoulli",
    "Newton inequalities improvement real-rooted",
    "Newton inequalities refinement coefficients",
    "real-rooted polynomial coefficients log-concavity variance",
    "Turan inequalities real-rooted polynomial coefficients",
    "Shepp-Olkin conjecture",
    "discrete log-Sobolev Bakry-Emery log-concave",
    "c-log-concave discrete curvature",
    "Poisson binomial mode ratio bound",
    "log-concave sequence quantitative curvature variance",
]
hits = {}
raw = {}
for q in QUERIES:
    url = "https://api.zbmath.org/v1/document/_search?" + urllib.parse.urlencode(
        {"search_string": q, "page": 0, "results_per_page": 100})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                d = json.load(r)
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                d = {"status": {"nr_total_results": 0}, "result": []}  # zbMATH answers 404 for no hits
                break
            if e.code < 500 or attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))
    raw[q] = d
    n = d["status"]["nr_total_results"]
    print(f"[{n:4d} total] {q}")
    for x in d.get("result", []):
        t = x.get("title", {})
        title = t.get("title") if isinstance(t, dict) else t
        authors = [a.get("name") for a in x.get("contributors", {}).get("authors", [])]
        key = x.get("id")
        hits.setdefault(key, {"year": x.get("year"), "title": title, "authors": authors, "queries": []})["queries"].append(q)
    time.sleep(1)
json.dump(raw, open("zbmath_raw.json", "w"))
json.dump(hits, open("zbmath_hits.json", "w"), indent=1)
print(len(hits), "distinct documents")
