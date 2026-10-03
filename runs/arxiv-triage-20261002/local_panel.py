"""Second-family panel: qwen3.8:27b (local Ollama) classifies each paper from an excerpt.

One call per paper, HTTP API, think off, JSON format, temperature 0, num_ctx 8192.
Output: JSON lines {id, verdict, reason} to the path given as argv[2].
"""
import json
import re
import sys
import urllib.request
from pathlib import Path

root = Path(__file__).resolve().parent
ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]
out_path = Path(sys.argv[2])
brief = (root / "BRIEF.md").read_text()
head = """Classify an arXiv paper for relevance to Erdos Problem #993: the independence
sequence i_k (number of independent sets of size k) of every tree/forest is unimodal.
Known: true for forests with >= N0 vertices (Fang et al. 2026, via a hard-core CLT);
log-concavity of i_k fails for some trees (order 26); tree independence polynomials
are not real-rooted in general. Live lanes: central-window log-concavity, defect /
matching-bag decompositions, Lorentzian methods, zeros of independence polynomials,
explicit tree families, TP2 and Turan-type margins.
Categories:
- relevant: about independence polynomials / independent-set counts / hard-core model
  of trees, forests or graph classes containing trees; or a general unimodality /
  log-concavity tool whose hypotheses tree independence sequences could meet.
- relevant?: unsure between relevant and peripheral.
- peripheral: log-concavity/unimodality/real-rootedness/Lorentzian results about other
  objects (matroids, Ehrhart, posets, permutations), no direct transfer.
- coincidence: keyword overlap only (log-concave densities, unimodal distributions, etc.).
"""
task = """
## Your task

Below is an excerpt (title, abstract, start of introduction) from one paper.
Classify it into one category from the brief using only the excerpt. Reply
with one JSON object:
{"verdict": "relevant" | "relevant?" | "peripheral" | "coincidence", "reason": "<one sentence>"}
"""
done = set()
if out_path.exists():
    done = {json.loads(l)["id"] for l in out_path.read_text().splitlines() if l.strip()}
with out_path.open("a") as fh:
    for pid in ids:
        if pid in done:
            continue
        body = (root / "papers" / f"{pid}.txt").read_text(errors="replace")
        body = re.sub(r"[ \t]+", " ", body)
        body = re.sub(r"\n{2,}", "\n", body)[:3000]
        prompt = head + task + f"\n===== PAPER {pid} =====\n{body}\n"
        req = urllib.request.Request(
            "http://localhost:11434/api/generate",
            data=json.dumps({
                "model": "qwen3.8:27b", "prompt": prompt, "stream": False,
                "think": False, "format": "json",
                "options": {"temperature": 0, "num_ctx": 8192},
            }).encode(),
            headers={"Content-Type": "application/json"},
        )
        try:
            resp = json.loads(urllib.request.urlopen(req, timeout=300).read())
            ans = json.loads(resp["response"])
            rec = {"id": pid, "verdict": ans.get("verdict"), "reason": ans.get("reason")}
        except Exception as exc:  # record the failure, keep going
            rec = {"id": pid, "verdict": None, "reason": f"ERROR: {exc}"}
        fh.write(json.dumps(rec) + "\n")
        fh.flush()
        print(pid, rec["verdict"], flush=True)
