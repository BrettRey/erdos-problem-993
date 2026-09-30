"""Write manifest.yaml: model, git HEAD, SHA-256 of inputs read and of every script/data file here."""
import hashlib, os, subprocess, datetime
root = os.path.dirname(os.path.abspath(__file__))
proj = os.path.abspath(os.path.join(root, '..', '..', '..'))
def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()
head = subprocess.check_output(['git', '-C', proj, 'rev-parse', 'HEAD']).decode().strip()
inputs = ['runs/attack-20260930/route-freecount/RETURN.json', 'runs/attack-20260930/route-freecount/pv_lib.py',
          'runs/attack-20260930/route-freecount/k5_repair_families.py', 'runs/attack-20260930/route-freecount/k6_star_of_hubs.py',
          'runs/attack-20260930/route-freecount/k4_hillclimb.py', 'runs/attack-20260930/route-freecount/k3_families.py',
          'runs/attack-20260930/route-freecount/data/k5_repair_families.jsonl', 'runs/attack-20260930/route-freecount/data/k6_star_of_hubs.jsonl',
          'runs/attack-20260930/census/families.py']
lines = ['# r1-adversarial manifest', f'generated: {datetime.datetime.now().isoformat(timespec="seconds")}',
         'model: claude-opus-5-5 (Claude Code subagent, workflow lane r1-adversarial)',
         'python: venv/bin/python (python-flint 0.9.0) for all exact computations',
         f'git_head: {head}', 'task: adversarial search against lemma R1 (closed-neighbourhood free-count average)',
         'inputs_read:']
for p in inputs:
    ap = os.path.join(proj, p)
    lines.append(f'  - path: {p}\n    sha256: {sha(ap)}')
lines.append('outputs:')
for dp, dn, fn in os.walk(root):
    if '__pycache__' in dp:
        continue
    for f in sorted(fn):
        if f == 'manifest.yaml':
            continue
        ap = os.path.join(dp, f)
        lines.append(f'  - path: {os.path.relpath(ap, proj)}\n    sha256: {sha(ap)}\n    bytes: {os.path.getsize(ap)}')
open(os.path.join(root, 'manifest.yaml'), 'w').write('\n'.join(lines) + '\n')
print('ok')
