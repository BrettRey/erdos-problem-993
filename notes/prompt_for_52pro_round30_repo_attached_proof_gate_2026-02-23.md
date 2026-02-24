# Round 30: Repo-Attached Canonical Proof Gate

You must run in the local repository for `Erdos_Problem_993`.

## Step 0 (mandatory context proof)
Run these commands first and include output exactly:

```bash
pwd
ls -l attack4_common.py scripts/verify_route1_moment_class_scan.py
python3 - <<'PY'
import attack4_common
import scripts.verify_route1_moment_class_scan as v
print('attack4_common:', attack4_common.__file__)
print('moment_scan:', v.__file__)
PY
```

If any command fails, stop and output exactly:
`INVALID_CONTEXT: repo not attached`

## Definitions (must match repository code)
Use only definitions from:
- `attack4_common.py`
- `scripts/verify_route1_moment_class_scan.py`

In particular:
- canonical validity gates:
  1) tree (`|E|=n-1`),
  2) `is_dleaf_le_1(...)`,
  3) `bridge_decomposition(..., require_dleaf=True) is not None`.
- `Threshold`, `B_max`, LP-dual step bounds, and class key `K4=(deg(P),m,lambda,mu1..mu4)` exactly as implemented.

## Required task (binary)
Return exactly one:

A) `SUCCESS`
- One concrete canonical symbolic lemma that implies `B_max >= Threshold` (or directly `R1_tail2`) for all canonical instances.
- Include full chain to `E_route1`.

B) `BLOCKED`
- Provide an explicit canonical pair `(T1,T2)` with:
  - both pass all three validity gates,
  - same `K4` (and any extra invariants you claim are required by your class),
  - different target side (`Threshold` differs OR sign of `B_max-Threshold` differs).
- Must include both g6 strings and exact numeric table.

If neither A nor B is achieved, output `INSUFFICIENT`.

## Mandatory verification payload
Include a reproducible Python snippet and output that prints:
- validity-gate booleans,
- `(m, lambda, Threshold, B_max, B_max-Threshold)`,
- for pair claims: both keys and equality check.

No heuristic analogies. No unrestricted perturbation arguments. No invented identities.

## Output format (strict)
1) `Canonical class inputs used`
2) `Lemma or explicit canonical pair`
3) `Derivation / verification` (with Step 0 output)
4) `Binary verdict: SUCCESS or BLOCKED or INSUFFICIENT`
5) `If BLOCKED/INSUFFICIENT: minimal next canonical invariant or class`
