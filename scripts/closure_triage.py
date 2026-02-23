#!/usr/bin/env python3
"""Closure triage wrapper for Route 1 / STRONG C2 lemma candidates.

Workflow:
1) Tight-witness check from existing artifacts (fast).
2) Full verifier run for the relevant lemma class.
3) Append a triage entry to notes/proof_closure_checklist_2026-02-21.md.

The script intentionally reuses existing verifiers/results and does not add
new mathematical checking logic.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LOG = REPO_ROOT / "notes" / "proof_closure_checklist_2026-02-21.md"
DEFAULT_TRIAGE_DIR = REPO_ROOT / "results" / "triage"


@dataclass
class CheckResult:
    ok: bool
    summary: str
    details: dict[str, Any]


@dataclass
class FullResult:
    ok: bool
    summary: str
    details: dict[str, Any]
    out_json: Path
    cmd: list[str]
    return_code: int


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def ensure_artifact(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing artifact: {path}")


def check_route1_tight(repo: Path) -> CheckResult:
    p = repo / "results" / "path1_direct_chain_from_existing_results_2026_02_19.json"
    ensure_artifact(p)
    data = load_json(p)
    derived = data.get("derived", {})
    lb = float(derived.get("lower_bound_from_chain", -1.0))
    implies = bool(derived.get("implies_muP_ge_m_minus_2", False))
    ok = implies and lb > 0.0
    return CheckResult(
        ok=ok,
        summary=f"lower_bound_from_chain={lb:.12f}, implies={implies}",
        details={
            "artifact": str(p.relative_to(repo)),
            "lower_bound_from_chain": lb,
            "implies_muP_ge_m_minus_2": implies,
        },
    )


def check_strong_c2_hard_ratio_tight(repo: Path) -> CheckResult:
    p = repo / "results" / "verify_strong_c2_rise_identity_2026_02_19.json"
    ensure_artifact(p)
    data = load_json(p)
    s = data.get("summary", {})
    hard_fail = int(s.get("hard_ratio_fail", 1))
    rise_fail = int(s.get("rise_fail", 1))
    comb_neg = int(s.get("combined_neg", 1))
    min_gap = float(s.get("min_hard_ratio_gap", -1.0))
    ok = hard_fail == 0 and rise_fail == 0 and comb_neg == 0 and min_gap > 0.0
    return CheckResult(
        ok=ok,
        summary=(
            f"hard_ratio_fail={hard_fail}, rise_fail={rise_fail}, "
            f"combined_neg={comb_neg}, min_hard_ratio_gap={min_gap:.12f}"
        ),
        details={
            "artifact": str(p.relative_to(repo)),
            "hard_ratio_fail": hard_fail,
            "rise_fail": rise_fail,
            "combined_neg": comb_neg,
            "min_hard_ratio_gap": min_gap,
        },
    )


def check_routeb_r_nonneg_tight(repo: Path) -> CheckResult:
    p = repo / "results" / "verify_strong_c2_route_b_pair_bounds_2026_02_19_n23.json"
    ensure_artifact(p)
    data = load_json(p)
    r_neg = int(data.get("R_neg", 1))
    r_min = int(data.get("R_min", -1))
    cabs_neg = int(data.get("cross_minus_abs_mismatch_neg", 1))
    cabs_min = int(data.get("cross_minus_abs_mismatch_min", -1))
    ok = r_neg == 0 and r_min >= 0 and cabs_neg == 0 and cabs_min >= 0
    return CheckResult(
        ok=ok,
        summary=(
            f"R_neg={r_neg}, R_min={r_min}, "
            f"cross_minus_abs_mismatch_neg={cabs_neg}, cross_minus_abs_mismatch_min={cabs_min}"
        ),
        details={
            "artifact": str(p.relative_to(repo)),
            "R_neg": r_neg,
            "R_min": r_min,
            "cross_minus_abs_mismatch_neg": cabs_neg,
            "cross_minus_abs_mismatch_min": cabs_min,
        },
    )


def check_routec_modep_tight(repo: Path) -> CheckResult:
    p_dom = repo / "results" / "prove_strong_c2_p_dominance_n20.json"
    qdrop = repo / "results" / "analyze_strong_c2_qdrop_witnesses_2026_02_19.json"
    ensure_artifact(p_dom)
    ensure_artifact(qdrop)
    d1 = load_json(p_dom)
    d2 = load_json(qdrop)

    total = int(d1.get("total", 0))
    p_ge_q = int(d1.get("p_ge_q", -1))
    min_gap = int(d1.get("min_gap", -1))
    witnesses = d2.get("witnesses", [])

    qdrop_ok = True
    worst_ratio = None
    for w in witnesses:
        ratio = float(w.get("p1_over_b1", 0.0))
        if worst_ratio is None or ratio < worst_ratio:
            worst_ratio = ratio
        if ratio <= 0.5:
            qdrop_ok = False

    ok = total > 0 and p_ge_q == total and min_gap >= 1 and qdrop_ok
    return CheckResult(
        ok=ok,
        summary=(
            f"p_ge_q={p_ge_q}/{total}, min_gap={min_gap}, "
            f"qdrop_worst_p1_over_b1={worst_ratio if worst_ratio is not None else 'NA'}"
        ),
        details={
            "artifact_p_dominance": str(p_dom.relative_to(repo)),
            "artifact_qdrop": str(qdrop.relative_to(repo)),
            "p_ge_q": p_ge_q,
            "total": total,
            "min_gap": min_gap,
            "qdrop_worst_p1_over_b1": worst_ratio,
            "qdrop_count": len(witnesses),
        },
    )


def parse_route1_full(out_json: Path) -> FullResult:
    data = load_json(out_json)
    s = data.get("summary", {})
    fail_keys = ["muP_fail", "sum_identity_fail", "chain_identity_fail", "transfer_cap_violation"]
    fails = {k: int(s.get(k, 0)) for k in fail_keys}
    ok = all(v == 0 for v in fails.values())
    summary = ", ".join(f"{k}={v}" for k, v in fails.items())
    return FullResult(
        ok=ok,
        summary=summary,
        details={"summary": fails, "checked": int(s.get("checked", 0))},
        out_json=out_json,
        cmd=[],
        return_code=0,
    )


def parse_strong_c2_hard_ratio_full(out_json: Path) -> FullResult:
    data = load_json(out_json)
    s = data.get("summary", {})
    fields = ["combined_neg", "rise_fail", "hard_ratio_fail", "identity_fail"]
    vals = {k: int(s.get(k, 0)) for k in fields}
    ok = all(vals[k] == 0 for k in ["combined_neg", "rise_fail", "hard_ratio_fail"])
    summary = ", ".join(f"{k}={v}" for k, v in vals.items())
    return FullResult(
        ok=ok,
        summary=summary,
        details={"summary": vals, "checked": int(s.get("checked", 0))},
        out_json=out_json,
        cmd=[],
        return_code=0,
    )


def parse_routeb_r_nonneg_full(out_json: Path) -> FullResult:
    data = load_json(out_json)
    vals = {
        "R_neg": int(data.get("R_neg", 0)),
        "R_min": int(data.get("R_min", -1)),
        "cross_minus_abs_mismatch_neg": int(data.get("cross_minus_abs_mismatch_neg", 0)),
    }
    ok = vals["R_neg"] == 0 and vals["R_min"] >= 0 and vals["cross_minus_abs_mismatch_neg"] == 0
    summary = ", ".join(f"{k}={v}" for k, v in vals.items())
    return FullResult(
        ok=ok,
        summary=summary,
        details={"summary": vals, "checked": int(data.get("checked", 0))},
        out_json=out_json,
        cmd=[],
        return_code=0,
    )


def parse_routec_modep_full(out_json: Path) -> FullResult:
    data = load_json(out_json)
    total = int(data.get("total", 0))
    mode_p = int(data.get("mode_P_ge_m_minus_1", -1))
    p_ge_q = int(data.get("p_ge_q", -1))
    ok = total > 0 and mode_p == total and p_ge_q == total
    summary = f"mode_P_ge_m_minus_1={mode_p}/{total}, p_ge_q={p_ge_q}/{total}"
    return FullResult(
        ok=ok,
        summary=summary,
        details={
            "summary": {
                "mode_P_ge_m_minus_1": mode_p,
                "p_ge_q": p_ge_q,
                "total": total,
            },
            "checked": total,
        },
        out_json=out_json,
        cmd=[],
        return_code=0,
    )


def ensure_triage_heading(log_path: Path) -> None:
    if not log_path.exists():
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("# Proof Closure Checklist\n\n", encoding="utf-8")
    txt = log_path.read_text(encoding="utf-8")
    if "## Triage Log" not in txt:
        if not txt.endswith("\n"):
            txt += "\n"
        txt += "\n## Triage Log\n"
        log_path.write_text(txt, encoding="utf-8")


def append_triage_log(
    log_path: Path,
    lemma_class: str,
    lemma_name: str,
    note: str | None,
    tight: CheckResult,
    full: FullResult | None,
) -> None:
    ensure_triage_heading(log_path)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")

    lines = [
        f"\n### {stamp} - {lemma_class} - {lemma_name}",
        f"- Tight check: {'PASS' if tight.ok else 'FAIL'} ({tight.summary})",
    ]
    if full is None:
        lines.append("- Full verifier: SKIPPED")
    else:
        try:
            out_ref = str(full.out_json.relative_to(REPO_ROOT))
        except ValueError:
            out_ref = str(full.out_json)
        lines.append(f"- Full verifier: {'PASS' if full.ok else 'FAIL'} ({full.summary})")
        lines.append(f"- Full artifact: `{out_ref}`")
    if note:
        lines.append(f"- Note: {note}")

    with log_path.open("a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def run_full(
    cmd: list[str],
    out_json: Path,
    parser: Callable[[Path], FullResult],
) -> FullResult:
    out_json.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(cmd, cwd=REPO_ROOT)
    if proc.returncode != 0:
        return FullResult(
            ok=False,
            summary=f"verifier exited with code {proc.returncode}",
            details={},
            out_json=out_json,
            cmd=cmd,
            return_code=proc.returncode,
        )
    if not out_json.exists():
        return FullResult(
            ok=False,
            summary=f"verifier completed but output missing: {out_json}",
            details={},
            out_json=out_json,
            cmd=cmd,
            return_code=proc.returncode,
        )
    parsed = parser(out_json)
    parsed.cmd = cmd
    parsed.return_code = proc.returncode
    return parsed


def main() -> int:
    ap = argparse.ArgumentParser(description="Run closure triage for a lemma candidate.")
    ap.add_argument(
        "--lemma-class",
        required=True,
        choices=["route1_exact_slack", "strong_c2_hard_ratio", "routeB_R_nonneg", "routeC_modeP"],
    )
    ap.add_argument("--lemma-name", required=True, help="Short identifier for the candidate lemma.")
    ap.add_argument("--note", default=None, help="Optional short note for triage log.")
    ap.add_argument("--python", default=sys.executable or "python3")
    ap.add_argument("--geng", default="/opt/homebrew/bin/geng")
    ap.add_argument("--skip-full", action="store_true", help="Only run tight-witness check.")
    ap.add_argument("--max-n", type=int, default=None, help="Override class default max n for full verifier.")
    ap.add_argument("--log-file", default=str(DEFAULT_LOG))
    ap.add_argument("--out-dir", default=str(DEFAULT_TRIAGE_DIR))
    args = ap.parse_args()

    log_path = Path(args.log_file)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Class-specific wiring.
    if args.lemma_class == "route1_exact_slack":
        tight = check_route1_tight(REPO_ROOT)
        max_n = args.max_n if args.max_n is not None else 23
        out_json = out_dir / f"triage_route1_exact_slack_{args.lemma_name}.json"
        cmd = [
            args.python,
            str(REPO_ROOT / "verify_muP_sum_of_means_2026_02_19.py"),
            "--min-n",
            "4",
            "--max-n",
            str(max_n),
            "--geng",
            args.geng,
            "--out",
            str(out_json),
            "--no-resume",
        ]
        parser = parse_route1_full
    elif args.lemma_class == "strong_c2_hard_ratio":
        tight = check_strong_c2_hard_ratio_tight(REPO_ROOT)
        max_n = args.max_n if args.max_n is not None else 23
        out_json = out_dir / f"triage_strong_c2_hard_ratio_{args.lemma_name}.json"
        cmd = [
            args.python,
            str(REPO_ROOT / "verify_strong_c2_rise_identity_2026_02_19.py"),
            "--min-n",
            "4",
            "--max-n",
            str(max_n),
            "--geng",
            args.geng,
            "--out",
            str(out_json),
        ]
        parser = parse_strong_c2_hard_ratio_full
    elif args.lemma_class == "routeB_R_nonneg":
        tight = check_routeb_r_nonneg_tight(REPO_ROOT)
        max_n = args.max_n if args.max_n is not None else 23
        out_json = out_dir / f"triage_routeB_R_nonneg_{args.lemma_name}.json"
        cmd = [
            args.python,
            str(REPO_ROOT / "verify_strong_c2_route_b_pair_bounds_2026_02_19.py"),
            "--min-n",
            "4",
            "--max-n",
            str(max_n),
            "--geng",
            args.geng,
            "--out",
            str(out_json),
        ]
        parser = parse_routeb_r_nonneg_full
    else:
        tight = check_routec_modep_tight(REPO_ROOT)
        max_n = args.max_n if args.max_n is not None else 20
        out_json = out_dir / f"triage_routeC_modeP_{args.lemma_name}.json"
        cmd = [
            args.python,
            str(REPO_ROOT / "verify_route_c_p_dominance_2026_02_19.py"),
            "--min-n",
            "4",
            "--max-n",
            str(max_n),
            "--geng",
            args.geng,
            "--out",
            str(out_json),
        ]
        parser = parse_routec_modep_full

    print(f"[triage] class={args.lemma_class} lemma={args.lemma_name}")
    print(f"[triage] tight check: {'PASS' if tight.ok else 'FAIL'} - {tight.summary}")

    full_res: FullResult | None = None
    if not args.skip_full:
        print(f"[triage] running full verifier: {' '.join(cmd)}")
        full_res = run_full(cmd, out_json, parser)
        print(f"[triage] full check: {'PASS' if full_res.ok else 'FAIL'} - {full_res.summary}")
        print(f"[triage] full artifact: {full_res.out_json}")
    else:
        print("[triage] full verifier skipped by --skip-full")

    append_triage_log(
        log_path=log_path,
        lemma_class=args.lemma_class,
        lemma_name=args.lemma_name,
        note=args.note,
        tight=tight,
        full=full_res,
    )
    print(f"[triage] log updated: {log_path}")

    if not tight.ok:
        return 2
    if full_res is not None and not full_res.ok:
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
