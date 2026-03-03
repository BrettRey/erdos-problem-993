#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


# Config aligned to orchestrator_v13 test/default settings.
DEFAULT_CONFIG = {
    "N_target": 25,
    "N_authoritative": 25,
    "epsilon_buffer": 0.0,
    "eps": 1e-12,
    "alpha_19": 0.2437206585182262,
    "alpha_star": 0.21034113597068071,
    "Delta_alpha": 0.033379522547545476,
    "rho_max": 1.0,
    "lambda_max": 0.05,
    "K_G": 16,
    "K_A": 16,
    "K_L": 16,
    "kappa_P": 1.0,
    "kappa_A": 1.0,
    "kappa_B": 1.0,
    "kappa_C": 1.0,
    "kappa_M": 1.0,
    "K0": 1000.0,
}


def _short_hash(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:12]


def _load_json(path: Path) -> dict[str, Any]:
    with path.open() as f:
        obj = json.load(f)
    if not isinstance(obj, dict):
        raise ValueError(f"expected object in {path}")
    return obj


def _iter_lambda_rows(obj: dict[str, Any], source_name: str):
    # Global lambda witness
    lw = obj.get("lambda_wit")
    if isinstance(lw, dict):
        yield source_name + ":lambda_wit", lw

    # Per-class pair witnesses
    by_pair = obj.get("by_pair")
    if isinstance(by_pair, dict):
        for pair_key in sorted(by_pair.keys()):
            info = by_pair[pair_key]
            if not isinstance(info, dict):
                continue
            wit = info.get("wit")
            if isinstance(wit, dict):
                yield source_name + f":by_pair:{pair_key}", wit


def _row_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        int(row["n"]),
        str(row["g6"]),
        int(row["root"]),
        int(row["step"]),
        int(row["k"]),
        int(row["a"]),
        int(row["b"]),
    )


def _to_witness_raw(row: dict[str, Any], source_tag: str) -> dict[str, Any]:
    # Semantics from scan_modal_lambda_frontier.py:
    #   X = Lambda - D
    #   R = R_shift
    #   sum_all = sum_s err_s
    #   need = max(0, sum_all - D)
    # For orchestrator fields we map exactly on these semantics.
    n = int(row["n"])
    a = int(row["a"])
    b = int(row["b"])
    k = int(row["k"])
    root = int(row["root"])
    g6 = str(row["g6"])
    x = float(row["X"])
    d = float(row["D"])
    r_shift = float(row["R"])
    sum_all = float(row["sum_all"])

    # Deterministic unique id across merged sources.
    wid = f"lam_n{n}_a{a}_b{b}_k{k}_r{root}_{_short_hash(g6)}"

    return {
        "witness_id": wid,
        "n": n,
        "a": a,
        "b": b,
        "k": k,
        # We only need prefix-valid rows for orchestrator filtering.
        "m": k + 1,
        "D": d,
        # In lambda frontier rows, sum_all is the computed total error mass.
        # We use it for both fields to keep Gamma_G consistent with lambda-needed geometry.
        "sum_err": sum_all,
        "C10": r_shift,
        "C01": 0.0,
        "C11": 0.0,
        "sum_all": sum_all,
        "Lambda": x + d,
        "Gamma_L": 0.0,
        "_source": source_tag,
    }


def build_dataset(lambda_files: list[Path], n25_snapshot: Path | None) -> dict[str, Any]:
    collected: list[tuple[str, dict[str, Any]]] = []

    for p in lambda_files:
        obj = _load_json(p)
        for src, row in _iter_lambda_rows(obj, p.name):
            collected.append((src, row))

    if n25_snapshot is not None:
        obj = _load_json(n25_snapshot)
        fw = obj.get("frontier", {})
        if isinstance(fw, dict):
            lw = fw.get("lambda_witness")
            if isinstance(lw, dict):
                collected.append((n25_snapshot.name + ":frontier:lambda_witness", lw))

    # Deduplicate by structural witness key; deterministic keep-first after source sort.
    collected.sort(key=lambda t: (t[0], _row_key(t[1])))
    seen: set[tuple[Any, ...]] = set()
    rows: list[dict[str, Any]] = []
    for src, row in collected:
        key = _row_key(row)
        if key in seen:
            continue
        seen.add(key)
        wr = _to_witness_raw(row, src)
        rows.append(wr)

    rows.sort(key=lambda r: r["witness_id"])

    # Drop internal source field from emitted raw_log, but preserve mapping in meta.
    source_map = {r["witness_id"]: r["_source"] for r in rows}
    raw_log = []
    for r in rows:
        out = dict(r)
        out.pop("_source", None)
        raw_log.append(out)

    n_values = sorted({int(r["n"]) for r in raw_log})
    cfg = dict(DEFAULT_CONFIG)
    if n_values:
        cfg["N_target"] = max(n_values)
        cfg["N_authoritative"] = max(n_values)

    return {
        "config": cfg,
        "raw_log": raw_log,
        "meta": {
            "builder": "build_orchestrator_input_from_lambda_frontiers.py",
            "sources": [str(p) for p in lambda_files]
            + ([str(n25_snapshot)] if n25_snapshot is not None else []),
            "witness_count": len(raw_log),
            "n_values": n_values,
            "source_map": source_map,
            "notes": [
                "Derived from lambda-frontier witness rows only.",
                "Exact script semantics used: X=Lambda-D, R=R_shift, sum_all from err-sum.",
                "C10,C01,C11 are represented as (R_shift,0,0) since orchestrator consumes only R_shift=C10+C01+C11.",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Build orchestrator_v13 input from lambda frontier witness JSONs")
    parser.add_argument(
        "--lambda-file",
        action="append",
        required=True,
        help="Path to lambda_frontier_modal_*.json (repeatable)",
    )
    parser.add_argument(
        "--n25-snapshot",
        default="results/n25_modal_frontier_authoritative_2026-03-03.json",
        help="Optional n=25 authoritative frontier snapshot with lambda_witness",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output JSON path",
    )
    args = parser.parse_args()

    lambda_files = [Path(p) for p in args.lambda_file]
    n25_snapshot = Path(args.n25_snapshot) if args.n25_snapshot else None
    if n25_snapshot is not None and not n25_snapshot.exists():
        n25_snapshot = None

    out_obj = build_dataset(lambda_files=lambda_files, n25_snapshot=n25_snapshot)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out_obj, indent=2))

    print(f"wrote {out_path}")
    print(f"witness_count={out_obj['meta']['witness_count']}")
    print(f"n_values={out_obj['meta']['n_values']}")


if __name__ == "__main__":
    main()
