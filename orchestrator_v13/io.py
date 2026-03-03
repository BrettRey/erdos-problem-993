from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .types import Config, WitnessRaw


def config_from_dict(d: dict[str, Any]) -> Config:
    return Config(
        N_target=int(d["N_target"]),
        N_authoritative=int(d["N_authoritative"]),
        epsilon_buffer=float(d["epsilon_buffer"]),
        eps=float(d["eps"]),
        alpha_19=float(d["alpha_19"]),
        alpha_star=float(d["alpha_star"]),
        Delta_alpha=float(d["Delta_alpha"]),
        rho_max=float(d["rho_max"]),
        lambda_max=float(d["lambda_max"]),
        K_G=int(d["K_G"]),
        K_A=int(d["K_A"]),
        K_L=int(d["K_L"]),
        kappa_P=float(d["kappa_P"]),
        kappa_A=float(d["kappa_A"]),
        kappa_B=float(d["kappa_B"]),
        kappa_C=float(d["kappa_C"]),
        kappa_M=float(d["kappa_M"]),
        K0=float(d.get("K0", 1000.0)),
    )


def witness_raw_from_dict(d: dict[str, Any]) -> WitnessRaw:
    return WitnessRaw(
        witness_id=str(d["witness_id"]),
        n=int(d["n"]),
        a=int(d["a"]),
        b=int(d["b"]),
        k=int(d.get("k", 0)),
        m=int(d.get("m", 1)),
        D=float(d["D"]),
        sum_err=float(d["sum_err"]),
        C10=float(d.get("C10", 1.0)),
        C01=float(d.get("C01", 0.0)),
        C11=float(d.get("C11", 0.0)),
        sum_all=float(d.get("sum_all", 0.0)),
        Lambda=float(d["Lambda"]),
        Gamma_L=float(d.get("Gamma_L", 0.0)),
    )


def load_fixture(path: Path) -> tuple[Config, list[WitnessRaw], dict[str, Any]]:
    obj = json.loads(path.read_text())
    config = config_from_dict(obj["config"])
    raw_log = [witness_raw_from_dict(x) for x in obj["raw_log"]]
    expected = dict(obj.get("expected", {}))
    return config, raw_log, expected
