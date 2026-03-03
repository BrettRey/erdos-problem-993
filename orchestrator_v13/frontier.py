from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Sequence

from .types import (
    Aggregates,
    AuthorityMode,
    ClassAgg,
    ClassKey,
    Config,
    Diagnostics,
    FrontierState,
    GlobalAgg,
    LineAgg,
    RegimeCore,
    RunMeta,
    Witness,
    WitnessRaw,
)
from .util import argmax_by, argmin_by, feq, pos_part, stable_sorted


def resolve_authority(config: Config) -> tuple[int, AuthorityMode]:
    n_auth = min(config.N_target, config.N_authoritative)
    mode = (
        AuthorityMode.AUTHORITATIVE
        if config.N_target <= config.N_authoritative
        else AuthorityMode.NON_AUTHORITATIVE
    )
    return n_auth, mode


def filter_prefix(records: Sequence[WitnessRaw], n_auth: int) -> list[WitnessRaw]:
    out: list[WitnessRaw] = []
    for r in stable_sorted(records, key=lambda x: x.witness_id):
        if r.n <= n_auth and 0 <= r.k < r.m:
            out.append(r)
    return out


def derive_scores(records: Sequence[WitnessRaw], config: Config) -> list[Witness]:
    out: list[Witness] = []
    for r in records:
        r_shift = r.C10 + r.C01 + r.C11
        defect = pos_part(r.sum_err - r.D)

        if r_shift > 0:
            gamma_g = defect / r_shift
        else:
            gamma_g = 0.0 if defect == 0 else math.inf

        t = r.Lambda + r.sum_all - r.D
        if r_shift > 0:
            alpha_req = t / r_shift
        else:
            alpha_req = math.inf if t >= 0 else -math.inf

        gamma_m_star = alpha_req - config.alpha_star

        if alpha_req == -math.inf:
            gamma_d = math.inf
        elif alpha_req == math.inf:
            gamma_d = 0.0
        else:
            gamma_d = pos_part(config.alpha_19 - alpha_req) / config.Delta_alpha

        out.append(
            Witness(
                raw=r,
                R_shift=r_shift,
                Defect=defect,
                Gamma_G=gamma_g,
                Alpha_req=alpha_req,
                Gamma_M_star=gamma_m_star,
                Gamma_D=gamma_d,
            )
        )
    return out


def _id(w: Witness) -> str:
    return w.raw.witness_id


def build_aggregates(ws: Sequence[Witness], config: Config) -> Aggregates:
    if not ws:
        raise ValueError("build_aggregates requires non-empty ws")

    ws_sorted = tuple(stable_sorted(ws, key=lambda w: w.raw.witness_id))

    w_alpha = argmin_by(ws_sorted, lambda w: w.Alpha_req, _id, config.eps)
    w_lambda = argmax_by(ws_sorted, lambda w: w.Gamma_G, _id, config.eps)
    w_local = argmax_by(ws_sorted, lambda w: w.raw.Gamma_L, _id, config.eps)

    global_ = GlobalAgg(
        alpha_front=w_alpha.Alpha_req,
        lambda_front=w_lambda.Gamma_G,
        local_front=w_local.raw.Gamma_L,
        gap=w_alpha.Alpha_req - w_lambda.Gamma_G,
        Delta_global=pos_part(w_lambda.Gamma_G - w_alpha.Alpha_req),
        w_alpha_id=w_alpha.raw.witness_id,
        w_lambda_id=w_lambda.raw.witness_id,
        w_local_id=w_local.raw.witness_id,
        c_alpha=ClassKey(w_alpha.raw.a, w_alpha.raw.b),
        c_lambda=ClassKey(w_lambda.raw.a, w_lambda.raw.b),
        c_local=ClassKey(w_local.raw.a, w_local.raw.b),
    )

    class_map: dict[ClassKey, list[Witness]] = defaultdict(list)
    for w in ws_sorted:
        class_map[ClassKey(w.raw.a, w.raw.b)].append(w)

    class_keys = stable_sorted(class_map.keys(), key=lambda ck: (ck.a, ck.b))
    classes: dict[ClassKey, ClassAgg] = {}

    for ck in class_keys:
        xs = class_map[ck]
        w_a = argmin_by(xs, lambda w: w.Alpha_req, _id, config.eps)
        w_g = argmax_by(xs, lambda w: w.Gamma_G, _id, config.eps)
        w_l = argmax_by(xs, lambda w: w.raw.Gamma_L, _id, config.eps)
        w_margin = argmin_by(xs, lambda w: w.Gamma_M_star, _id, config.eps)
        w_drift = argmax_by(xs, lambda w: w.Gamma_D, _id, config.eps)
        classes[ck] = ClassAgg(
            key=ck,
            alpha_c=w_a.Alpha_req,
            lambda_c=w_g.Gamma_G,
            local_c=w_l.raw.Gamma_L,
            margin_c_star=w_margin.Gamma_M_star,
            drift_c=w_drift.Gamma_D,
            w_alpha_id=w_a.raw.witness_id,
            w_lambda_id=w_g.raw.witness_id,
            w_local_id=w_l.raw.witness_id,
        )

    line_to_cks: dict[int, list[ClassKey]] = defaultdict(list)
    for ck in class_keys:
        line_to_cks[ck.a].append(ck)

    line_keys = stable_sorted(line_to_cks.keys(), key=lambda a: a)
    lines: dict[int, LineAgg] = {}

    for a in line_keys:
        cks = line_to_cks[a]

        def alpha_key(ck: ClassKey) -> tuple[float, int, int]:
            return (classes[ck].alpha_c, ck.a, ck.b)

        def lambda_key(ck: ClassKey) -> tuple[float, int, int]:
            return (-classes[ck].lambda_c, ck.a, ck.b)

        def local_key(ck: ClassKey) -> tuple[float, int, int]:
            return (-classes[ck].local_c, ck.a, ck.b)

        c_alpha = min(cks, key=alpha_key)
        c_lambda = min(cks, key=lambda_key)
        c_local = min(cks, key=local_key)

        alpha_line = classes[c_alpha].alpha_c
        lambda_line = classes[c_lambda].lambda_c
        local_line = classes[c_local].local_c
        delta_line = pos_part(max(lambda_line, local_line) - alpha_line)

        lines[a] = LineAgg(
            a=a,
            alpha_line=alpha_line,
            lambda_line=lambda_line,
            local_line=local_line,
            Delta_line=delta_line,
            c_alpha_line=c_alpha,
            c_lambda_line=c_lambda,
            c_local_line=c_local,
        )

    delta_withinclass_max = 0.0
    for ck in class_keys:
        class_delta = pos_part(max(classes[ck].lambda_c, classes[ck].local_c) - classes[ck].alpha_c)
        if class_delta > delta_withinclass_max and not feq(class_delta, delta_withinclass_max, config.eps):
            delta_withinclass_max = class_delta

    delta_line_max = 0.0
    a_bad = None
    for a in line_keys:
        d = lines[a].Delta_line
        if d > delta_line_max and not feq(d, delta_line_max, config.eps):
            delta_line_max = d
            a_bad = a
        elif feq(d, delta_line_max, config.eps) and d > 0 and a_bad is not None and a < a_bad:
            a_bad = a

    if global_.Delta_global <= config.eps:
        regime = RegimeCore.NO_BREAK
        a_bad = None
    else:
        if delta_line_max <= config.eps:
            regime = RegimeCore.CROSS_LINE_MISMATCH
            a_bad = None
        else:
            regime = RegimeCore.INTERNAL_LINE_DEFICIT

    diagnostics = Diagnostics(
        regime_core=regime,
        a_bad=a_bad,
        Delta_line_max=delta_line_max,
        Delta_withinclass_max=delta_withinclass_max,
    )

    return Aggregates(global_=global_, classes=classes, lines=lines, diagnostics=diagnostics)


def build_frontier_state(
    run_meta: RunMeta,
    agg: Aggregates,
    provisional: dict[str, object] | None,
) -> FrontierState:
    line_summary = tuple(agg.lines[a] for a in sorted(agg.lines.keys()))
    return FrontierState(
        run_meta=run_meta,
        global_frontier=agg.global_,
        line_summary=line_summary,
        diagnostics=agg.diagnostics,
        provisional=provisional,
    )
