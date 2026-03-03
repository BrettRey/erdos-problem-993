from __future__ import annotations

import math
from pathlib import Path

from .frontier import (
    build_aggregates,
    build_frontier_state,
    derive_scores,
    filter_prefix,
    resolve_authority,
)
from .obligations import build_v13_obligations
from .partition import build_partition_plan, build_partitions, eval_partition, select_partition
from .queue import build_repair_queue
from .types import (
    AuthorityMode,
    BucketStats,
    ClassKey,
    ClosureStatus,
    Config,
    Diagnostics,
    FrontierState,
    GlobalAgg,
    PartitionId,
    PartitionPlan,
    RegimeCore,
    RunMeta,
    RunOutputs,
    SelectionKey,
    V13Obligations,
    WitnessRaw,
)
from .util import (
    argmax_by,
    argmin_by,
    canonical_json_dumps,
    hash_objects,
    require_not_nan,
    stable_sorted,
)


def _validate_config(config: Config) -> None:
    numeric = {
        "epsilon_buffer": config.epsilon_buffer,
        "eps": config.eps,
        "alpha_19": config.alpha_19,
        "alpha_star": config.alpha_star,
        "Delta_alpha": config.Delta_alpha,
        "rho_max": config.rho_max,
        "lambda_max": config.lambda_max,
        "kappa_P": config.kappa_P,
        "kappa_A": config.kappa_A,
        "kappa_B": config.kappa_B,
        "kappa_C": config.kappa_C,
        "kappa_M": config.kappa_M,
        "K0": config.K0,
    }
    for k, v in numeric.items():
        require_not_nan(v, k)

    if config.N_authoritative < 1:
        raise ValueError("N_authoritative must be >=1")
    if config.eps < 0:
        raise ValueError("eps must be >=0")
    if config.epsilon_buffer < 0:
        raise ValueError("epsilon_buffer must be >=0")
    if config.Delta_alpha <= 0:
        raise ValueError("Delta_alpha must be >0")
    if config.rho_max < 0 or config.lambda_max < 0:
        raise ValueError("rho_max and lambda_max must be >=0")
    if min(config.kappa_P, config.kappa_A, config.kappa_B, config.kappa_C, config.kappa_M) <= 0:
        raise ValueError("all kappa_* must be >0")


def _build_run_meta(
    raw_sorted: list[WitnessRaw],
    config: Config,
    n_auth_used: int,
    authority_mode: AuthorityMode,
    selected_partition: PartitionId | None,
) -> RunMeta:
    seed = {
        "version": "v13",
        "raw_log": raw_sorted,
        "config": config,
        "n_auth_used": n_auth_used,
        "authority_mode": authority_mode.value,
        "selected_partition": selected_partition.value if selected_partition is not None else None,
    }
    run_id = hash_objects(seed)[:32]
    return RunMeta(
        version="v13",
        run_id=run_id,
        N_target=config.N_target,
        N_authoritative=config.N_authoritative,
        N_auth_used=n_auth_used,
        authority_mode=authority_mode,
    )


def _make_no_data_outputs(
    run_meta: RunMeta,
    config: Config,
) -> RunOutputs:
    dummy_frontier = GlobalAgg(
        alpha_front=math.inf,
        lambda_front=0.0,
        local_front=0.0,
        gap=math.inf,
        Delta_global=0.0,
        w_alpha_id="",
        w_lambda_id="",
        w_local_id="",
        c_alpha=ClassKey(0, 0),
        c_lambda=ClassKey(0, 0),
        c_local=ClassKey(0, 0),
    )

    frontier_state = FrontierState(
        run_meta=run_meta,
        global_frontier=dummy_frontier,
        line_summary=tuple(),
        diagnostics=Diagnostics(
            regime_core=RegimeCore.NO_BREAK,
            a_bad=None,
            Delta_line_max=0.0,
            Delta_withinclass_max=0.0,
        ),
        provisional=None,
    )

    # rebuild as concrete dataclasses to avoid string enums above
    from .types import BucketId, ModuleTag

    stats_tuple = (
        BucketStats(BucketId.S1, ModuleTag.ENV, 0, math.inf, 0.0, 0.0, 0.0, None, None, None),
        BucketStats(BucketId.S2, ModuleTag.ENV, 0, math.inf, 0.0, 0.0, 0.0, None, None, None),
        BucketStats(BucketId.S3, ModuleTag.ENV, 0, math.inf, 0.0, 0.0, 0.0, None, None, None),
    )

    partition_plan = PartitionPlan(
        run_meta=run_meta,
        selected_partition_id=PartitionId.PI0,
        selection_key=SelectionKey(Phi=0.0, Complexity=0, PartitionIdOrder=0),
        buckets=stats_tuple,
        a2_calibration=None,
        hardstep_summary=None,
    )

    from .types import BucketObligation, GlobalClosure, ModuleTag

    obligations = V13Obligations(
        run_meta=run_meta,
        partition_id=PartitionId.PI0,
        buckets=(
            BucketObligation(BucketId.S1, ModuleTag.ENV, ClosureStatus.CLOSED, 0.0, None, None, None, None, None, None),
            BucketObligation(BucketId.S2, ModuleTag.ENV, ClosureStatus.CLOSED, 0.0, None, None, None, None, None, None),
            BucketObligation(BucketId.S3, ModuleTag.ENV, ClosureStatus.CLOSED, 0.0, None, None, None, None, None, None),
        ),
        global_closure=GlobalClosure(
            status=ClosureStatus.NO_DATA,
            Phi=0.0,
            epsilon_buffer=config.epsilon_buffer,
            reasons=("NO_DATA",),
        ),
    )

    return RunOutputs(
        frontier_state=frontier_state,
        partition_plan=partition_plan,
        repair_queue=tuple(),
        v13_obligations=obligations,
    )


def run_orchestrator_v13(raw_log: list[WitnessRaw], config: Config) -> RunOutputs:
    _validate_config(config)

    raw_sorted = stable_sorted(raw_log, key=lambda r: r.witness_id)
    n_auth, mode = resolve_authority(config)
    auth_records = filter_prefix(raw_sorted, n_auth)

    if not auth_records:
        run_meta = _build_run_meta(raw_sorted, config, n_auth, mode, None)
        return _make_no_data_outputs(run_meta, config)

    ws = derive_scores(auth_records, config)
    agg = build_aggregates(ws, config)

    parts = build_partitions(
        ws,
        c_alpha_a=agg.global_.c_alpha.a,
        c_lambda_a=agg.global_.c_lambda.a,
        c_lambda=(agg.global_.c_lambda.a, agg.global_.c_lambda.b),
    )
    evals = [eval_partition(p, ws, config) for p in parts]
    selected = select_partition(evals, config)

    run_meta = _build_run_meta(raw_sorted, config, n_auth, mode, selected.partition.id)

    provisional = None
    if mode == AuthorityMode.NON_AUTHORITATIVE:
        prov_records = filter_prefix(raw_sorted, config.N_target)
        if prov_records:
            ws_prov = derive_scores(prov_records, config)
            alpha_prov = argmin_by(
                ws_prov,
                lambda w: w.Alpha_req,
                lambda w: w.raw.witness_id,
                config.eps,
            ).Alpha_req
            lambda_prov = argmax_by(
                ws_prov,
                lambda w: w.Gamma_G,
                lambda w: w.raw.witness_id,
                config.eps,
            ).Gamma_G
            provisional = {
                "alpha_front_prov": alpha_prov,
                "lambda_front_prov": lambda_prov,
                "gap_prov": alpha_prov - lambda_prov,
                "notes": "PROVISIONAL (non-authoritative): computed using n<=N_target; not used for selection",
            }

    frontier_state = build_frontier_state(run_meta, agg, provisional)
    partition_plan = build_partition_plan(run_meta, selected, config)
    repair_queue = build_repair_queue(selected, evals, agg, config)
    obligations = build_v13_obligations(run_meta, selected, config, mode)

    # hard postconditions for authority gating
    if mode == AuthorityMode.NON_AUTHORITATIVE and obligations.global_closure.status == ClosureStatus.CLOSED:
        raise AssertionError("non-authoritative runs must never report CLOSED")

    if mode == AuthorityMode.AUTHORITATIVE and obligations.global_closure.status == ClosureStatus.CLOSED:
        if selected.Phi > config.epsilon_buffer + config.eps:
            raise AssertionError("authoritative CLOSED requires Phi <= epsilon_buffer + eps")

    return RunOutputs(
        frontier_state=frontier_state,
        partition_plan=partition_plan,
        repair_queue=repair_queue,
        v13_obligations=obligations,
    )


def emit_artifacts(out: RunOutputs, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {
        "frontier_state.json": out.frontier_state,
        "partition_plan.json": out.partition_plan,
        "repair_queue.json": out.repair_queue,
        "v13_obligations.json": out.v13_obligations,
    }
    for name, obj in files.items():
        (out_dir / name).write_bytes(canonical_json_dumps(obj))
