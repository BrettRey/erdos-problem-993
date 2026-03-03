from __future__ import annotations

import math
from collections.abc import Sequence

from .a2 import calibrate_a2_bucket
from .hardstep import build_hardstep_list
from .types import (
    Bucket,
    BucketId,
    BucketStats,
    Config,
    ModuleTag,
    PARTITION_ID_ORDER,
    PartitionEval,
    PartitionId,
    PartitionPlan,
    PartitionSpec,
    SelectionKey,
    Witness,
)
from .util import argmax_by, argmin_by, bucket_order, cmp_phi, feq, pos_part


def build_partitions(ws: Sequence[Witness], c_alpha_a: int, c_lambda_a: int, c_lambda: tuple[int, int]) -> list[PartitionSpec]:
    has_a2 = any(w.raw.a == 2 for w in ws)

    pi0 = PartitionSpec(
        id=PartitionId.PI0,
        pred_S1=lambda _w: True,
        pred_S2=lambda _w: False,
        module_S1=ModuleTag.ENV,
        module_S2=ModuleTag.ENV,
        module_S3=ModuleTag.ENV,
        pred_S1_str="all",
        pred_S2_str="empty",
        pred_S3_str="empty",
    )

    pi1 = PartitionSpec(
        id=PartitionId.PI1,
        pred_S1=lambda w, aa=c_alpha_a: w.raw.a == aa,
        pred_S2=lambda w, la=c_lambda_a, aa=c_alpha_a: (w.raw.a == la) and (w.raw.a != aa),
        module_S1=ModuleTag.ENV,
        module_S2=ModuleTag.ENV,
        module_S3=ModuleTag.ENV,
        pred_S1_str=f"a=={c_alpha_a}",
        pred_S2_str=f"a=={c_lambda_a} and not S1",
        pred_S3_str="rest",
    )

    pi2 = PartitionSpec(
        id=PartitionId.PI2,
        pred_S1=lambda w: w.raw.a == 2,
        pred_S2=lambda w: w.raw.a != 2,
        module_S1=ModuleTag.A2_LOCAL,
        module_S2=ModuleTag.ENV,
        module_S3=ModuleTag.ENV,
        pred_S1_str="a==2",
        pred_S2_str="a!=2",
        pred_S3_str="empty",
    )

    pi3 = PartitionSpec(
        id=PartitionId.PI3,
        pred_S1=lambda w, cl=c_lambda: (w.raw.a == cl[0] and w.raw.b == cl[1]),
        pred_S2=lambda w, cl=c_lambda: not (w.raw.a == cl[0] and w.raw.b == cl[1]),
        module_S1=ModuleTag.HARDSTEP,
        module_S2=ModuleTag.ENV,
        module_S3=ModuleTag.ENV,
        pred_S1_str=f"class==({c_lambda[0]},{c_lambda[1]})",
        pred_S2_str="not S1",
        pred_S3_str="empty",
    )

    pi4 = PartitionSpec(
        id=PartitionId.PI4,
        pred_S1=lambda w: w.raw.a == 2,
        pred_S2=lambda w, cl=c_lambda: (w.raw.a == cl[0] and w.raw.b == cl[1]) and (w.raw.a != 2),
        module_S1=ModuleTag.A2_LOCAL,
        module_S2=ModuleTag.HARDSTEP,
        module_S3=ModuleTag.ENV,
        pred_S1_str="a==2",
        pred_S2_str=f"class==({c_lambda[0]},{c_lambda[1]}) and not S1",
        pred_S3_str="rest",
    )

    out = [pi0, pi1]
    if has_a2:
        out.append(pi2)
        if c_lambda[0] != 2:
            out.append(pi4)
    out.append(pi3)

    return sorted(out, key=lambda p: PARTITION_ID_ORDER[p.id])


def assign_buckets(part: PartitionSpec, ws: Sequence[Witness]) -> tuple[Bucket, Bucket, Bucket]:
    s1: list[Witness] = []
    s2: list[Witness] = []
    s3: list[Witness] = []
    for w in ws:
        if part.pred_S1(w):
            s1.append(w)
        elif part.pred_S2(w):
            s2.append(w)
        else:
            s3.append(w)

    return (
        Bucket(id=BucketId.S1, module=part.module_S1, records=tuple(s1)),
        Bucket(id=BucketId.S2, module=part.module_S2, records=tuple(s2)),
        Bucket(id=BucketId.S3, module=part.module_S3, records=tuple(s3)),
    )


def _id(w: Witness) -> str:
    return w.raw.witness_id


def _env_deficit(alpha_s: float, lambda_s: float, local_s: float) -> float:
    # Pathological infinities should remain visible to partition selection.
    if math.isnan(alpha_s) or math.isnan(lambda_s) or math.isnan(local_s):
        raise ValueError("NaN in ENV frontier statistics")
    if math.isinf(lambda_s) or math.isinf(local_s):
        return math.inf
    if alpha_s == math.inf:
        return 0.0
    if alpha_s == -math.inf:
        return math.inf
    return pos_part(max(lambda_s, local_s) - alpha_s)


def compute_bucket_stats(bucket: Bucket, config: Config) -> BucketStats:
    xs = list(bucket.records)
    if not xs:
        return BucketStats(
            bucket_id=bucket.id,
            module=bucket.module,
            size=0,
            alpha_s=math.inf,
            lambda_s=0.0,
            local_s=0.0,
            deficit_s=0.0,
            w_alpha_s=None,
            w_lambda_s=None,
            w_local_s=None,
        )

    w_a = argmin_by(xs, lambda w: w.Alpha_req, _id, config.eps)
    w_g = argmax_by(xs, lambda w: w.Gamma_G, _id, config.eps)
    w_l = argmax_by(xs, lambda w: w.raw.Gamma_L, _id, config.eps)

    alpha_s = w_a.Alpha_req
    lambda_s = w_g.Gamma_G
    local_s = w_l.raw.Gamma_L

    if bucket.module == ModuleTag.ENV:
        deficit = _env_deficit(alpha_s, lambda_s, local_s)
    else:
        deficit = 0.0

    return BucketStats(
        bucket_id=bucket.id,
        module=bucket.module,
        size=len(xs),
        alpha_s=alpha_s,
        lambda_s=lambda_s,
        local_s=local_s,
        deficit_s=deficit,
        w_alpha_s=w_a.raw.witness_id,
        w_lambda_s=w_g.raw.witness_id,
        w_local_s=w_l.raw.witness_id,
    )


def eval_partition(part: PartitionSpec, ws: Sequence[Witness], config: Config) -> PartitionEval:
    buckets = assign_buckets(part, ws)
    stats = [compute_bucket_stats(b, config) for b in buckets]

    a2 = None
    hs = None

    for i, b in enumerate(buckets):
        if b.module == ModuleTag.A2_LOCAL:
            a2 = calibrate_a2_bucket(b, config)
            deficit = pos_part(a2.lambda_hat - config.lambda_max)
            s = stats[i]
            stats[i] = BucketStats(
                bucket_id=s.bucket_id,
                module=s.module,
                size=s.size,
                alpha_s=s.alpha_s,
                lambda_s=s.lambda_s,
                local_s=s.local_s,
                deficit_s=deficit,
                w_alpha_s=s.w_alpha_s,
                w_lambda_s=s.w_lambda_s,
                w_local_s=s.w_local_s,
            )

    for b in buckets:
        if b.module == ModuleTag.HARDSTEP:
            hs = build_hardstep_list(b, config)
            break

    phi = 0.0
    for s in stats:
        phi = max(phi, s.deficit_s)

    bottleneck = tuple(
        sorted(
            [s.bucket_id for s in stats if feq(s.deficit_s, phi, config.eps) and phi > 0],
            key=bucket_order,
        )
    )

    nonempty = sum(1 for s in stats if s.size > 0)
    module_count = sum(1 for s in stats if s.size > 0 and s.module != ModuleTag.ENV)
    hs_size = len(hs.entries) if hs is not None else 0
    complexity = (nonempty - 1) + module_count + hs_size

    return PartitionEval(
        partition=part,
        buckets=(buckets[0], buckets[1], buckets[2]),
        stats=(stats[0], stats[1], stats[2]),
        a2=a2,
        hardstep=hs,
        Phi=phi,
        Complexity=complexity,
        bottleneck_bucket_ids=bottleneck,
    )


def select_partition(evals: Sequence[PartitionEval], config: Config) -> PartitionEval:
    best = evals[0]
    for e in evals[1:]:
        c = cmp_phi(e.Phi, best.Phi, config.eps)
        if c < 0:
            best = e
        elif c == 0:
            if e.Complexity < best.Complexity:
                best = e
            elif e.Complexity == best.Complexity:
                if PARTITION_ID_ORDER[e.partition.id] < PARTITION_ID_ORDER[best.partition.id]:
                    best = e
    return best


def build_partition_plan(run_meta, selected: PartitionEval, config: Config) -> PartitionPlan:
    hardstep_summary = None
    if selected.hardstep is not None:
        hs = selected.hardstep
        hardstep_summary = {
            "K_G": config.K_G,
            "K_A": config.K_A,
            "K_L": config.K_L,
            "hardstep_list_size": len(hs.entries),
            "must_include": list(hs.must_include),
        }

    key = SelectionKey(
        Phi=selected.Phi,
        Complexity=selected.Complexity,
        PartitionIdOrder=PARTITION_ID_ORDER[selected.partition.id],
    )

    return PartitionPlan(
        run_meta=run_meta,
        selected_partition_id=selected.partition.id,
        selection_key=key,
        buckets=selected.stats,
        a2_calibration=selected.a2,
        hardstep_summary=hardstep_summary,
    )
