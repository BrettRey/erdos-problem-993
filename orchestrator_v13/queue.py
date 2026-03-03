from __future__ import annotations

import math
from collections.abc import Sequence

from .types import (
    ACTION_TYPE_ORDER,
    ActionType,
    Aggregates,
    BucketId,
    ClassKey,
    Config,
    ModuleTag,
    PARTITION_ID_ORDER,
    PartitionEval,
    RepairQueueItem,
    RepairTarget,
)
from .util import bucket_order, pos_part


def _classes_in_bucket(eval_: PartitionEval, bucket_index: int) -> list[ClassKey]:
    b = eval_.buckets[bucket_index]
    uniq = {(w.raw.a, w.raw.b) for w in b.records}
    return [ClassKey(a, b) for a, b in sorted(uniq)]


def _sort_queue(items: list[RepairQueueItem], selected: PartitionEval, config: Config) -> list[RepairQueueItem]:
    bottleneck = set(selected.bottleneck_bucket_ids)

    def key(it: RepairQueueItem):
        p = it.priority
        if math.isnan(p):
            p = -math.inf

        positive = it.estimated_phi_reduction > config.eps
        targets_b = it.bucket_id in bottleneck if it.bucket_id is not None else False

        fr = 9
        if it.action_type in {ActionType.P, ActionType.M_A2}:
            fr = 0
        elif it.action_type == ActionType.B:
            fr = 0
        elif it.action_type == ActionType.A:
            fr = 1
        elif it.action_type == ActionType.C:
            fr = 2

        part_ord = 99
        if it.target.partition_id is not None:
            part_ord = PARTITION_ID_ORDER[it.target.partition_id]

        bucket_ord = 9
        if it.bucket_id is not None:
            bucket_ord = bucket_order(it.bucket_id)

        cls_ord = (999, 999)
        if it.target.class_key is not None:
            cls_ord = (it.target.class_key.a, it.target.class_key.b)

        witness_ord = min(it.witness_justification) if it.witness_justification else ""

        return (
            -p,
            -int(positive),
            -int(targets_b),
            ACTION_TYPE_ORDER[it.action_type],
            fr,
            it.estimated_cost,
            part_ord,
            bucket_ord,
            cls_ord,
            witness_ord,
        )

    return sorted(items, key=key)


def build_repair_queue(
    selected: PartitionEval,
    evals_all: Sequence[PartitionEval],
    agg: Aggregates,
    config: Config,
) -> tuple[RepairQueueItem, ...]:
    current_phi = selected.Phi
    bottleneck = set(selected.bottleneck_bucket_ids)

    items: list[RepairQueueItem] = []

    # P actions: exact partition switches
    for e in evals_all:
        if e.partition.id == selected.partition.id:
            continue
        red = pos_part(current_phi - e.Phi)
        cost = config.kappa_P * (1.0 + float(e.Complexity))
        pr = red / cost if cost > 0 else math.inf
        just = sorted({agg.global_.w_alpha_id, agg.global_.w_lambda_id, agg.global_.w_local_id})
        items.append(
            RepairQueueItem(
                rank=0,
                action_type=ActionType.P,
                target=RepairTarget(partition_id=e.partition.id),
                bucket_id=None,
                estimated_phi_reduction=red,
                estimated_cost=cost,
                priority=pr,
                witness_justification=tuple(just),
                notes="switch partition",
            )
        )

    # ENV local actions
    for i, st in enumerate(selected.stats):
        if st.module != ModuleTag.ENV:
            continue
        delta_bucket = st.deficit_s
        if delta_bucket <= config.epsilon_buffer + config.eps:
            continue

        cks = _classes_in_bucket(selected, i)
        if not cks:
            continue

        c_alpha = min(cks, key=lambda ck: (agg.classes[ck].alpha_c, ck.a, ck.b))
        c_lambda = min(cks, key=lambda ck: (-agg.classes[ck].lambda_c, ck.a, ck.b))

        alpha_bucket = st.alpha_s
        lambda_bucket = st.lambda_s
        local_bucket = st.local_s

        red_b = pos_part(lambda_bucket - alpha_bucket)
        cost_b = config.kappa_B * (1.0 + (lambda_bucket / max(lambda_bucket, 1e-12)))
        pr_b = red_b / cost_b if cost_b > 0 else math.inf
        items.append(
            RepairQueueItem(
                rank=0,
                action_type=ActionType.B,
                target=RepairTarget(class_key=c_lambda),
                bucket_id=st.bucket_id,
                estimated_phi_reduction=min(red_b, delta_bucket),
                estimated_cost=cost_b,
                priority=pr_b,
                witness_justification=(agg.classes[c_lambda].w_lambda_id,),
                notes="lower lambda in bucket frontier class",
            )
        )

        red_a = delta_bucket
        cost_a = config.kappa_A * (1.0 + agg.classes[c_alpha].drift_c)
        pr_a = red_a / cost_a if cost_a > 0 else math.inf
        items.append(
            RepairQueueItem(
                rank=0,
                action_type=ActionType.A,
                target=RepairTarget(class_key=c_alpha),
                bucket_id=st.bucket_id,
                estimated_phi_reduction=min(red_a, delta_bucket),
                estimated_cost=cost_a,
                priority=pr_a,
                witness_justification=(agg.classes[c_alpha].w_alpha_id,),
                notes="raise alpha in bucket frontier class",
            )
        )

        for ck in cks:
            lam = agg.classes[ck].lambda_c
            loc = agg.classes[ck].local_c
            class_def = pos_part(max(lam, loc) - alpha_bucket)
            if class_def <= config.eps:
                continue
            cost_c = config.kappa_C * (
                1.0 + max(lam, loc) / max(max(lambda_bucket, local_bucket), 1e-12)
            )
            pr_c = class_def / cost_c if cost_c > 0 else math.inf
            just = tuple(sorted({agg.classes[ck].w_lambda_id, agg.classes[ck].w_alpha_id}))
            items.append(
                RepairQueueItem(
                    rank=0,
                    action_type=ActionType.C,
                    target=RepairTarget(class_key=ck),
                    bucket_id=st.bucket_id,
                    estimated_phi_reduction=min(class_def, delta_bucket),
                    estimated_cost=cost_c,
                    priority=pr_c,
                    witness_justification=just,
                    notes="add/modify reserve channels in class",
                )
            )

    # A2 bound adjustment action
    for st in selected.stats:
        if st.module != ModuleTag.A2_LOCAL:
            continue
        if st.deficit_s <= config.epsilon_buffer + config.eps:
            continue

        targets_b = st.bucket_id in bottleneck
        red = min(st.deficit_s, current_phi) if targets_b else 0.0
        cost = config.kappa_M * (1.0 + float(st.size) / max(config.K0, 1.0))
        pr = red / cost if cost > 0 else math.inf
        active = tuple(selected.a2.active_witnesses) if selected.a2 is not None else tuple()
        items.append(
            RepairQueueItem(
                rank=0,
                action_type=ActionType.M_A2,
                target=RepairTarget(bucket_id=st.bucket_id),
                bucket_id=st.bucket_id,
                estimated_phi_reduction=red,
                estimated_cost=cost,
                priority=pr,
                witness_justification=active,
                notes="A2 bucket infeasible within bounds; adjust rho_max/lambda_max or strengthen local lemma",
            )
        )

    ordered = _sort_queue(items, selected, config)
    ranked = [
        RepairQueueItem(
            rank=i + 1,
            action_type=it.action_type,
            target=it.target,
            bucket_id=it.bucket_id,
            estimated_phi_reduction=it.estimated_phi_reduction,
            estimated_cost=it.estimated_cost,
            priority=it.priority,
            witness_justification=it.witness_justification,
            notes=it.notes,
        )
        for i, it in enumerate(ordered)
    ]
    return tuple(ranked)
