from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

from .errors import InputLogError
from .ordering import argmax_uid, argmin_uid
from .ratio import canon_ratio, cmp_ratio
from .transcript import BASE_SPLIT, apply_split, build_initial_partition
from .types import (
    ABClass,
    BucketCert,
    DegenerateFailure,
    InstanceRecord,
    PiFailClass,
    PiFailDegenerate,
    PiPass,
    PiResult,
    Ratio,
    SplitStep,
    Transcript,
)


@dataclass(frozen=True)
class Prepared:
    log_sorted: tuple[InstanceRecord, ...]
    r0: tuple[InstanceRecord, ...]
    rp: tuple[InstanceRecord, ...]
    clsmap: dict[str, ABClass]
    rlower: dict[str, Ratio]
    rupper: dict[str, Ratio]


def prepare_log(log: Sequence[InstanceRecord]) -> Prepared:
    log_sorted = tuple(sorted(log, key=lambda x: x.uid))
    if len({x.uid for x in log_sorted}) != len(log_sorted):
        raise InputLogError("uids must be unique")

    r0_list: list[InstanceRecord] = []
    rp_list: list[InstanceRecord] = []
    clsmap: dict[str, ABClass] = {}
    rlower: dict[str, Ratio] = {}
    rupper: dict[str, Ratio] = {}

    for x in log_sorted:
        clsmap[x.uid] = x.cls
        if x.R_shift == 0:
            r0_list.append(x)
            continue
        rp_list.append(x)
        rlower[x.uid] = Ratio(num=x.Lambda + x.sum_all - x.D, den=x.R_shift)
        rupper[x.uid] = Ratio(num=x.sum_all - x.D, den=x.R_shift)

    return Prepared(
        log_sorted=log_sorted,
        r0=tuple(r0_list),
        rp=tuple(rp_list),
        clsmap=clsmap,
        rlower=rlower,
        rupper=rupper,
    )


def compute_degenerate_failures(r0: Sequence[InstanceRecord]) -> tuple[DegenerateFailure, ...]:
    failures: list[DegenerateFailure] = []
    for x in sorted(r0, key=lambda y: y.uid):
        violations: list[str] = []
        if x.sum_all > x.D:
            violations.append("sum_all > D")
        if x.Lambda < x.D - x.sum_all:
            violations.append("Lambda < D - sum_all")
        if violations:
            failures.append(
                DegenerateFailure(
                    uid=x.uid,
                    cls=x.cls,
                    Lambda=x.Lambda,
                    D=x.D,
                    sum_all=x.sum_all,
                    violations=tuple(violations),
                )
            )
    return tuple(failures)


def minimal_degenerate_failures(
    failures: tuple[DegenerateFailure, ...],
) -> tuple[DegenerateFailure, ...]:
    if not failures:
        return failures
    smallest = min(failures, key=lambda f: f.uid)
    return (smallest,)


def compute_bucket_cert(bucket: tuple[str, ...], prepared: Prepared) -> BucketCert:
    y = argmin_uid(bucket, prepared.rlower, prepared.clsmap)
    x = argmax_uid(bucket, prepared.rupper, prepared.clsmap)
    alpha = canon_ratio(prepared.rlower[y])
    lambda_ = canon_ratio(prepared.rupper[x])
    return BucketCert(alpha=alpha, lambda_=lambda_, witness_lower=y, witness_upper=x)


def evaluate_partition(
    partition: tuple[tuple[str, ...], ...],
    prepared: Prepared,
) -> tuple[tuple[BucketCert, ...], Optional[int]]:
    certs: list[BucketCert] = []
    first_bad: Optional[int] = None
    for i, bucket in enumerate(partition):
        cert = compute_bucket_cert(bucket, prepared)
        certs.append(cert)
        if first_bad is None and cmp_ratio(cert.alpha, cert.lambda_) < 0:
            first_bad = i
    return tuple(certs), first_bad


def compute_pi(log: Sequence[InstanceRecord], input_digest: Optional[str] = None) -> PiResult:
    prepared = prepare_log(log)
    failures = compute_degenerate_failures(prepared.r0)
    if failures:
        return PiFailDegenerate(
            status="FAIL_DEGENERATE",
            failures=minimal_degenerate_failures(failures),
            input_digest=input_digest,
        )

    sorted_positive_uids = [x.uid for x in prepared.rp]
    partition = build_initial_partition(sorted_positive_uids, prepared.clsmap)

    splits: list[SplitStep] = []

    while True:
        bucket_certs, first_bad = evaluate_partition(partition, prepared)

        if first_bad is None:
            return PiPass(
                status="PASS",
                partition=partition,
                bucket_certs=bucket_certs,
                transcript=Transcript(base=BASE_SPLIT, splits=tuple(splits)),
                input_digest=input_digest,
            )

        bad_cert = bucket_certs[first_bad]
        cls_star = prepared.clsmap[bad_cert.witness_upper]
        bad_bucket = partition[first_bad]

        if all(prepared.clsmap[uid] == cls_star for uid in bad_bucket):
            return PiFailClass(
                status="FAIL_CLASS",
                obstruction_class=cls_star,
                obstruction_bucket_index=first_bad,
                bucket_cert=bad_cert,
                transcript=Transcript(base=BASE_SPLIT, splits=tuple(splits)),
                input_digest=input_digest,
            )

        partition = apply_split(partition, first_bad, cls_star, prepared.clsmap)
        splits.append(SplitStep(bucket_index=first_bad, cls=cls_star))
