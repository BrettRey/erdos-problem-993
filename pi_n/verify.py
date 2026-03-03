from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from .algorithm import (
    compute_degenerate_failures,
    evaluate_partition,
    minimal_degenerate_failures,
    prepare_log,
)
from .errors import CertificateError
from .ratio import cmp_ratio
from .serialize import compute_input_digest, loads_certificate
from .transcript import apply_split
from .types import PiFailClass, PiFailDegenerate, PiPass, PiResult


def _verify_digest_if_present(log: Sequence[Any], input_digest: str | None) -> None:
    if input_digest is None:
        return
    actual = compute_input_digest(log)
    if actual != input_digest:
        raise CertificateError("input digest mismatch")


def _cert_as_type(cert: Mapping[str, Any] | PiResult, expected_type: type[PiResult]) -> PiResult:
    parsed = loads_certificate(cert)
    if not isinstance(parsed, expected_type):
        raise CertificateError(
            f"certificate has wrong status/type: expected {expected_type.__name__}, got {type(parsed).__name__}"
        )
    return parsed


def verify_pass_certificate(log: Sequence[Any], cert: Mapping[str, Any] | PiPass) -> None:
    parsed = _cert_as_type(cert, PiPass)
    _verify_digest_if_present(log, parsed.input_digest)

    prepared = prepare_log(log)
    failures = compute_degenerate_failures(prepared.r0)
    if failures:
        raise CertificateError("degenerate failures exist; PASS certificate invalid")

    partition = tuple(parsed.partition)

    # Reconstruct from canonical base rather than trusting cert.partition.
    partition = tuple()
    sorted_positive_uids = [x.uid for x in prepared.rp]
    from .transcript import build_initial_partition

    partition = build_initial_partition(sorted_positive_uids, prepared.clsmap)

    for idx, step in enumerate(parsed.transcript.splits):
        bucket_certs, first_bad = evaluate_partition(partition, prepared)
        if first_bad is None:
            raise CertificateError(f"transcript has extra split at step {idx}")
        if step.bucket_index != first_bad:
            raise CertificateError(
                f"split {idx} has non-canonical bucket index: got {step.bucket_index}, expected {first_bad}"
            )

        cls_star = prepared.clsmap[bucket_certs[first_bad].witness_upper]
        if step.cls != cls_star:
            raise CertificateError(
                f"split {idx} has non-canonical class: got {step.cls}, expected {cls_star}"
            )

        bad_bucket = partition[first_bad]
        if all(prepared.clsmap[uid] == cls_star for uid in bad_bucket):
            raise CertificateError(
                f"split {idx} attempts split on singleton-class failing bucket; expected FAIL_CLASS"
            )

        partition = apply_split(partition, step.bucket_index, step.cls, prepared.clsmap)

    final_certs, first_bad = evaluate_partition(partition, prepared)
    if first_bad is not None:
        raise CertificateError("final partition still has negative-gap bucket")

    if parsed.partition != partition:
        raise CertificateError("partition mismatch against replayed transcript")

    if parsed.bucket_certs != final_certs:
        raise CertificateError("bucket certificates mismatch against replayed partition")


def verify_fail_class(log: Sequence[Any], cert: Mapping[str, Any] | PiFailClass) -> None:
    parsed = _cert_as_type(cert, PiFailClass)
    _verify_digest_if_present(log, parsed.input_digest)

    prepared = prepare_log(log)
    failures = compute_degenerate_failures(prepared.r0)
    if failures:
        raise CertificateError("degenerate failures exist; FAIL_CLASS certificate invalid")

    from .transcript import build_initial_partition

    partition = build_initial_partition([x.uid for x in prepared.rp], prepared.clsmap)

    for idx, step in enumerate(parsed.transcript.splits):
        bucket_certs, first_bad = evaluate_partition(partition, prepared)
        if first_bad is None:
            raise CertificateError(f"transcript has extra split at step {idx}")
        if step.bucket_index != first_bad:
            raise CertificateError(
                f"split {idx} has non-canonical bucket index: got {step.bucket_index}, expected {first_bad}"
            )
        cls_star = prepared.clsmap[bucket_certs[first_bad].witness_upper]
        if step.cls != cls_star:
            raise CertificateError(
                f"split {idx} has non-canonical class: got {step.cls}, expected {cls_star}"
            )

        bad_bucket = partition[first_bad]
        if all(prepared.clsmap[uid] == cls_star for uid in bad_bucket):
            raise CertificateError(
                f"split {idx} attempts split on singleton-class failing bucket; expected FAIL_CLASS"
            )

        partition = apply_split(partition, step.bucket_index, step.cls, prepared.clsmap)

    final_certs, first_bad = evaluate_partition(partition, prepared)
    if first_bad is None:
        raise CertificateError("replayed partition does not fail; FAIL_CLASS certificate invalid")

    if parsed.obstruction_bucket_index != first_bad:
        raise CertificateError(
            f"obstruction bucket index mismatch: got {parsed.obstruction_bucket_index}, expected {first_bad}"
        )

    failing_cert = final_certs[first_bad]
    cls_star = prepared.clsmap[failing_cert.witness_upper]
    if cls_star != parsed.obstruction_class:
        raise CertificateError("obstruction class mismatch")

    bad_bucket = partition[first_bad]
    if not all(prepared.clsmap[uid] == parsed.obstruction_class for uid in bad_bucket):
        raise CertificateError("failing bucket is not singleton-class at claimed obstruction")

    if cmp_ratio(failing_cert.alpha, failing_cert.lambda_) >= 0:
        raise CertificateError("claimed obstruction does not have negative gap")

    if parsed.bucket_cert != failing_cert:
        raise CertificateError("bucket_cert mismatch against replayed failing bucket")


def verify_fail_degenerate(log: Sequence[Any], cert: Mapping[str, Any] | PiFailDegenerate) -> None:
    parsed = _cert_as_type(cert, PiFailDegenerate)
    _verify_digest_if_present(log, parsed.input_digest)

    prepared = prepare_log(log)
    actual_all = compute_degenerate_failures(prepared.r0)
    if not actual_all:
        raise CertificateError("no degenerate failures exist")

    expected = minimal_degenerate_failures(actual_all)
    if parsed.failures != expected:
        raise CertificateError("degenerate failure payload is not canonical minimal witness set")

    lookup = {x.uid: x for x in prepared.r0}
    for f in parsed.failures:
        if f.uid not in lookup:
            raise CertificateError(f"degenerate witness uid not present: {f.uid}")
        rec = lookup[f.uid]
        if rec.R_shift != 0:
            raise CertificateError(f"degenerate witness has nonzero R_shift: {f.uid}")
        if rec.cls != f.cls or rec.Lambda != f.Lambda or rec.D != f.D or rec.sum_all != f.sum_all:
            raise CertificateError(f"degenerate witness payload mismatch: {f.uid}")

        expected_violations = []
        if rec.sum_all > rec.D:
            expected_violations.append("sum_all > D")
        if rec.Lambda < rec.D - rec.sum_all:
            expected_violations.append("Lambda < D - sum_all")
        if tuple(expected_violations) != tuple(f.violations):
            raise CertificateError(f"degenerate witness violations mismatch: {f.uid}")
