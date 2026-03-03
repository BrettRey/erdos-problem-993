from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any, Union

from .errors import CertificateError
from .ratio import canon_ratio
from .transcript import BASE_SPLIT
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


def _int_as_str(v: int) -> str:
    return str(v)


def _ratio_to_obj(r: Ratio) -> dict[str, str]:
    rc = canon_ratio(r)
    return {"num": _int_as_str(rc.num), "den": _int_as_str(rc.den)}


def _class_to_obj(cls: ABClass) -> list[str]:
    return [_int_as_str(cls.a), _int_as_str(cls.b)]


def _bucket_cert_to_obj(cert: BucketCert) -> dict[str, Any]:
    return {
        "alpha": _ratio_to_obj(cert.alpha),
        "lambda": _ratio_to_obj(cert.lambda_),
        "witness_lower": cert.witness_lower,
        "witness_upper": cert.witness_upper,
    }


def _transcript_to_obj(t: Transcript) -> dict[str, Any]:
    return {
        "base": t.base,
        "splits": [
            {"i": step.bucket_index, "class": _class_to_obj(step.cls)}
            for step in t.splits
        ],
    }


def to_certificate_obj(result: PiResult) -> dict[str, Any]:
    if isinstance(result, PiPass):
        out: dict[str, Any] = {
            "status": result.status,
            "partition": [list(bucket) for bucket in result.partition],
            "bucket_certs": [_bucket_cert_to_obj(c) for c in result.bucket_certs],
            "transcript": _transcript_to_obj(result.transcript),
        }
        if result.input_digest is not None:
            out["input_digest"] = result.input_digest
        return out

    if isinstance(result, PiFailClass):
        out = {
            "status": result.status,
            "obstruction_class": _class_to_obj(result.obstruction_class),
            "obstruction_bucket_index": result.obstruction_bucket_index,
            "bucket_cert": _bucket_cert_to_obj(result.bucket_cert),
            "transcript": _transcript_to_obj(result.transcript),
        }
        if result.input_digest is not None:
            out["input_digest"] = result.input_digest
        return out

    if isinstance(result, PiFailDegenerate):
        out = {
            "status": result.status,
            "failures": [
                {
                    "uid": f.uid,
                    "class": _class_to_obj(f.cls),
                    "Lambda": _int_as_str(f.Lambda),
                    "D": _int_as_str(f.D),
                    "sum_all": _int_as_str(f.sum_all),
                    "violations": list(f.violations),
                }
                for f in result.failures
            ],
        }
        if result.input_digest is not None:
            out["input_digest"] = result.input_digest
        return out

    raise CertificateError(f"unsupported result type: {type(result)}")


def dumps_certificate(result: PiResult) -> str:
    obj = to_certificate_obj(result)
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _parse_int(v: Any) -> int:
    if isinstance(v, int):
        return v
    if isinstance(v, str):
        try:
            return int(v)
        except ValueError as exc:
            raise CertificateError(f"invalid integer string: {v!r}") from exc
    raise CertificateError(f"expected integer or integer-string, got {type(v)}")


def _parse_class(v: Any) -> ABClass:
    if not isinstance(v, Sequence) or len(v) != 2:
        raise CertificateError("class must be length-2 sequence")
    return ABClass(a=_parse_int(v[0]), b=_parse_int(v[1]))


def _parse_ratio(v: Any) -> Ratio:
    if not isinstance(v, Mapping):
        raise CertificateError("ratio must be an object")
    r = Ratio(num=_parse_int(v.get("num")), den=_parse_int(v.get("den")))
    return canon_ratio(r)


def _parse_bucket_cert(v: Any) -> BucketCert:
    if not isinstance(v, Mapping):
        raise CertificateError("bucket cert must be an object")
    return BucketCert(
        alpha=_parse_ratio(v.get("alpha")),
        lambda_=_parse_ratio(v.get("lambda")),
        witness_lower=str(v.get("witness_lower")),
        witness_upper=str(v.get("witness_upper")),
    )


def _parse_transcript(v: Any) -> Transcript:
    if not isinstance(v, Mapping):
        raise CertificateError("transcript must be an object")
    base = v.get("base")
    if base != BASE_SPLIT:
        raise CertificateError(f"unsupported transcript base: {base!r}")
    raw_splits = v.get("splits", [])
    if not isinstance(raw_splits, Sequence):
        raise CertificateError("transcript.splits must be an array")

    splits: list[SplitStep] = []
    for raw in raw_splits:
        if not isinstance(raw, Mapping):
            raise CertificateError("split step must be object")
        i = _parse_int(raw.get("i"))
        cls = _parse_class(raw.get("class"))
        splits.append(SplitStep(bucket_index=i, cls=cls))

    return Transcript(base=BASE_SPLIT, splits=tuple(splits))


def _parse_failures(v: Any) -> tuple[DegenerateFailure, ...]:
    if not isinstance(v, Sequence):
        raise CertificateError("failures must be an array")
    out: list[DegenerateFailure] = []
    for row in v:
        if not isinstance(row, Mapping):
            raise CertificateError("failure entry must be object")
        raw_violations = row.get("violations", [])
        if not isinstance(raw_violations, Sequence):
            raise CertificateError("violations must be array")
        violations = tuple(str(x) for x in raw_violations)
        out.append(
            DegenerateFailure(
                uid=str(row.get("uid")),
                cls=_parse_class(row.get("class")),
                Lambda=_parse_int(row.get("Lambda")),
                D=_parse_int(row.get("D")),
                sum_all=_parse_int(row.get("sum_all")),
                violations=violations,
            )
        )
    return tuple(out)


def _parse_pass_cert(obj: Mapping[str, Any]) -> PiPass:
    raw_partition = obj.get("partition")
    if not isinstance(raw_partition, Sequence):
        raise CertificateError("partition must be array")
    partition = tuple(tuple(str(uid) for uid in bucket) for bucket in raw_partition)

    raw_bucket_certs = obj.get("bucket_certs")
    if not isinstance(raw_bucket_certs, Sequence):
        raise CertificateError("bucket_certs must be array")
    bucket_certs = tuple(_parse_bucket_cert(c) for c in raw_bucket_certs)

    transcript = _parse_transcript(obj.get("transcript"))
    digest = obj.get("input_digest")
    if digest is not None:
        digest = str(digest)

    return PiPass(
        status="PASS",
        partition=partition,
        bucket_certs=bucket_certs,
        transcript=transcript,
        input_digest=digest,
    )


def _parse_fail_class_cert(obj: Mapping[str, Any]) -> PiFailClass:
    digest = obj.get("input_digest")
    if digest is not None:
        digest = str(digest)
    return PiFailClass(
        status="FAIL_CLASS",
        obstruction_class=_parse_class(obj.get("obstruction_class")),
        obstruction_bucket_index=_parse_int(obj.get("obstruction_bucket_index")),
        bucket_cert=_parse_bucket_cert(obj.get("bucket_cert")),
        transcript=_parse_transcript(obj.get("transcript")),
        input_digest=digest,
    )


def _parse_fail_degenerate_cert(obj: Mapping[str, Any]) -> PiFailDegenerate:
    digest = obj.get("input_digest")
    if digest is not None:
        digest = str(digest)
    return PiFailDegenerate(
        status="FAIL_DEGENERATE",
        failures=_parse_failures(obj.get("failures")),
        input_digest=digest,
    )


def _parse_obj(obj: Mapping[str, Any]) -> PiResult:
    status = obj.get("status")
    if status == "PASS":
        return _parse_pass_cert(obj)
    if status == "FAIL_CLASS":
        return _parse_fail_class_cert(obj)
    if status == "FAIL_DEGENERATE":
        return _parse_fail_degenerate_cert(obj)
    raise CertificateError(f"unsupported certificate status: {status!r}")


def loads_certificate(src: Union[str, bytes, Mapping[str, Any], PiResult]) -> PiResult:
    if isinstance(src, (PiPass, PiFailClass, PiFailDegenerate)):
        return src
    if isinstance(src, bytes):
        src = src.decode("utf-8")
    if isinstance(src, str):
        obj = json.loads(src)
        if not isinstance(obj, Mapping):
            raise CertificateError("certificate JSON must decode to object")
        return _parse_obj(obj)
    if isinstance(src, Mapping):
        return _parse_obj(src)
    raise CertificateError(f"unsupported certificate source type: {type(src)}")


def canonical_log_bytes(log: Sequence[InstanceRecord]) -> bytes:
    rows: list[str] = []
    for x in sorted(log, key=lambda r: r.uid):
        row = {
            "uid": x.uid,
            "class": _class_to_obj(x.cls),
            "Lambda": _int_as_str(x.Lambda),
            "D": _int_as_str(x.D),
            "R_shift": _int_as_str(x.R_shift),
            "sum_all": _int_as_str(x.sum_all),
        }
        rows.append(json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False))
    return ("\n".join(rows) + "\n").encode("utf-8")


def compute_input_digest(log: Sequence[InstanceRecord]) -> str:
    h = hashlib.sha256()
    h.update(canonical_log_bytes(log))
    return h.hexdigest()
