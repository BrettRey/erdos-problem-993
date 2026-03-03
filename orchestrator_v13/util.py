from __future__ import annotations

import dataclasses
import hashlib
import json
import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from enum import Enum
from typing import Any, Optional, TypeVar

from .types import BucketId

T = TypeVar("T")


def pos_part(x: float) -> float:
    return x if x > 0.0 else 0.0


def feq(a: float, b: float, eps: float) -> bool:
    if math.isfinite(a) and math.isfinite(b):
        return abs(a - b) <= eps
    return a == b


def cmp_phi(a: float, b: float, eps: float) -> int:
    if feq(a, b, eps):
        return 0
    return -1 if a < b else 1


def require_not_nan(x: float, label: str) -> None:
    if math.isnan(x):
        raise ValueError(f"NaN is forbidden: {label}")


def stable_sorted(items: Iterable[T], key: Callable[[T], Any]) -> list[T]:
    return sorted(items, key=key)


def argmin_by(items: Sequence[T], key_fn: Callable[[T], float], id_fn: Callable[[T], str], eps: float) -> T:
    if not items:
        raise ValueError("argmin_by requires non-empty sequence")
    best = items[0]
    best_v = key_fn(best)
    for x in items[1:]:
        v = key_fn(x)
        if (v < best_v) and (not feq(v, best_v, eps)):
            best = x
            best_v = v
        elif feq(v, best_v, eps) and id_fn(x) < id_fn(best):
            best = x
            best_v = v
    return best


def argmax_by(items: Sequence[T], key_fn: Callable[[T], float], id_fn: Callable[[T], str], eps: float) -> T:
    if not items:
        raise ValueError("argmax_by requires non-empty sequence")
    best = items[0]
    best_v = key_fn(best)
    for x in items[1:]:
        v = key_fn(x)
        if (v > best_v) and (not feq(v, best_v, eps)):
            best = x
            best_v = v
        elif feq(v, best_v, eps) and id_fn(x) < id_fn(best):
            best = x
            best_v = v
    return best


def bucket_order(bucket_id: BucketId) -> int:
    if bucket_id == BucketId.S1:
        return 0
    if bucket_id == BucketId.S2:
        return 1
    return 2


def _json_real(x: float) -> object:
    if math.isnan(x):
        raise ValueError("NaN cannot be serialized")
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x


def canonicalize_obj(obj: Any) -> Any:
    if dataclasses.is_dataclass(obj):
        obj = dataclasses.asdict(obj)
    if isinstance(obj, Enum):
        return obj.value
    if isinstance(obj, float):
        return _json_real(obj)
    if isinstance(obj, Mapping):
        out: dict[str, Any] = {}
        for k in sorted(obj.keys(), key=lambda z: str(z)):
            out[str(k)] = canonicalize_obj(obj[k])
        return out
    if isinstance(obj, (list, tuple)):
        return [canonicalize_obj(x) for x in obj]
    return obj


def canonical_json_dumps(obj: Any) -> bytes:
    canon = canonicalize_obj(obj)
    s = json.dumps(canon, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    return (s + "\n").encode("utf-8")


def canonical_json_loads(blob: bytes | str) -> Any:
    if isinstance(blob, bytes):
        blob = blob.decode("utf-8")
    return json.loads(blob)


def hash_blob(blob: bytes) -> str:
    h = hashlib.sha256()
    h.update(blob)
    return h.hexdigest()


def hash_objects(*objs: Any) -> str:
    h = hashlib.sha256()
    for obj in objs:
        h.update(canonical_json_dumps(obj))
    return h.hexdigest()
