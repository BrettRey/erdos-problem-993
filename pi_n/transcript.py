from __future__ import annotations

from .errors import CertificateError
from .types import ABClass, Bucket, Partition

BASE_SPLIT = "LINE_SPLIT_A_LE_3"


def apply_split(
    partition: Partition,
    bucket_index: int,
    cls: ABClass,
    clsmap: dict[str, ABClass],
) -> Partition:
    if bucket_index < 0 or bucket_index >= len(partition):
        raise CertificateError(f"split bucket index out of range: {bucket_index}")
    target = partition[bucket_index]
    left = tuple(uid for uid in target if clsmap[uid] == cls)
    right = tuple(uid for uid in target if clsmap[uid] != cls)
    if not left or not right:
        raise CertificateError("invalid split: class side or residual side is empty")
    return partition[:bucket_index] + (left, right) + partition[bucket_index + 1 :]


def build_initial_partition(sorted_positive_uids: list[str], clsmap: dict[str, ABClass]) -> Partition:
    s_members = tuple(uid for uid in sorted_positive_uids if clsmap[uid].a <= 3)
    l_members = tuple(uid for uid in sorted_positive_uids if clsmap[uid].a >= 4)
    buckets: list[Bucket] = []
    if s_members:
        buckets.append(s_members)
    if l_members:
        buckets.append(l_members)
    return tuple(buckets)
