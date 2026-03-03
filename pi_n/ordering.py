from __future__ import annotations

from .ratio import cmp_ratio
from .types import ABClass, Ratio, UID


def cls_lt(c1: ABClass, c2: ABClass) -> bool:
    return (c1.a, c1.b) < (c2.a, c2.b)


def _cls_key(cls: ABClass) -> tuple[int, int]:
    return (cls.a, cls.b)


def argmin_uid(
    members: tuple[UID, ...],
    ratio_map: dict[UID, Ratio],
    clsmap: dict[UID, ABClass],
) -> UID:
    best = members[0]
    for uid in members[1:]:
        c = cmp_ratio(ratio_map[uid], ratio_map[best])
        if c < 0:
            best = uid
            continue
        if c == 0:
            if _cls_key(clsmap[uid]) < _cls_key(clsmap[best]):
                best = uid
                continue
            if _cls_key(clsmap[uid]) == _cls_key(clsmap[best]) and uid < best:
                best = uid
    return best


def argmax_uid(
    members: tuple[UID, ...],
    ratio_map: dict[UID, Ratio],
    clsmap: dict[UID, ABClass],
) -> UID:
    best = members[0]
    for uid in members[1:]:
        c = cmp_ratio(ratio_map[uid], ratio_map[best])
        if c > 0:
            best = uid
            continue
        if c == 0:
            if _cls_key(clsmap[uid]) < _cls_key(clsmap[best]):
                best = uid
                continue
            if _cls_key(clsmap[uid]) == _cls_key(clsmap[best]) and uid < best:
                best = uid
    return best
