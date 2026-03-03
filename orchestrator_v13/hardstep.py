from __future__ import annotations

import math
from collections.abc import Sequence

from .types import Bucket, Config, HardstepEntry, HardstepList, Witness
from .util import argmax_by, argmin_by


def _id(w: Witness) -> str:
    return w.raw.witness_id


def build_hardstep_list(bucket: Bucket, config: Config) -> HardstepList:
    xs = list(bucket.records)
    must: set[str] = set()

    if xs:
        w_a = argmin_by(xs, lambda w: w.Alpha_req, _id, config.eps)
        w_g = argmax_by(xs, lambda w: w.Gamma_G, _id, config.eps)
        w_l = argmax_by(xs, lambda w: w.raw.Gamma_L, _id, config.eps)
        must.update([w_a.raw.witness_id, w_g.raw.witness_id, w_l.raw.witness_id])

    for w in xs:
        if w.Gamma_G == math.inf:
            must.add(w.raw.witness_id)
        if w.raw.Gamma_L == math.inf:
            must.add(w.raw.witness_id)
        if w.Alpha_req == -math.inf:
            must.add(w.raw.witness_id)

    must_list = tuple(sorted(must))

    top_g = sorted(xs, key=lambda w: (-w.Gamma_G, w.raw.witness_id))[: config.K_G]
    top_a = sorted(xs, key=lambda w: (w.Alpha_req, w.raw.witness_id))[: config.K_A]
    top_l = sorted(xs, key=lambda w: (-w.raw.Gamma_L, w.raw.witness_id))[: config.K_L]

    hardset: set[str] = set(must)
    hardset.update(w.raw.witness_id for w in top_g)
    hardset.update(w.raw.witness_id for w in top_a)
    hardset.update(w.raw.witness_id for w in top_l)

    top_g_ids = {w.raw.witness_id for w in top_g}
    top_a_ids = {w.raw.witness_id for w in top_a}
    top_l_ids = {w.raw.witness_id for w in top_l}

    def tier(wid: str) -> int:
        if wid in must:
            return 0
        if wid in top_g_ids:
            return 1
        if wid in top_a_ids:
            return 2
        if wid in top_l_ids:
            return 3
        return 4

    ordered_ids = sorted(hardset, key=lambda wid: (tier(wid), wid))
    index = {w.raw.witness_id: w for w in xs}

    entries: list[HardstepEntry] = []
    for wid in ordered_ids:
        w = index[wid]
        flags: list[str] = []
        if w.Gamma_G == math.inf:
            flags.append("Gamma_G=+inf")
        if w.raw.Gamma_L == math.inf:
            flags.append("Gamma_L=+inf")
        if w.Alpha_req == -math.inf:
            flags.append("Alpha_req=-inf")

        entries.append(
            HardstepEntry(
                witness_id=wid,
                a=w.raw.a,
                b=w.raw.b,
                n=w.raw.n,
                k=w.raw.k,
                Gamma_G=w.Gamma_G,
                Alpha_req=w.Alpha_req,
                Gamma_L=w.raw.Gamma_L,
                flags=tuple(flags),
            )
        )

    return HardstepList(bucket_id=bucket.id, must_include=must_list, entries=tuple(entries))
