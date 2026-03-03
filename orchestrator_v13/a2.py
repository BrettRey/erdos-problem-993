from __future__ import annotations

import math
from collections.abc import Sequence

from .types import A2Calibration, Bucket, BucketId, Config, Witness
from .util import feq, pos_part, stable_sorted


def calibrate_a2_bucket(bucket: Bucket, config: Config) -> A2Calibration:
    xs = list(bucket.records)

    rho_min = 0.0
    for w in xs:
        if w.raw.a != 2:
            rho_min = math.inf

    if rho_min != math.inf:
        for w in xs:
            if w.R_shift == 0:
                if w.raw.D == 0 and w.raw.sum_err > 0:
                    rho_min = math.inf
                elif w.raw.D > 0 and w.raw.sum_err > w.raw.D:
                    req = (w.raw.sum_err - w.raw.D) / w.raw.D
                    rho_min = max(rho_min, req)

    rho_max = config.rho_max
    lambda_max = config.lambda_max

    feasible_domain = math.isfinite(rho_min) and (rho_min <= rho_max + config.eps)
    if not feasible_domain:
        return A2Calibration(
            bucket_id=bucket.id,
            rho_min=rho_min,
            rho_max=rho_max,
            lambda_max=lambda_max,
            rho_hat=rho_min,
            lambda_hat=math.inf,
            active_witnesses=tuple(),
            feasible_domain=False,
        )

    pairs: list[tuple[Witness, float, float]] = []
    for w in xs:
        if w.R_shift > 0:
            ai = (w.raw.sum_err - w.raw.D) / w.R_shift
            bi = w.raw.D / w.R_shift
            pairs.append((w, ai, bi))

    candidates: set[float] = {rho_min, rho_max}

    for _, ai, bi in pairs:
        if bi > 0 and ai > 0:
            rho0 = ai / bi
            rho0 = min(max(rho0, rho_min), rho_max)
            candidates.add(rho0)

    for i in range(len(pairs)):
        _, ai, bi = pairs[i]
        for j in range(i + 1, len(pairs)):
            _, aj, bj = pairs[j]
            if not feq(bi, bj, config.eps):
                rho_ij = (ai - aj) / (bi - bj)
                if rho_min - config.eps <= rho_ij <= rho_max + config.eps:
                    rho_ij = min(max(rho_ij, rho_min), rho_max)
                    candidates.add(rho_ij)

    c_list = stable_sorted(candidates, key=lambda x: x)

    def lambda_req_for(w: Witness, rho: float) -> float:
        if w.R_shift > 0:
            return pos_part(w.raw.sum_err - (1.0 + rho) * w.raw.D) / w.R_shift
        if w.raw.sum_err <= (1.0 + rho) * w.raw.D:
            return 0.0
        return math.inf

    best_rho = c_list[0]
    best_val = math.inf

    for rho in c_list:
        if rho < rho_min - config.eps or rho > rho_max + config.eps:
            continue
        cur = 0.0
        for w in xs:
            cur = max(cur, lambda_req_for(w, rho))
        if cur < best_val and not feq(cur, best_val, config.eps):
            best_val = cur
            best_rho = rho
        elif feq(cur, best_val, config.eps) and rho < best_rho:
            best_val = cur
            best_rho = rho

    active = sorted(
        [w.raw.witness_id for w in xs if feq(lambda_req_for(w, best_rho), best_val, config.eps)]
    )

    return A2Calibration(
        bucket_id=bucket.id,
        rho_min=rho_min,
        rho_max=rho_max,
        lambda_max=lambda_max,
        rho_hat=best_rho,
        lambda_hat=best_val,
        active_witnesses=tuple(active),
        feasible_domain=True,
    )
