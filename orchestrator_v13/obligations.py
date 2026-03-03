from __future__ import annotations

from collections.abc import Sequence

from .types import (
    AuthorityMode,
    BucketObligation,
    ClosureStatus,
    Config,
    GlobalClosure,
    ModuleTag,
    PartitionEval,
    RunMeta,
    V13Obligations,
)


def _sorted_reasons(reasons: Sequence[str]) -> tuple[str, ...]:
    return tuple(sorted(reasons))


def check_closure(
    authority_mode: AuthorityMode,
    selected: PartitionEval,
    buckets_out: Sequence[BucketObligation],
    config: Config,
) -> tuple[ClosureStatus, float, tuple[str, ...]]:
    phi = selected.Phi
    reasons: list[str] = []

    env_ok = True
    a2_ok = True
    hs_ok = True

    for st in selected.stats:
        if st.module == ModuleTag.ENV and st.deficit_s > config.epsilon_buffer + config.eps:
            env_ok = False
            reasons.append(f"ENV bucket {st.bucket_id.value} deficit > buffer")

        if st.module == ModuleTag.A2_LOCAL:
            if selected.a2 is None or (not selected.a2.feasible_domain):
                a2_ok = False
                reasons.append("A2 bucket infeasible domain (rho_min>rho_max or invalid a)")
            elif selected.a2.lambda_hat > config.lambda_max + config.epsilon_buffer + config.eps:
                a2_ok = False
                reasons.append("A2 bucket lambda_hat > lambda_max + buffer")

        if st.module == ModuleTag.HARDSTEP and st.size > 0:
            hs = selected.hardstep
            if hs is None or len(hs.entries) == 0:
                hs_ok = False
                reasons.append("HARDSTEP bucket empty hardstep_list")
            else:
                ids = {e.witness_id for e in hs.entries}
                for wid in hs.must_include:
                    if wid not in ids:
                        hs_ok = False
                        reasons.append(f"HARDSTEP bucket missing must_include witness {wid}")

    phi_ok = phi <= config.epsilon_buffer + config.eps
    if not phi_ok:
        reasons.append("Phi(partition) > buffer")

    if authority_mode == AuthorityMode.NON_AUTHORITATIVE:
        if not phi_ok:
            return ClosureStatus.OPEN, phi, _sorted_reasons(reasons)
        rs = list(reasons)
        rs.append("Non-authoritative run: cannot assert CLOSED")
        return ClosureStatus.UNKNOWN, phi, _sorted_reasons(rs)

    if env_ok and a2_ok and hs_ok and phi_ok:
        return ClosureStatus.CLOSED, phi, _sorted_reasons(reasons)
    return ClosureStatus.OPEN, phi, _sorted_reasons(reasons)


def build_v13_obligations(
    run_meta: RunMeta,
    selected: PartitionEval,
    config: Config,
    authority_mode: AuthorityMode,
) -> V13Obligations:
    out: list[BucketObligation] = []

    for st in selected.stats:
        status = ClosureStatus.OPEN
        if st.size == 0 or st.deficit_s <= config.epsilon_buffer + config.eps:
            status = ClosureStatus.CLOSED

        constants = None
        frontier_witnesses = None
        a2_params = None
        active_witnesses = None
        hardstep_list = None
        must_include = None

        if st.module == ModuleTag.ENV and st.size > 0:
            constants = {
                "alpha": st.alpha_s,
                "lambda": st.lambda_s,
                "local_cap": st.local_s,
            }
            frontier_witnesses = {
                "w_alpha": st.w_alpha_s or "",
                "w_lambda": st.w_lambda_s or "",
                "w_local": st.w_local_s or "",
            }

        if st.module == ModuleTag.A2_LOCAL and selected.a2 is not None:
            a2 = selected.a2
            a2_params = {
                "rho_hat": a2.rho_hat,
                "lambda_hat": a2.lambda_hat,
                "rho_min": a2.rho_min,
                "rho_max": a2.rho_max,
                "lambda_max": a2.lambda_max,
            }
            active_witnesses = a2.active_witnesses

        if st.module == ModuleTag.HARDSTEP and selected.hardstep is not None:
            hs = selected.hardstep
            hardstep_list = hs.entries
            must_include = hs.must_include
            if st.size > 0:
                ids = {e.witness_id for e in hs.entries}
                if len(hs.entries) == 0 or any(wid not in ids for wid in hs.must_include):
                    status = ClosureStatus.OPEN

        out.append(
            BucketObligation(
                bucket_id=st.bucket_id,
                module=st.module,
                status=status,
                deficit=st.deficit_s,
                constants=constants,
                frontier_witnesses=frontier_witnesses,
                a2_params=a2_params,
                active_witnesses=active_witnesses,
                hardstep_list=hardstep_list,
                must_include=must_include,
            )
        )

    global_status, phi, reasons = check_closure(authority_mode, selected, out, config)

    return V13Obligations(
        run_meta=run_meta,
        partition_id=selected.partition.id,
        buckets=(out[0], out[1], out[2]),
        global_closure=GlobalClosure(
            status=global_status,
            Phi=phi,
            epsilon_buffer=config.epsilon_buffer,
            reasons=reasons,
        ),
    )
