import tempfile
import unittest
from pathlib import Path

from orchestrator_v13 import emit_artifacts, run_orchestrator_v13, verify_replay
from orchestrator_v13.types import ClosureStatus, Config, PartitionId, WitnessRaw


def wr(
    witness_id,
    n,
    a,
    b,
    D,
    sum_err,
    Lambda,
    Gamma_L=0.0,
    k=0,
    m=1,
    C10=1.0,
    C01=0.0,
    C11=0.0,
    sum_all=0.0,
):
    return WitnessRaw(
        witness_id=witness_id,
        n=n,
        a=a,
        b=b,
        k=k,
        m=m,
        D=D,
        sum_err=sum_err,
        C10=C10,
        C01=C01,
        C11=C11,
        sum_all=sum_all,
        Lambda=Lambda,
        Gamma_L=Gamma_L,
    )


def default_config(N_target=24, N_authoritative=24):
    return Config(
        N_target=N_target,
        N_authoritative=N_authoritative,
        epsilon_buffer=0.0,
        eps=1e-12,
        alpha_19=0.2437206585182262,
        alpha_star=0.21034113597068071,
        Delta_alpha=0.033379522547545476,
        rho_max=1.0,
        lambda_max=0.05,
        K_G=1,
        K_A=1,
        K_L=1,
        kappa_P=1.0,
        kappa_A=1.0,
        kappa_B=1.0,
        kappa_C=1.0,
        kappa_M=1.0,
        K0=1000.0,
    )


class TestOrchestratorSmoke(unittest.TestCase):
    def test_t1_no_break_pi0(self):
        config = default_config()
        log = [wr("w1", 10, 3, 13, D=1.0, sum_err=1.1, Lambda=1.4)]
        out = run_orchestrator_v13(log, config)
        self.assertEqual(out.partition_plan.selected_partition_id, PartitionId.PI0)
        self.assertEqual(out.v13_obligations.global_closure.status, ClosureStatus.CLOSED)

        with tempfile.TemporaryDirectory() as td:
            out_dir = Path(td)
            emit_artifacts(out, out_dir)
            errors = verify_replay(log, config, out_dir)
            self.assertEqual(errors, [])

    def test_t6_nonauthoritative_unknown(self):
        config = default_config(N_target=25, N_authoritative=24)
        log = [
            wr("w1", 10, 3, 13, D=1.0, sum_err=1.1, Lambda=1.4),
            # This one is filtered out from authoritative selection.
            wr("w2", 25, 4, 18, D=1.0, sum_err=1.4, Lambda=1.1),
        ]
        out = run_orchestrator_v13(log, config)
        self.assertEqual(out.v13_obligations.global_closure.status, ClosureStatus.UNKNOWN)


if __name__ == "__main__":
    unittest.main()
