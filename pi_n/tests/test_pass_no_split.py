import unittest

from pi_n.algorithm import compute_pi
from pi_n.types import ABClass, InstanceRecord
from pi_n.verify import verify_pass_certificate


def mk(uid, a, b, Lambda, D, R_shift, sum_all):
    return InstanceRecord(
        uid=uid,
        cls=ABClass(a=a, b=b),
        Lambda=Lambda,
        D=D,
        R_shift=R_shift,
        sum_all=sum_all,
    )


class TestPassNoSplit(unittest.TestCase):
    def test_pass_without_split(self):
        log = [
            mk("u1", 2, 18, 1, 0, 1, 0),
            mk("u2", 4, 18, 2, 0, 1, 0),
        ]
        res = compute_pi(log)
        self.assertEqual(res.status, "PASS")
        self.assertEqual(res.transcript.splits, ())
        verify_pass_certificate(log, res)


if __name__ == "__main__":
    unittest.main()
