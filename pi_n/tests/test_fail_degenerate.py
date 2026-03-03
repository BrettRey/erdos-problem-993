import unittest

from pi_n.algorithm import compute_pi
from pi_n.types import ABClass, InstanceRecord
from pi_n.verify import verify_fail_degenerate


def mk(uid, a, b, Lambda, D, R_shift, sum_all):
    return InstanceRecord(
        uid=uid,
        cls=ABClass(a=a, b=b),
        Lambda=Lambda,
        D=D,
        R_shift=R_shift,
        sum_all=sum_all,
    )


class TestFailDegenerate(unittest.TestCase):
    def test_fail_degenerate(self):
        log = [mk("z", 3, 14, 0, 0, 0, 1)]
        res = compute_pi(log)
        self.assertEqual(res.status, "FAIL_DEGENERATE")
        self.assertEqual(res.failures[0].uid, "z")
        self.assertIn("sum_all > D", res.failures[0].violations)
        verify_fail_degenerate(log, res)


if __name__ == "__main__":
    unittest.main()
