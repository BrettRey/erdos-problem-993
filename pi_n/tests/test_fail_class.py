import unittest

from pi_n.algorithm import compute_pi
from pi_n.types import ABClass, InstanceRecord
from pi_n.verify import verify_fail_class


def mk(uid, a, b, Lambda, D, R_shift, sum_all):
    return InstanceRecord(
        uid=uid,
        cls=ABClass(a=a, b=b),
        Lambda=Lambda,
        D=D,
        R_shift=R_shift,
        sum_all=sum_all,
    )


class TestFailClass(unittest.TestCase):
    def test_fail_class(self):
        log = [
            mk("a", 2, 18, 0, 100, 1, 0),
            mk("b", 2, 18, 0, 0, 1, 100),
        ]
        res = compute_pi(log)
        self.assertEqual(res.status, "FAIL_CLASS")
        self.assertEqual(res.obstruction_class, ABClass(2, 18))
        self.assertEqual(res.obstruction_bucket_index, 0)
        verify_fail_class(log, res)


if __name__ == "__main__":
    unittest.main()
