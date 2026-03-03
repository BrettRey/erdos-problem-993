import unittest

from pi_n.algorithm import compute_pi
from pi_n.types import ABClass, InstanceRecord


def mk(uid, a, b, Lambda, D, R_shift, sum_all):
    return InstanceRecord(
        uid=uid,
        cls=ABClass(a=a, b=b),
        Lambda=Lambda,
        D=D,
        R_shift=R_shift,
        sum_all=sum_all,
    )


class TestTieBreaks(unittest.TestCase):
    def test_tie_break_prefers_class_before_uid(self):
        log = [
            mk("a", 2, 19, 0, 0, 1, 0),
            mk("b", 2, 18, 0, 0, 1, 0),
        ]
        res = compute_pi(log)
        self.assertEqual(res.status, "PASS")
        cert = res.bucket_certs[0]
        self.assertEqual(cert.witness_lower, "b")
        self.assertEqual(cert.witness_upper, "b")

    def test_tie_break_same_class_uses_uid(self):
        log = [
            mk("z", 2, 18, 0, 0, 1, 0),
            mk("a", 2, 18, 0, 0, 1, 0),
        ]
        res = compute_pi(log)
        self.assertEqual(res.status, "PASS")
        cert = res.bucket_certs[0]
        self.assertEqual(cert.witness_lower, "a")
        self.assertEqual(cert.witness_upper, "a")


if __name__ == "__main__":
    unittest.main()
