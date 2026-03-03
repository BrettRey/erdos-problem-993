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


class TestFirstFailingBucketOrder(unittest.TestCase):
    def test_s_bucket_split_happens_before_l_bucket_fail(self):
        log = [
            # S bucket: fails but not singleton-class, so should split first.
            mk("a", 2, 18, 0, 0, 1, 10),
            mk("b", 2, 19, 0, 100, 1, 0),
            # L bucket: failing singleton class if considered immediately.
            mk("c", 4, 18, 0, 100, 1, 0),
            mk("d", 4, 18, 0, 0, 1, 100),
        ]
        res = compute_pi(log)
        self.assertEqual(res.status, "FAIL_CLASS")
        self.assertGreaterEqual(len(res.transcript.splits), 1)
        self.assertEqual(res.transcript.splits[0].bucket_index, 0)
        self.assertEqual(res.transcript.splits[0].cls, ABClass(2, 18))


if __name__ == "__main__":
    unittest.main()
