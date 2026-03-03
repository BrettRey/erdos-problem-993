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


class TestPassOneSplit(unittest.TestCase):
    def test_pass_with_single_split(self):
        log = [
            mk("a", 2, 18, 0, 0, 1, 10),
            mk("b", 2, 19, 0, 100, 1, 0),
        ]
        res = compute_pi(log)
        self.assertEqual(res.status, "PASS")
        self.assertEqual(len(res.transcript.splits), 1)
        self.assertEqual(res.transcript.splits[0].bucket_index, 0)
        self.assertEqual(res.transcript.splits[0].cls, ABClass(2, 18))
        verify_pass_certificate(log, res)


if __name__ == "__main__":
    unittest.main()
