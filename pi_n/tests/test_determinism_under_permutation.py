import itertools
import unittest

from pi_n.algorithm import compute_pi
from pi_n.serialize import dumps_certificate
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


class TestDeterminismUnderPermutation(unittest.TestCase):
    def test_output_is_identical_for_permuted_input(self):
        base = [
            mk("a", 2, 18, 0, 0, 1, 10),
            mk("b", 2, 19, 0, 100, 1, 0),
            mk("c", 4, 18, 1, 0, 1, 0),
        ]
        dumps = {
            dumps_certificate(compute_pi(list(perm)))
            for perm in itertools.permutations(base)
        }
        self.assertEqual(len(dumps), 1)


if __name__ == "__main__":
    unittest.main()
