"""Small exact cross-checks; not a full W7 distance computation."""
import itertools
import random
import unittest
from exact_distance import solve_projection_task


def cycle_columns(n):
    return [(1 << i ^ 1 << ((i + 1) % n), int(i == 0)) for i in range(n)]


def brute(columns, observable, start, stop, weight):
    for size in range(1, min(weight, len(columns)) + 1):
        for ids in itertools.combinations(range(len(columns)), size):
            if not start <= min(ids) < stop:
                continue
            d = o = 0
            for i in ids:
                d ^= columns[i][0]
                o ^= columns[i][1]
            if d == 0 and o & 1 << observable:
                return True
    return False


class ExactSearchTests(unittest.TestCase):
    def test_distance_six_and_seven(self):
        six = solve_projection_task(cycle_columns(6), observable=0, start=0, stop=6, seconds=2)
        seven = solve_projection_task(cycle_columns(7), observable=0, start=0, stop=7, seconds=2)
        self.assertEqual(six['status'], 'SAT')
        self.assertEqual(six['weight'], 6)
        self.assertEqual(seven['status'], 'UNSAT')

    def test_zero_budget_is_unknown(self):
        out = solve_projection_task(cycle_columns(7), observable=0, start=0, stop=7, seconds=0)
        self.assertEqual(out['status'], 'UNKNOWN')

    def test_shards_and_at_most_not_exactly_six(self):
        rng = random.Random(89430)
        for case in range(25):
            n = rng.randrange(3, 11)
            columns = [(rng.randrange(32), rng.randrange(4)) for _ in range(n)]
            weight = rng.randrange(1, 7)
            for observable in (0, 1):
                for shard in range(3):
                    start, stop = n * shard // 3, n * (shard + 1) // 3
                    out = solve_projection_task(columns, observable=observable, start=start,
                                                stop=stop, max_weight=weight, seconds=2)
                    self.assertEqual(out['status'] == 'SAT',
                                     brute(columns, observable, start, stop, weight))
                    self.assertIn(out['status'], ('SAT', 'UNSAT'))


if __name__ == '__main__':
    unittest.main()
