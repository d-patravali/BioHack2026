"""Check ranking semantics on tiny examples with exactly known answers."""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from explore_inference import expected_metrics


class RankingTests(unittest.TestCase):
    def setUp(self):
        self.harmonic = np.r_[0., np.cumsum(1 / np.arange(1, 5))]

    def test_unique_top(self):
        self.assertEqual(expected_metrics(np.array([3., 2., 1.]), 0, set(),
                                          self.harmonic), [1., 1., 1., 1.])

    def test_uniform_ties(self):
        values = expected_metrics(np.zeros(3), 0, set(), self.harmonic)
        np.testing.assert_allclose(values, [(1 + 1/2 + 1/3)/3, 1/3, 1., 1.])

    def test_filter_other_positives_but_keep_target(self):
        values = expected_metrics(np.array([5., 2., 1.]), 1, {0, 1}, self.harmonic)
        self.assertEqual(values, [1., 1., 1., 1.])

    def test_tie_below_higher_score(self):
        values = expected_metrics(np.array([5., 2., 2.]), 1, set(), self.harmonic)
        np.testing.assert_allclose(values, [(1/2 + 1/3)/2, 0., 1., 1.])


if __name__ == "__main__":
    unittest.main()
