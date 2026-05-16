import unittest
import sys
from unittest.mock import MagicMock

# Robust mocking to handle submodules and third-party dependencies
mda = MagicMock()
mda_lib = MagicMock()
mda_lib_distances = MagicMock()
mda_lib.distances = mda_lib_distances
mda.lib = mda_lib

sys.modules['MDAnalysis'] = mda
sys.modules['MDAnalysis.lib'] = mda_lib
sys.modules['MDAnalysis.lib.distances'] = mda_lib_distances

numpy_mock = MagicMock()
# Make numpy mocks behave reasonably for compute_persistence
sys.modules['numpy'] = numpy_mock
numpy_mock.mean.side_effect = lambda x, **kwargs: float(sum(x) / len(x)) if x else 0.0
numpy_mock.max.side_effect = lambda x, **kwargs: int(max(x)) if x else 0

sys.modules['networkx'] = MagicMock()
sys.modules['scipy'] = MagicMock()
sys.modules['scipy.spatial'] = MagicMock()
sys.modules['scipy.spatial.distance'] = MagicMock()
sys.modules['scipy.cluster'] = MagicMock()
sys.modules['scipy.cluster.hierarchy'] = MagicMock()
sys.modules['similaritymeasures'] = MagicMock()

# Now import the function to test
from gephyra.analysis import compute_persistence

class TestComputePersistence(unittest.TestCase):
    def test_empty_list(self):
        mean_p, max_p = compute_persistence([], 10, 1)
        self.assertEqual(mean_p, 0.0)
        self.assertEqual(max_p, 0)

    def test_single_frame(self):
        mean_p, max_p = compute_persistence([5], 10, 1)
        self.assertEqual(mean_p, 1.0)
        self.assertEqual(max_p, 1)

    def test_continuous_sequence(self):
        # Frame indices: 1, 2, 3. Stride: 1.
        # Run: [1, 2, 3] -> length 3
        mean_p, max_p = compute_persistence([1, 2, 3], 10, 1)
        self.assertEqual(mean_p, 3.0)
        self.assertEqual(max_p, 3)

    def test_discontinuous_sequence(self):
        # Frame indices: 1, 2, 4, 5, 6. Stride: 1.
        # Runs: [1, 2] (len 2), [4, 5, 6] (len 3)
        # Mean: (2+3)/2 = 2.5
        # Max: 3
        mean_p, max_p = compute_persistence([1, 2, 4, 5, 6], 10, 1)
        self.assertEqual(mean_p, 2.5)
        self.assertEqual(max_p, 3)

    def test_stride_greater_than_one(self):
        # Frame indices: 0, 2, 4. Stride: 2.
        # Run: [0, 2, 4] -> length 3
        mean_p, max_p = compute_persistence([0, 2, 4], 10, 2)
        self.assertEqual(mean_p, 3.0)
        self.assertEqual(max_p, 3)

    def test_discontinuous_with_stride(self):
        # Frame indices: 0, 2, 6, 8, 10. Stride: 2.
        # Runs: [0, 2] (len 2), [6, 8, 10] (len 3)
        # Mean: 2.5, Max: 3
        mean_p, max_p = compute_persistence([0, 2, 6, 8, 10], 10, 2)
        self.assertEqual(mean_p, 2.5)
        self.assertEqual(max_p, 3)

    def test_all_single_frame_runs(self):
        # Frame indices: 0, 2, 4. Stride: 1.
        # Runs: [0] (len 1), [2] (len 1), [4] (len 1)
        # Mean: 1.0, Max: 1
        mean_p, max_p = compute_persistence([0, 2, 4], 10, 1)
        self.assertEqual(mean_p, 1.0)
        self.assertEqual(max_p, 1)

if __name__ == '__main__':
    unittest.main()
