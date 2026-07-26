import unittest

import numpy as np

from hydrogenic import dislocation


class TestDislocation(unittest.TestCase):

    def test_orientation(self):
        for disloc in ["edge", "screw"]:
            for gp in ["x", "y"]:
                orientation = dislocation.get_orientation(disloc, gp)
                self.assertAlmostEqual(np.linalg.det(orientation), np.prod(np.linalg.norm(orientation, axis=-1)))

if __name__ == "__main__":
    unittest.main()
