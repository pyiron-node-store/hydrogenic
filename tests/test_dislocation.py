import unittest
from unittest.mock import patch

import numpy as np
from ase import Atoms

from hydrogenic import dislocation


class TestDislocation(unittest.TestCase):

    def test_orientation(self):
        for disloc in ["edge", "screw"]:
            for gp in ["x", "y"]:
                orientation = dislocation.get_orientation(disloc, gp)
                self.assertAlmostEqual(np.linalg.det(orientation), np.prod(np.linalg.norm(orientation, axis=-1)))

    def test_burgers_vector_is_fcc(self):
        burgers_vector = dislocation.get_burgers_vector(
            lattice_parameter=4.0, dislocation_type="edge"
        )

        np.testing.assert_allclose(burgers_vector, [2.0, -2.0, 0.0])

    def test_lattice_parameter_normalizes_primitive_cell_volume(self):
        structure = Atoms("Ni")
        with (
            patch.object(dislocation, "bulk", return_value=structure) as bulk_mock,
            patch.object(dislocation, "get_tasks_for_energy_volume_curve", return_value={}),
            patch.object(dislocation, "evaluate_with_lammpslib", return_value={}),
            patch.object(
                dislocation,
                "analyse_results_for_energy_volume_curve",
                return_value={"volume_eq": 100.0},
            ),
        ):
            lattice_parameter = dislocation.get_lattice_parameter(
                element="Ni", potential_dataframe=None, cubic=False
            )

        bulk_mock.assert_called_once_with("Ni", cubic=False)
        self.assertAlmostEqual(lattice_parameter, np.cbrt(400.0))

if __name__ == "__main__":
    unittest.main()
