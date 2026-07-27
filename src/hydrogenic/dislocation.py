from typing import Annotated

import flowrep as fr
import numpy as np
import pandas as pd
from ase import Atom, Atoms
from ase.build import bulk
from atomistics.calculators import (
    calc_static_with_lammpslib,
    evaluate_with_lammpslib,
    get_potential_by_name,
    optimize_positions_with_lammpslib,
)
from atomistics.workflows import (
    analyse_results_for_elastic_matrix,
    analyse_results_for_energy_volume_curve,
    get_tasks_for_elastic_matrix,
    get_tasks_for_energy_volume_curve,
)
from elaston import LinearElasticity, tools
from elaston.orientation import get_dislocation_orientation, get_shockley_partials
from pint import UnitRegistry


def get_orientation(dislocation_type: str = "screw", glide_plane: str = "y") -> list:
    """Return an FCC dislocation orientation aligned to a glide plane.

    Args:
        dislocation_type: Dislocation character, either ``"edge"`` or ``"screw"``.
        glide_plane: Box axis to which the glide-plane normal is aligned.

    Returns:
        Orientation matrix with rows corresponding to the box axes.
    """
    assert glide_plane in ["x", "y"]
    assert dislocation_type in ["edge", "screw"]
    orient = get_dislocation_orientation(
        dislocation_type=dislocation_type, crystal="fcc"
    )
    perpend = np.cross(orient["glide_plane"], orient["dislocation_line"])
    if glide_plane == "x":
        return np.array([orient["glide_plane"], -perpend, orient["dislocation_line"]])
    else:
        return np.array([perpend, orient["glide_plane"], orient["dislocation_line"]])


def get_lattice_parameter(
    element: str, potential_dataframe: pd.DataFrame, cubic=True
) -> Annotated[float, {"units": "angstrom"}]:
    """Calculate an element's equilibrium lattice parameter from an energy-volume fit.

    Args:
        element: Chemical symbol of the elemental bulk structure.
        potential_dataframe: LAMMPS potential definition used for the calculations.
        cubic: Whether to create a conventional cubic bulk cell.

    Returns:
        Equilibrium lattice parameter in angstrom.
    """
    structure = bulk(element, cubic=cubic)
    task_dict = get_tasks_for_energy_volume_curve(
        structure=structure,
        num_points=11,
        vol_range=0.05,
        axes=("x", "y", "z"),
    )
    result_dict = evaluate_with_lammpslib(
        task_dict=task_dict,
        potential_dataframe=potential_dataframe,
    )
    fit_dict = analyse_results_for_energy_volume_curve(
        output_dict=result_dict,
        task_dict=task_dict,
        fit_type="polynomial",
        fit_order=3,
    )
    lattice_parameter = (4 * fit_dict["volume_eq"] / len(structure)) ** (1 / 3)
    return lattice_parameter


def get_burgers_vector(
    lattice_parameter: Annotated[float, {"units": "angstrom"}],
    dislocation_type: str = "edge",
):
    """Return the FCC Burgers vector in crystal coordinates.

    Args:
        lattice_parameter: Equilibrium lattice parameter in angstrom.
        dislocation_type: Dislocation character, either ``"edge"`` or ``"screw"``.

    Returns:
        Three-component Burgers vector in angstrom.
    """
    direction = get_dislocation_orientation(dislocation_type, crystal="fcc")[
        "burgers_vector"
    ]
    burgers_vector = (
        lattice_parameter
        * np.asarray(direction)
        / np.linalg.norm(direction)
        / np.sqrt(2)
    )
    return burgers_vector


def get_elastic_matrix(
    fit_dict: dict,
) -> Annotated[
    np.ndarray, {"shape": (6, 6), "units": "gigapascal", "label": "elastic_matrix"}
]:
    """Extract the elastic matrix from fitted elastic-constant results.

    Args:
        fit_dict: Results dictionary returned by elastic-matrix analysis.

    Returns:
        Six-by-six elastic matrix in gigapascal.
    """
    return fit_dict["elastic_matrix"]


def evaluate_lammps_for_elastic_matrix(
    structure: Atoms,
    potential_name="1995--Angelo-J-E--Ni-Al-H--LAMMPS--ipr1",
    num_point=5,
    eps_range=0.005,
) -> Annotated[np.ndarray, {"shape": (6, 6), "units": "gigapascal"}]:
    """Calculate a structure's elastic matrix with the selected LAMMPS potential.

    Args:
        structure: Atomic structure for which to calculate elastic constants.
        potential_name: Name of the LAMMPS potential to use.
        num_point: Number of strain points for each elastic deformation.
        eps_range: Maximum strain magnitude for the elastic fit.

    Returns:
        Six-by-six elastic matrix in gigapascal.
    """
    potential_dataframe = get_potential_by_name(potential_name=potential_name)
    task_dict, sym_dict = get_tasks_for_elastic_matrix(
        structure=structure,
        num_of_point=num_point,
        eps_range=eps_range,
    )
    result_dict = evaluate_with_lammpslib(
        task_dict=task_dict,
        potential_dataframe=potential_dataframe,
    )
    fit_dict, sym_dict = analyse_results_for_elastic_matrix(
        output_dict=result_dict, sym_dict=sym_dict
    )
    elastic_matrix = get_elastic_matrix(fit_dict)
    return elastic_matrix


def rotate_elastic_tensor(elastic_matrix: np.ndarray, orientation: np.ndarray) -> np.ndarray:
    """Rotate an elastic matrix from crystal to box coordinates.

    Args:
        elastic_matrix: Six-by-six elastic matrix in gigapascal.
        orientation: Crystal-to-box orientation matrix.

    Returns:
        Six-by-six elastic matrix in box coordinates in gigapascal.
    """
    medium = LinearElasticity(C_tensor=elastic_matrix, orientation=orientation)
    elastic_matrix_rotated = medium.get_elastic_tensor(voigt=True, rotate=True)
    return elastic_matrix_rotated


@fr.workflow
def get_elastic_tensor(
    element="Ni",
    cubic=True,
    potential_name="1995--Angelo-J-E--Ni-Al-H--LAMMPS--ipr1",
    num_point=5,
    eps_range=0.005,
    orientation: list | np.ndarray | None = None,
) -> Annotated[
    np.ndarray, {"shape": (6, 6), "units": "gigapascal", "label": "elastic_matrix"}
]:
    """Calculate the elastic matrix for a bulk elemental structure.

    Args:
        element: Chemical symbol of the elemental bulk structure.
        cubic: Whether to create a conventional cubic bulk cell.
        potential_name: Name of the LAMMPS potential to use.
        num_point: Number of strain points for each elastic deformation.
        eps_range: Maximum strain magnitude for the elastic fit.

    Returns:
        Six-by-six elastic matrix in gigapascal.
    """
    structure = bulk(element, cubic=cubic)
    C = evaluate_lammps_for_elastic_matrix(
        structure=structure,
        potential_name=potential_name,
        num_point=num_point,
        eps_range=eps_range,
    )
    elastic_matrix = rotate_elastic_tensor(C, orientation)
    return elastic_matrix


def get_partial_burgers_vectors(
    burgers_vector: Annotated[np.ndarray, {"shape": (3,), "units": "angstrom"}],
    orientation: np.ndarray | list,
) -> Annotated[np.ndarray, {"shape": (2, 3), "units": "angstrom"}]:
    """Split a Burgers vector into Shockley partials in box coordinates.

    Args:
        burgers_vector: Perfect-dislocation Burgers vector in angstrom.
        orientation: Crystal-to-box orientation matrix.

    Returns:
        Two three-component partial Burgers vectors in angstrom.
    """
    shockley_partials = get_shockley_partials(burgers_vector)
    burgers_vectors = tools.crystal_to_box(shockley_partials, orientation=orientation)
    return burgers_vectors


def get_dislocation_distance(
    elastic_matrix: Annotated[
        np.ndarray, {"shape": (6, 6), "units": "gigapascal"}
        ],
    burgers_vectors: Annotated[np.ndarray, {"shape": (2, 3), "units": "angstrom"}],
    x_min: Annotated[float, {"units": "angstrom"}] = -10,
    x_max: Annotated[float, {"units": "angstrom"}] = 10,
    n_x: int = 100,
    sfe: Annotated[float, {"units": "millijoule / meter**2"}] = 90,
) -> Annotated[float, {"units": "angstrom", "label": "dislocation_distance"}]:
    """Find the partial separation where glide force equals stacking-fault energy.

    Args:
        medium: Linear-elastic medium used to calculate stress and force.
        burgers_vectors: Partial Burgers vectors in angstrom.
        x_min: Lower bound of the separation search interval in angstrom.
        x_max: Upper bound of the separation search interval in angstrom.
        n_x: Number of separation values to evaluate.
        sfe: Stacking-fault energy in millijoule per square meter.

    Returns:
        Partial-dislocation separation in angstrom.
    """
    ureg = UnitRegistry()
    medium = LinearElasticity(C_tensor=elastic_matrix * ureg.gigapascal)
    x = np.linspace(x_min, x_max, n_x)[:, None] * [1, 0] * ureg.angstrom
    stress = medium.get_dislocation_stress(x, burgers_vectors[0] * ureg.angstrom)
    f = (
        medium.get_dislocation_force(
            stress,
            glide_plane=[0, 1, 0],
            burgers_vector=burgers_vectors[1] * ureg.angstrom,
        )
        .to("millijoule/meter**2")
        .magnitude
    )
    F = (f - sfe)[:, 0]
    X = x.magnitude[:, 0]
    return X[np.abs(F).argmin()]


def get_hydrogen_structure(element="Ni", n_repeat=3):
    """Build a repeated cubic elemental cell with one interstitial hydrogen atom.

    Args:
        element: Chemical symbol of the host crystal.
        n_repeat: Number of repetitions along each cell axis.

    Returns:
        Atomic structure containing the host crystal and one hydrogen atom.
    """
    bulk_structure = bulk(element, cubic=True)
    structure = bulk_structure.repeat(n_repeat) + Atom(
        symbol="H", position=[0, 0, 0.5 * bulk_structure.cell[0, 0]]
    )
    return structure


def get_dipole_tensor(
    structure, potential_dataframe
) -> Annotated[np.ndarray, {"shape": (3, 3), "units": "eV"}]:
    """Relax a hydrogen-containing structure and calculate its elastic dipole tensor.

    Args:
        structure: Atomic structure containing the hydrogen defect.
        potential_dataframe: LAMMPS potential definition used for relaxation.

    Returns:
        Three-by-three elastic dipole tensor in electron volts.
    """
    relaxed_structure = optimize_positions_with_lammpslib(
        structure, potential_dataframe=potential_dataframe
    )
    result = calc_static_with_lammpslib(
        relaxed_structure, potential_dataframe=potential_dataframe
    )
    ureg = UnitRegistry()
    dipole_tensor = (
        (result["stress"] * ureg.bar * result["volume"] * ureg.angstrom**3)
        .to("eV")
        .magnitude
    )
    return dipole_tensor


def linspace(
    x_min: Annotated[float, {"units": "angstrom"}],
    x_max: Annotated[float, {"units": "angstrom"}],
    n_points: int,
) -> Annotated[np.ndarray, {"units": "angstrom"}]:
    """Return evenly spaced positions between two distances.

    Args:
        x_min: First position in angstrom.
        x_max: Last position in angstrom.
        n_points: Number of positions to generate.

    Returns:
        One-dimensional array of positions in angstrom.
    """
    return np.linspace(x_min, x_max, n_points)


def create_mesh(
    x: Annotated[np.ndarray, {"units": "angstrom"}],
) -> Annotated[np.ndarray, {"units": "angstrom"}]:
    """Create a flattened two-dimensional Cartesian mesh from one coordinate axis.

    Args:
        x: One-dimensional coordinate axis in angstrom.

    Returns:
        Flattened array of two-dimensional mesh coordinates in angstrom.
    """
    mesh = np.meshgrid(x, x, indexing="ij")
    mesh = np.stack(mesh, axis=-1).reshape(-1, 2)
    return mesh


def get_strain_field(
    elastic_matrix: Annotated[np.ndarray, {"shape": (6, 6), "units": "gigapascal"}],
    mesh: Annotated[np.ndarray, {"units": "angstrom"}],
    d_dislocations: Annotated[np.ndarray, {"units": "angstrom"}],
    burgers_vectors: Annotated[np.ndarray, {"shape": (2, 3), "units": "angstrom"}],
) -> np.ndarray:
    """Calculate the combined strain field of two partial dislocations.

    Args:
        medium: Linear-elastic medium used to calculate strain.
        mesh: Two-dimensional coordinates at which to evaluate strain in angstrom.
        d_dislocations: Separation vector between the partial dislocations in angstrom.
        burgers_vectors: Partial Burgers vectors in angstrom.

    Returns:
        Strain tensor evaluated at every mesh coordinate.
    """
    ureg = UnitRegistry()
    medium = LinearElasticity(C_tensor=elastic_matrix * ureg.gigapascal)
    strain = medium.get_dislocation_strain(
        (mesh - np.array([0.5, 0]) * d_dislocations) * ureg.angstrom,
        burgers_vector=burgers_vectors[0] * ureg.angstrom
    )
    strain += medium.get_dislocation_strain(
        (mesh + np.array([0.5, 0]) * d_dislocations) * ureg.angstrom,
        burgers_vector=burgers_vectors[1] * ureg.angstrom
    )
    return strain


def get_binding_energy_field(
    dipole_tensor: Annotated[np.ndarray, {"shape": (3, 3), "units": "eV"}],
    strain: np.ndarray,
):
    """Calculate and reshape the hydrogen binding-energy field from strain.

    Args:
        dipole_tensor: Three-by-three elastic dipole tensor in electron volts.
        strain: Strain tensor evaluated at each point of a square mesh.

    Returns:
        Two-dimensional hydrogen binding-energy field in electron volts.
    """
    binding_energy = -(dipole_tensor * strain).sum(axis=(-1, -2))
    n_x = int(np.sqrt(len(binding_energy)))
    binding_energy = binding_energy.reshape(n_x, n_x)
    return binding_energy


@fr.workflow
def get_hydrogen_binding(
    element: str = "Ni",
    cubic: bool = True,
    dislocation_type: str = "edge",
    potential_name="1995--Angelo-J-E--Ni-Al-H--LAMMPS--ipr1",
    n_repeat: int = 3,
    x_min: Annotated[float, {"units": "angstrom"}] = -10,
    x_max: Annotated[float, {"units": "angstrom"}] = 10,
    n_x: int = 100,
):
    """Calculate the hydrogen binding-energy field around an FCC dislocation.

    Args:
        element: Chemical symbol of the FCC host crystal.
        cubic: Whether to create a conventional cubic bulk cell.
        dislocation_type: Dislocation character, either ``"edge"`` or ``"screw"``.
        potential_name: Name of the LAMMPS potential to use.
        n_repeat: Number of repetitions along each hydrogen-structure cell axis.
        x_min: Lower bound of the field coordinate range in angstrom.
        x_max: Upper bound of the field coordinate range in angstrom.
        n_x: Number of coordinates along each field axis.

    Returns:
        Coordinate array in angstrom and a two-dimensional binding-energy field in
        electron volts.
    """
    structure = get_hydrogen_structure(element=element, n_repeat=n_repeat)
    orientation = get_orientation(dislocation_type=dislocation_type)
    potential_dataframe = get_potential_by_name(potential_name=potential_name)
    dipole_tensor = get_dipole_tensor(structure, potential_dataframe)
    elastic_matrix = get_elastic_tensor(
        element=element, cubic=cubic, potential_name=potential_name, orientation=orientation
    )
    lattice_parameter = get_lattice_parameter(element, potential_dataframe, cubic=cubic)
    burgers_vector = get_burgers_vector(lattice_parameter, dislocation_type)
    burgers_vectors = get_partial_burgers_vectors(
        burgers_vector, orientation=orientation
    )
    d_dislocations = get_dislocation_distance(
        elastic_matrix, burgers_vectors, x_min, x_max, n_x
    )
    x = linspace(x_min, x_max, n_x)
    mesh = create_mesh(x)
    strain_field = get_strain_field(elastic_matrix, mesh, d_dislocations, burgers_vectors)
    binding_energy_field = get_binding_energy_field(dipole_tensor, strain_field)
    return x, binding_energy_field
