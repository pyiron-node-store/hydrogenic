from typing import Annotated
import flowrep as fr
from elaston import LinearElasticity, tools
from elaston.orientation import get_shockley_partials, get_dislocation_orientation
import numpy as np
from pint import UnitRegistry
from ase import Atom, Atoms
from atomistics.calculators import (
    optimize_positions_with_lammpslib,
    calc_static_with_lammpslib,
)
import os
import pandas as pd
from ase.build import bulk
from atomistics.calculators import evaluate_with_lammpslib, get_potential_by_name
from atomistics.workflows import (
    analyse_results_for_elastic_matrix,
    analyse_results_for_energy_volume_curve,
    get_tasks_for_elastic_matrix,
    get_tasks_for_energy_volume_curve,
)


def get_orientation(dislocation_type: str = "screw", glide_plane: str = "y") -> list:
    assert glide_plane in ["x", "y"]
    assert dislocation_type in ["edge", "screw"]
    orient = get_dislocation_orientation(dislocation_type=dislocation_type, crystal="fcc")
    perpend = np.cross(orient["glide_plane"], orient["dislocation_line"])
    if glide_plane == "x":
        return np.array([orient["glide_plane"], -perpend, orient["dislocation_line"]])
    else:
        return np.array([perpend, orient["glide_plane"], orient["dislocation_line"]])


def get_lattice_parameter(
    element: str,
    potential_dataframe: pd.DataFrame,
    cubic=True
) -> Annotated[float, {"units": "angstrom"}]:
    structure = bulk(element, cubic=True)
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
    lattice_parameter = fit_dict["volume_eq"]**(1 / 3)
    return lattice_parameter


def get_burgers_vector(
    lattice_parameter: Annotated[float, {"units": "angstrom"}],
    dislocation_type: str = "edge",
):
    direction = get_dislocation_orientation(dislocation_type)["burgers_vector"]
    burgers_vector = (
        lattice_parameter
        * np.asarray(direction)
        / np.linalg.norm(direction)
        / np.sqrt(2)
    )
    return burgers_vector


def _as_quantity(value, unit, ureg: UnitRegistry):
    return value if hasattr(value, "to") else value * unit


def get_elastic_matrix(
    fit_dict: dict,
) -> Annotated[
    np.ndarray, {"shape": (6, 6), "units": "gigapascal", "label": "elastic_matrix"}
]:
    return fit_dict["elastic_matrix"]


def evaluate_lammps_for_elastic_matrix(
    structure: Atoms,
    potential_name="1995--Angelo-J-E--Ni-Al-H--LAMMPS--ipr1",
    num_point=5,
    eps_range=0.005,
) -> Annotated[np.ndarray, {"shape": (6, 6), "units": "gigapascal"}]:
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


@fr.workflow
def get_elastic_tensor(
    element="Ni",
    cubic=True,
    potential_name="1995--Angelo-J-E--Ni-Al-H--LAMMPS--ipr1",
    num_point=5,
    eps_range=0.005,
) -> Annotated[
    np.ndarray, {"shape": (6, 6), "units": "gigapascal", "label": "elastic_matrix"}
]:
    structure = bulk(element, cubic=cubic)
    elastic_matrix = evaluate_lammps_for_elastic_matrix(
        structure=structure,
        potential_name=potential_name,
        num_point=num_point,
        eps_range=eps_range,
    )   
    return elastic_matrix


def get_partial_burgers_vectors(
    burgers_vector: Annotated[np.ndarray, {"shape": (3,), "units": "angstrom"}],
    orientation: np.ndarray | list,
) -> Annotated[np.ndarray, {"shape": (2, 3), "units": "angstrom"}]:
    shockley_partials = get_shockley_partials(burgers_vector)
    burgers_vectors = tools.crystal_to_box(shockley_partials, orientation=orientation)
    return burgers_vectors


def get_dislocation_distance(
    medium,
    burgers_vectors: Annotated[np.ndarray, {"shape": (2, 3), "units": "angstrom"}],
    x_min: Annotated[float, {"units": "angstrom"}] = -10,
    x_max: Annotated[float, {"units": "angstrom"}] = 10,
    n_x: int = 100,
    sfe: Annotated[float, {"units": "millijoule / meter**2"}] = 90,
) -> Annotated[float, {"units": "angstrom", "label": "dislocation_distance"}]:
    ureg = UnitRegistry()
    x = _as_quantity(np.linspace(x_min, x_max, n_x)[:, None] * [1, 0], ureg.angstrom, ureg)
    stress = _as_quantity(
        medium.get_dislocation_stress(
            x,
            _as_quantity(burgers_vectors[0], ureg.angstrom, ureg),
        ),
        ureg.gigapascal,
        ureg,
    )
    f = _as_quantity(
        medium.get_dislocation_force(
            stress,
            glide_plane=[0, 1, 0],
            burgers_vector=_as_quantity(burgers_vectors[1], ureg.angstrom, ureg),
        ),
        ureg.millijoule / ureg.meter**2,
        ureg,
    ).to("millijoule/meter**2").magnitude
    F = (f - sfe)[:, 0]
    X = x.magnitude[:, 0]
    return X[np.abs(F).argmin()]


def get_hydrogen_structure(element="Ni", n_repeat=3):
    bulk_structure = bulk(element, cubic=True)
    structure = bulk_structure.repeat(n_repeat) + Atom(
        symbol="H", position=[0, 0, 0.5 * bulk_structure.cell[0, 0]]
    )
    return structure


def get_dipole_tensor(
    structure, potential_dataframe
) -> Annotated[np.ndarray, {"shape": (3, 3), "units": "eV"}]:
    relaxed_structure = optimize_positions_with_lammpslib(
        structure, potential_dataframe=potential_dataframe
    )
    result = calc_static_with_lammpslib(
        relaxed_structure, potential_dataframe=potential_dataframe
    )
    ureg = UnitRegistry()
    dipole_tensor = (
        (
            _as_quantity(result["stress"], ureg.bar, ureg)
            * _as_quantity(result["volume"], ureg.angstrom**3, ureg)
        ).to("eV").magnitude
    )
    return dipole_tensor


def linspace(
    x_min: Annotated[float, {"units": "angstrom"}],
    x_max: Annotated[float, {"units": "angstrom"}],
    n_points: int,
) -> Annotated[np.ndarray, {"units": "angstrom"}]:
    return np.linspace(x_min, x_max, n_points)


def create_mesh(
    x: Annotated[np.ndarray, {"units": "angstrom"}],
) -> Annotated[np.ndarray, {"units": "angstrom"}]:
    mesh = np.meshgrid(x, x, indexing="ij")
    mesh = np.stack(mesh, axis=-1).reshape(-1, 2)
    return mesh


def get_strain_field(
    medium,
    mesh: Annotated[np.ndarray, {"units": "angstrom"}],
    d_dislocations: Annotated[np.ndarray, {"units": "angstrom"}],
    burgers_vectors: Annotated[np.ndarray, {"shape": (2, 3), "units": "angstrom"}],
) -> np.ndarray:
    strain = medium.get_dislocation_strain(
        mesh - np.array([0.5, 0]) * d_dislocations, burgers_vector=burgers_vectors[0]
    )
    strain += medium.get_dislocation_strain(
        mesh + np.array([0.5, 0]) * d_dislocations, burgers_vector=burgers_vectors[1]
    )
    return strain


def get_binding_energy_field(
    dipole_tensor: Annotated[np.ndarray, {"shape": (3, 3), "units": "eV"}],
    strain: np.ndarray,
):
    binding_energy = -(dipole_tensor * strain).sum(axis=(-1, -2))
    n_x = int(np.sqrt(len(binding_energy)))
    binding_energy = binding_energy.reshape(n_x, n_x)
    return binding_energy


def get_medium(
    elastic_matrix: Annotated[
        np.ndarray, {"shape": (6, 6), "units": "gigapascal", "label": "elastic_matrix"}
    ],
    orientation: list | np.ndarray | None,
):
    medium = LinearElasticity(
        C_tensor=elastic_matrix, orientation=np.asarray(orientation)
    )
    return medium


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
    structure = get_hydrogen_structure(element=element, n_repeat=n_repeat)
    orientation = get_orientation(dislocation_type=dislocation_type)
    potential_dataframe = get_potential_by_name(potential_name=potential_name)
    dipole_tensor = get_dipole_tensor(structure, potential_dataframe)
    elastic_matrix = get_elastic_tensor(
        element=element, cubic=cubic, potential_name=potential_name
    )
    medium = get_medium(elastic_matrix, orientation=orientation)
    lattice_parameter = get_lattice_parameter(element, potential_dataframe, cubic=cubic)
    burgers_vector = get_burgers_vector(lattice_parameter, dislocation_type)
    burgers_vectors = get_partial_burgers_vectors(
        burgers_vector, orientation=orientation
    )
    d_dislocations = get_dislocation_distance(
        medium, burgers_vectors, x_min, x_max, n_x
    )
    x = linspace(x_min, x_max, n_x)
    mesh = create_mesh(x)
    strain_field = get_strain_field(medium, mesh, d_dislocations, burgers_vectors)
    binding_energy_field = get_binding_energy_field(dipole_tensor, strain_field)
    return x, binding_energy_field
