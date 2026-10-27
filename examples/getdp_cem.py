"""Independent GetDP complete-electrode reference (SI units).

NumPy/Gmsh handle validation, geometry and output only. All variational assembly,
electrode coupling, field evaluation and power integration run in GetDP.
There is no DOLFINx import or reuse of its assembled matrices in this module.
"""

from collections.abc import Mapping
from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import gmsh
import numpy as np


def find_getdp(executable=None):
    candidate = executable or os.environ.get("GETDP_EXECUTABLE") or shutil.which("getdp")
    local = Path(__file__).resolve().parents[1] / ".tools/getdp-git-MacOSARM/bin/getdp"
    if not candidate and local.is_file():
        candidate = str(local)
    if not candidate:
        raise FileNotFoundError("Install GetDP or set GETDP_EXECUTABLE")
    resolved = shutil.which(str(candidate))
    if not resolved:
        raise FileNotFoundError(f"GetDP executable not found: {candidate}")
    return str(Path(resolved).resolve())


def _real(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real")
    result = np.array(raw, dtype=float, copy=True)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def _scalar(value, name):
    value = _real(value, name)
    if value.shape != ():
        raise ValueError(f"{name} must be scalar")
    return float(value)


def _tensor(value):
    value = _real(value, "conductivity")
    if value.shape == ():
        value = np.eye(3) * value
    if value.shape != (3, 3) or not np.allclose(value, value.T, rtol=1e-12, atol=0):
        raise ValueError("Conductivity must be scalar or symmetric 3x3")
    if np.min(np.linalg.eigvalsh(value)) <= 0:
        raise ValueError("Conductivity must be positive definite")
    return value


def inspect_mesh(mesh_path):
    """Read physical tags and triangle areas from a linear tetrahedral MSH mesh."""
    if gmsh.isInitialized():
        raise RuntimeError("Close the existing Gmsh session before loading this standalone reference")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(mesh_path))
        ids, xyz, _ = gmsh.model.mesh.getNodes()
        sort = np.argsort(ids)
        ids, xyz = ids[sort], xyz.reshape(-1, 3)[sort]
        materials = set()
        areas = {}
        for dim, tag in gmsh.model.getPhysicalGroups():
            if dim not in (2, 3):
                continue
            total = 0.
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                types, _, nodes = gmsh.model.mesh.getElements(dim, int(entity))
                for kind, connectivity in zip(types, nodes):
                    if int(kind) != (2 if dim == 2 else 4):
                        raise ValueError("Only linear triangular facets and tetrahedral geometry are supported")
                    vertices = xyz[np.searchsorted(ids, connectivity)].reshape(-1, 3 if dim == 2 else 4, 3)
                    if dim == 2:
                        total += float(np.linalg.norm(np.cross(vertices[:, 1] - vertices[:, 0],
                                                               vertices[:, 2] - vertices[:, 0]), axis=1).sum() / 2)
            if dim == 3:
                materials.add(int(tag))
            else:
                areas[int(tag)] = total
        if not materials:
            raise ValueError("Mesh needs physical material regions")
        return materials, areas
    finally:
        gmsh.finalize()


@dataclass
class GetDPSolution:
    """Independent GetDP outputs, sampled at tetrahedral barycenters in SI units."""

    sample_coordinates_m: np.ndarray
    potential_v: np.ndarray
    electric_field_v_m: np.ndarray
    current_density_a_m2: np.ndarray
    applied_current_a: dict[int, float]
    electrode_potential_v: dict[int, float]
    electrode_current_a: dict[int, float]
    bulk_power_w: float
    contact_power_w: float
    supplied_power_w: float
    solver_info: str


def _number(x):
    return format(float(x), ".17g")


def _write_problem(root, materials, currents, contacts, areas, degree, reference):
    bulk = ",".join(map(str, sorted(materials)))
    electrode = ",".join(map(str, sorted(currents)))
    groups = [f"Bulk = Region[{{{bulk}}}];", f"Electrodes = Region[{{{electrode}}}];",
              "Domain = Region[{Bulk, Electrodes}];"]
    groups += [f"E{t} = Region[{t}];" for t in currents]
    functions = []
    for t, sigma in materials.items():
        functions.append(f"sigma[Region[{t}]] = Tensor[{','.join(_number(v) for v in sigma.ravel())}];")
    for t in currents:
        functions += [f"contact[E{t}] = {_number(contacts[t])};",
                      f"injection[E{t}] = {_number(currents[t])};", f"area[E{t}] = {_number(areas[t])};"]
    parameters = (f"TissueDegree = {degree};\nGroup {{\n" + "\n".join(groups) + "\n}\nFunction {\n"
                  + "\n".join(functions) + "\n}\nConstraint { { Name Gauge; Type Assign; Case {"
                  + f" {{ Region E{reference}; Value 0.; }}" + " } } }\n")
    (root / "parameters.pro").write_text(parameters)
    quantities, outputs = [], []
    for t in currents:
        for name, expression in ((f"voltage{t}", "{U} / area[]"),
                                 (f"current{t}", "({U} - {u}) / contact[]")):
            quantities.append(f"{{ Name {name}; Value {{ Integral {{ [ {expression} ]; In E{t}; "
                              "Jacobian Surface; Integration Gauss; } } }")
            outputs.append(f'Print[{name}[E{t}], OnGlobal, Format SimpleTable, File "{name}.txt"];')
    template = Path(__file__).with_suffix(".pro").read_text()
    (root / "problem.pro").write_text(template.replace("// @ELECTRODE_QUANTITIES@", "\n".join(quantities))
                                     .replace("// @ELECTRODE_OUTPUTS@", "\n".join(outputs)))


def _execute(root, mesh_path, exe):
    info = subprocess.run([exe, "-info"], capture_output=True, text=True, timeout=30, check=True)
    solver_info = "GetDP\n" + info.stdout + info.stderr
    (root / "solver-info.txt").write_text(solver_info)
    command = [exe, "problem.pro", "-msh", str(mesh_path), "-solve", "CEM", "-pos", "Export",
               "-v", "3", "-ksp_type", "preonly", "-pc_type", "lu"]
    completed = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=120)
    log = completed.stdout + completed.stderr
    (root / "solver.log").write_text(log)
    if completed.returncode:
        raise RuntimeError(f"GetDP failed with exit {completed.returncode}:\n{log[-6000:]}")
    return solver_info


def _table(path, components):
    table = np.loadtxt(path, ndmin=2)
    if table.shape[1] != 3 + components or not np.all(np.isfinite(table)):
        raise RuntimeError(f"Unexpected/nonfinite GetDP field table: {path.name}, shape {table.shape}")
    return table[:, :3], table[:, 3:]


def _global(path):
    values = np.loadtxt(path, ndmin=1)
    value = float(values.ravel()[-1])
    if not np.isfinite(value):
        raise RuntimeError(f"Nonfinite GetDP quantity: {path.name}")
    return value


def solve_getdp(mesh_path, conductivity_s_m, electrode_currents_a,
                contact_impedance_ohm_m2, *, degree=1, reference_tag=None,
                output_dir=None, executable=None):
    """Run GetDP's own CEM assembly/solve and return its fields and diagnostics.

    Input geometry is linear tetrahedral MSH with physical material/facet tags.
    Conductivity is scalar, SPD tensor, or a complete material-tag mapping.
    Currents are inward-positive and sum to zero; contact impedance is positive
    in ohm m^2. Tissue degree 1/2 uses nodal/hierarchical GetDP basis functions.
    One independent region-constant voltage unknown is used per metal electrode.
    Unlisted facets are insulated; the selected electrode gauge is zero volts.
    An output directory retains native formulation, result, fields and logs.
    """
    if isinstance(degree, (bool, np.bool_)) or not isinstance(degree, (int, np.integer)) or degree not in (1, 2):
        raise ValueError("degree must be integer 1 or 2")
    if not isinstance(electrode_currents_a, Mapping) or len(electrode_currents_a) < 2:
        raise ValueError("At least two electrode currents are required")
    tags = sorted(electrode_currents_a)
    if any(isinstance(t, bool) or not isinstance(t, (int, np.integer)) or t <= 0 for t in tags):
        raise ValueError("Electrode tags must be positive integers")
    if not isinstance(contact_impedance_ohm_m2, Mapping) or set(tags) != set(contact_impedance_ohm_m2):
        raise ValueError("Current and contact tags must match")
    currents = {int(t): _scalar(electrode_currents_a[t], "current") for t in tags}
    contacts = {int(t): _scalar(contact_impedance_ohm_m2[t], "contact") for t in tags}
    if any(z <= 0 for z in contacts.values()):
        raise ValueError("Contact impedance must be positive")
    if abs(sum(currents.values())) > 1e-12 * sum(abs(i) for i in currents.values()):
        raise ValueError("Electrode currents must sum to zero")
    reference = tags[-1] if reference_tag is None else reference_tag
    if reference not in tags:
        raise ValueError("Reference must identify an electrode")
    mesh_path = Path(mesh_path).resolve(strict=True)
    material_tags, areas = inspect_mesh(mesh_path)
    if any(t not in areas or not np.isfinite(areas[t]) or areas[t] <= 0 for t in tags):
        raise ValueError("Every electrode must have a finite positive surface area")
    if isinstance(conductivity_s_m, Mapping):
        if set(conductivity_s_m) != material_tags:
            raise ValueError("Conductivity mapping must cover exactly the material tags")
        materials = {int(t): _tensor(conductivity_s_m[t]) for t in sorted(material_tags)}
    else:
        materials = {t: _tensor(conductivity_s_m) for t in sorted(material_tags)}
    exe = find_getdp(executable)

    def run(root):
        root.mkdir(parents=True, exist_ok=True)
        _write_problem(root, materials, currents, contacts, areas, degree, reference)
        info = _execute(root, mesh_path, exe)
        coordinates, potential = _table(root / "potential.txt", 1)
        e_coordinates, electric = _table(root / "electric.txt", 3)
        j_coordinates, current = _table(root / "current.txt", 3)
        if not np.array_equal(coordinates, e_coordinates) or not np.array_equal(coordinates, j_coordinates):
            raise RuntimeError("GetDP field sample locations differ")
        voltages = {t: _global(root / f"voltage{t}.txt") for t in tags}
        measured = {t: _global(root / f"current{t}.txt") for t in tags}
        return GetDPSolution(coordinates, potential[:, 0], electric, current, currents, voltages,
                              measured, _global(root / "bulk_power.txt"), _global(root / "contact_power.txt"),
                              sum(voltages[t] * currents[t] for t in tags), info)

    if output_dir is not None:
        return run(Path(output_dir).resolve())
    with tempfile.TemporaryDirectory(prefix="getdp-cem-") as directory:
        return run(Path(directory))
