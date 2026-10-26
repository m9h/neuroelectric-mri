"""DOLFINx complete-electrode reference on tagged tetrahedra, in SI units.

CG1/CG2 potential, piecewise-constant symmetric positive-definite conductivity,
and one scalar unknown per equipotential electrode. The small reference system
is assembled by DOLFINx/UFL and solved with SciPy in serial. No PETSc/MPI scaling
or locally conservative H(div) current reconstruction is claimed here.
"""

from collections.abc import Mapping
from dataclasses import dataclass

import gmsh
import numpy as np
from mpi4py import MPI
from dolfinx import fem, mesh
from dolfinx.io import gmsh as gmshio
import ufl
from scipy.sparse import bmat, csr_matrix, diags
from scipy.sparse.linalg import spsolve


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
        raise ValueError(f"{name} must be a scalar")
    return float(value)


def build_tank_mesh(lengths_m=(.12, .08, .06), mesh_size_m=.02,
                    patch_fraction=1., layered=False, mesh_path=None):
    """Build a conforming OpenCASCADE/Gmsh tank on MPI.COMM_SELF.

    Parameters are in metres. ``patch_fraction`` scales BOTH transverse
    electrode dimensions; area scales as its square. Electrodes 101/102 are
    centered on x-/x+, remaining exterior facets have tag 199. Material 1
    fills the tank unless ``layered=True``, when the right half has tag 2.
    Optional ``mesh_path`` writes the tagged MSH 4.1 mesh, before import.
    Finite electrode edges use a distance size field: h/8 at the edge, grading
    to h over distance h. Linear geometry is exact for this planar CAD; solution order
    is chosen separately in ``solve_cem``. The builder owns its Gmsh session.
    """
    lengths = _real(lengths_m, "lengths_m")
    h = _scalar(mesh_size_m, "mesh_size_m")
    fraction = _scalar(patch_fraction, "patch_fraction")
    if lengths.shape != (3,) or np.any(lengths <= 0) or h <= 0 or not 0 < fraction <= 1:
        raise ValueError("Need three positive lengths, positive mesh size, and 0 < patch_fraction <= 1")
    if gmsh.isInitialized():
        raise RuntimeError("Close the existing Gmsh session before using this standalone builder")
    lx, ly, lz = lengths
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.model.add("complete-electrode-tank")
        occ = gmsh.model.occ
        if layered:
            volumes = [(3, occ.addBox(0, 0, 0, lx / 2, ly, lz)),
                       (3, occ.addBox(lx / 2, 0, 0, lx / 2, ly, lz))]
        else:
            volumes = [(3, occ.addBox(0, 0, 0, lx, ly, lz))]
        patches = []
        if fraction < 1:
            y0, y1 = ly * (1 - fraction) / 2, ly * (1 + fraction) / 2
            z0, z1 = lz * (1 - fraction) / 2, lz * (1 + fraction) / 2
            for x in (0., lx):
                points = [occ.addPoint(x, y, z) for y, z in
                          ((y0, z0), (y1, z0), (y1, z1), (y0, z1))]
                lines = [occ.addLine(points[i], points[(i + 1) % 4]) for i in range(4)]
                patches.append((2, occ.addPlaneSurface([occ.addCurveLoop(lines)])))
        patch_entities = {}
        if len(volumes) > 1 or patches:
            # Boolean fragments imprint electrode edges and material interfaces.
            _, fragment_map = occ.fragment(volumes[:1], volumes[1:] + patches)
            # Preserve CAD ancestry: a patch and its surrounding face can
            # have identical centroids and even identical areas.
            for i in range(len(patches)):
                patch_entities[101 + i] = {t for d, t in fragment_map[len(volumes) + i] if d == 2}
        occ.synchronize()
        materials = {1: [], 2: []}
        for _, tag in gmsh.model.getEntities(3):
            material = 2 if layered and occ.getCenterOfMass(3, tag)[0] > lx / 2 else 1
            materials[material].append(tag)
        for tag, entities in materials.items():
            if entities:
                gmsh.model.addPhysicalGroup(3, entities, tag, name=f"material_{tag}")
        groups = {101: [], 102: [], 199: []}
        expected_area = ly * lz * fraction**2
        exterior = gmsh.model.getBoundary(gmsh.model.getEntities(3), combined=True, oriented=False)
        for dim, tag in exterior:
            center = occ.getCenterOfMass(dim, tag)
            area = occ.getMass(dim, tag)
            electrode = None
            for x, label in ((0., 101), (lx, 102)):
                if fraction < 1:
                    if tag in patch_entities[label]:
                        electrode = label
                elif (np.allclose(center, [x, ly / 2, lz / 2], rtol=0, atol=1e-10)
                        and np.isclose(area, expected_area, rtol=1e-8, atol=0)):
                    electrode = label
            groups[199 if electrode is None else electrode].append(tag)
        if any(not groups[t] for t in (101, 102)):
            raise RuntimeError("CAD electrode surface identification failed")
        for tag, entities in groups.items():
            if entities:
                gmsh.model.addPhysicalGroup(2, entities, tag, name=f"boundary_{tag}")
        gmsh.option.setNumber("Mesh.MeshSizeMin", h / 8 if fraction < 1 else h)
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        if fraction < 1:
            # Resolve the contact/insulation transition rather than spending
            # the same fine mesh throughout the otherwise smooth tank.
            electrode_surfaces = [(2, t) for label in (101, 102) for t in groups[label]]
            edges = gmsh.model.getBoundary(electrode_surfaces, combined=False, oriented=False)
            field = gmsh.model.mesh.field
            distance = field.add("Distance")
            field.setNumbers(distance, "CurvesList", sorted({t for _, t in edges}))
            field.setNumber(distance, "Sampling", 100)
            threshold = field.add("Threshold")
            for key, value in (("InField", distance), ("SizeMin", h / 8),
                               ("SizeMax", h), ("DistMin", 0.), ("DistMax", h)):
                field.setNumber(threshold, key, value)
            field.setAsBackgroundMesh(threshold)
            gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(3)
        if mesh_path is not None:
            gmsh.write(str(mesh_path))
        return gmshio.model_to_mesh(gmsh.model, MPI.COMM_SELF, rank=0)
    finally:
        gmsh.finalize()


@dataclass
class CEMSolution:
    """SI fields and weak-form electrode diagnostics; currents are inward-positive.

    E/J are discontinuous polynomial fields from the CG potential gradient.
    They are NOT an equilibrated, locally conservative H(div) flux. Electrode
    current is integrated from the contact law (U-u)/z, not raw grad(u) traces.
    """

    potential_v: fem.Function
    electric_field_v_m: fem.Function
    current_density_a_m2: fem.Function
    conductivity_s_m: fem.Function
    applied_current_a: dict[int, float]
    electrode_potential_v: dict[int, float]
    electrode_current_a: dict[int, float]
    electrode_area_m2: dict[int, float]
    bulk_power_w: float
    contact_power_w: float
    supplied_power_w: float
    relative_residual: float


def _tensor(value):
    value = _real(value, "conductivity_s_m")
    if value.shape == ():
        value = np.eye(3) * value
    if value.shape != (3, 3) or not np.allclose(value, value.T, rtol=1e-12, atol=0):
        raise ValueError("Conductivity must be a scalar or symmetric 3x3 tensor")
    if np.min(np.linalg.eigvalsh(value)) <= 0:
        raise ValueError("Conductivity must be positive definite")
    return value


def solve_cem(mesh_data, conductivity_s_m, electrode_currents_a,
              contact_impedance_ohm_m2, *, degree=1, reference_tag=None):
    """Solve div(sigma grad(u))=0 with complete-electrode boundary conditions.

    ``mesh_data`` has a 3D tetrahedral mesh, cell tags and exterior facet tags.
    Conductivity is a scalar/tensor or a mapping from EVERY material tag to a
    scalar/tensor, in S/m. Electrode/current and electrode/contact mappings must
    have identical keys. Current is positive INTO the tank; its sum must be
    zero. Contact impedance is strictly positive in ohm m^2. Unlisted boundary
    tags are insulated. A listed zero-current electrode floats equipotentially.

    The reference electrode (default largest tag) is fixed at 0 V. Its current
    equation is redundant and checked in the full residual after the solve.
    The singular z=0 ideal-contact limit is not supported: use a different
    constrained formulation instead of silently replacing zero by epsilon.
    This reference supports CG1/CG2 and serial execution only.
    """
    domain = mesh_data.mesh
    if domain.comm.size != 1:
        raise ValueError("This SciPy reference is serial; use a one-rank communicator")
    if domain.topology.dim != 3 or domain.geometry.dim != 3:
        raise ValueError("A three-dimensional mesh is required")
    if degree not in (1, 2) or not isinstance(degree, (int, np.integer)):
        raise ValueError("degree must be 1 or 2")
    if not isinstance(electrode_currents_a, Mapping) or len(electrode_currents_a) < 2:
        raise ValueError("At least two electrode currents are required")
    if (not isinstance(contact_impedance_ohm_m2, Mapping)
            or set(electrode_currents_a) != set(contact_impedance_ohm_m2)):
        raise ValueError("Current and contact-impedance electrode tags must match")
    tags = sorted(electrode_currents_a)
    applied = {t: _scalar(electrode_currents_a[t], "current_a") for t in tags}
    contacts = {t: _scalar(contact_impedance_ohm_m2[t], "contact_impedance_ohm_m2") for t in tags}
    if any(v <= 0 for v in contacts.values()):
        raise ValueError("Contact impedance must be strictly positive")
    current_scale = sum(abs(v) for v in applied.values())
    if abs(sum(applied.values())) > 1e-12 * current_scale:
        raise ValueError("Electrode currents must sum to zero")
    reference = tags[-1] if reference_tag is None else reference_tag
    if reference not in tags:
        raise ValueError("reference_tag must identify an electrode")
    cell_tags, facet_tags = mesh_data.cell_tags, mesh_data.facet_tags
    if cell_tags is None or facet_tags is None:
        raise ValueError("Cell and boundary facet tags are required")
    ncells = domain.topology.index_map(3).size_local
    if cell_tags.dim != 3 or not np.array_equal(np.sort(cell_tags.indices), np.arange(ncells)):
        raise ValueError("Every cell must have exactly one material tag")
    domain.topology.create_connectivity(2, 3)
    exterior = mesh.exterior_facet_indices(domain.topology)
    if facet_tags.dim != 2:
        raise ValueError("Boundary tags must identify facets")
    for tag in tags:
        facets = facet_tags.find(tag)
        if not len(facets) or not np.all(np.isin(facets, exterior)):
            raise ValueError(f"Electrode {tag} must tag nonempty exterior facets")
    material_tags = set(int(t) for t in cell_tags.values)
    if isinstance(conductivity_s_m, Mapping):
        if set(conductivity_s_m) != material_tags:
            raise ValueError("Conductivity mapping must cover exactly the material tags")
        materials = {t: _tensor(conductivity_s_m[t]) for t in material_tags}
    else:
        tensor = _tensor(conductivity_s_m)
        materials = dict.fromkeys(material_tags, tensor)
    tensor_space = fem.functionspace(domain, ("DG", 0, (3, 3)))
    sigma = fem.Function(tensor_space, name="conductivity_S_per_m")
    values = sigma.x.array.reshape(-1, 9)
    for cell, tag in zip(cell_tags.indices, cell_tags.values):
        values[tensor_space.dofmap.cell_dofs(cell)[0]] = materials[tag].ravel()
    sigma.x.scatter_forward()

    space = fem.functionspace(domain, ("Lagrange", degree))
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    dx = ufl.Measure("dx", domain=domain)
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
    bulk_form = ufl.inner(ufl.dot(sigma, ufl.grad(u)), ufl.grad(v)) * dx
    form = bulk_form
    coupling, areas = [], {}
    for tag in tags:
        z = contacts[tag]
        form += u * v / z * ds(tag)
        areas[tag] = float(fem.assemble_scalar(fem.form(1 * ds(tag))))
        if not np.isfinite(areas[tag]) or areas[tag] <= 0:
            raise ValueError(f"Electrode {tag} has no finite positive area")
        coupling.append(-fem.assemble_vector(fem.form(v / z * ds(tag))).array.copy())
    assembled = fem.assemble_matrix(fem.form(form))
    assembled.scatter_reverse()
    bulk_matrix = assembled.to_scipy()
    cross = csr_matrix(np.column_stack(coupling))
    electrode_matrix = diags([areas[t] / contacts[t] for t in tags])
    matrix = bmat([[bulk_matrix, cross], [cross.T, electrode_matrix]], format="csr")
    n = space.dofmap.index_map.size_local
    rhs = np.concatenate([np.zeros(n), [applied[t] for t in tags]])
    free = np.delete(np.arange(len(rhs)), n + tags.index(reference))
    coefficients = np.zeros_like(rhs)
    coefficients[free] = spsolve(matrix[free][:, free], rhs[free])
    if not np.all(np.isfinite(coefficients)):
        raise RuntimeError("CEM solve produced nonfinite potentials")
    residual = float(np.linalg.norm(matrix @ coefficients - rhs) / (np.linalg.norm(rhs) or 1.))
    if residual > 1e-8:
        raise RuntimeError(f"CEM residual is too large: {residual:g}")
    potential = fem.Function(space, name="potential_V")
    potential.x.array[:] = coefficients[:n]
    potential.x.scatter_forward()
    voltages = dict(zip(tags, map(float, coefficients[n:])))
    currents = {t: float(fem.assemble_scalar(fem.form((voltages[t] - potential) / contacts[t] * ds(t))))
                for t in tags}
    bulk_power = float(fem.assemble_scalar(fem.form(
        ufl.inner(ufl.dot(sigma, ufl.grad(potential)), ufl.grad(potential)) * dx)))
    contact_power = sum(float(fem.assemble_scalar(fem.form(
        (potential - voltages[t])**2 / contacts[t] * ds(t)))) for t in tags)
    vector_space = fem.functionspace(domain, ("DG", degree - 1, (3,)))
    electric = fem.Function(vector_space, name="E_V_per_m")
    current = fem.Function(vector_space, name="J_A_per_m2")
    points = vector_space.element.interpolation_points
    electric.interpolate(fem.Expression(-ufl.grad(potential), points))
    current.interpolate(fem.Expression(-ufl.dot(sigma, ufl.grad(potential)), points))
    return CEMSolution(potential, electric, current, sigma, applied, voltages, currents, areas,
                       bulk_power, contact_power, sum(applied[t] * voltages[t] for t in tags), residual)
