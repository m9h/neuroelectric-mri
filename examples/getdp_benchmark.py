"""Same-mesh GetDP/DOLFINx CEM comparison and reproducible synthetic release."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
from types import SimpleNamespace

import numpy as np
from dolfinx import mesh as dmesh
from dolfinx.io import XDMFFile, gmsh as gmshio
from mpi4py import MPI

if __package__:
    from .fem_benchmark import _sample
    from .fem_tank import solve_cem
    from .fem_flux import reconstruct_current
    from .getdp_cem import find_getdp, solve_getdp
else:
    from fem_benchmark import _sample
    from fem_tank import solve_cem
    from fem_flux import reconstruct_current
    from getdp_cem import find_getdp, solve_getdp

ROOT = Path(__file__).resolve().parents[1]
CURRENTS = {101: .001, 102: -.001}
CONTACTS = {101: .002, 102: .005}


def write_mesh(data, path):
    """Export serial linear tetrahedra and physical labels as ASCII Gmsh 2.2.

    Preserve coordinates at round-trip double precision. Map topology entities
    to geometric DOFs explicitly; local topology IDs need not be geometry IDs.
    Reorient tetrahedra to positive determinants without changing connectivity.
    Triangle orientation is immaterial to the scalar CEM contact integrals.
    """
    domain = data.mesh
    if domain.comm.size != 1 or domain.topology.cell_name() != "tetrahedron":
        raise ValueError("Export requires a serial tetrahedral mesh")
    n = domain.topology.index_map(3).size_local
    if domain.geometry.dofmaps[0].shape != (n, 4):
        raise ValueError("Export requires linear geometry")
    if not np.array_equal(data.cell_tags.indices, np.arange(n)):
        raise ValueError("Every tetrahedron needs exactly one physical material tag")
    if np.any(data.cell_tags.values <= 0) or np.any(data.facet_tags.values <= 0):
        raise ValueError("Physical tags must be positive")
    xyz = domain.geometry.x
    if not np.all(np.isfinite(xyz)):
        raise ValueError("Geometry must be finite")
    tets = dmesh.entities_to_geometry(domain, 3, data.cell_tags.indices).copy()
    triangles = dmesh.entities_to_geometry(domain, 2, data.facet_tags.indices)
    vertices = xyz[tets]
    determinants = np.linalg.det(vertices[:, 1:] - vertices[:, :1])
    if np.any(determinants == 0):
        raise ValueError("Degenerate tetrahedron")
    negative = determinants < 0
    tets[negative, 1], tets[negative, 2] = tets[negative, 2].copy(), tets[negative, 1].copy()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write(f"$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n{len(xyz)}\n")
        for node, point in enumerate(xyz, 1):
            stream.write(f"{node} " + " ".join(format(float(x), ".17g") for x in point) + "\n")
        stream.write(f"$EndNodes\n$Elements\n{len(tets) + len(triangles)}\n")
        element = 0
        for kind, entities, tags in ((2, triangles, data.facet_tags.values), (4, tets, data.cell_tags.values)):
            for entity, tag in zip(entities, tags):
                element += 1
                stream.write(f"{element} {kind} 2 {tag} {tag} " + " ".join(str(int(v) + 1) for v in entity) + "\n")
        stream.write("$EndElements\n")


def _relative(actual, reference):
    return float(np.linalg.norm(np.asarray(actual) - reference) / np.linalg.norm(reference))


def compare_solvers(mesh_path, sigma, *, degree=2, output_dir=None):
    """Compare each solver's own fields at the same tetrahedral barycenters.

    Relative RMS is an unweighted sample norm, not a volume L2 error or a
    discretization-error bound. GetDP and DOLFINx assemble and integrate
    independently, sharing only geometry, physical parameters and conventions.
    """
    data = gmshio.read_from_msh(str(mesh_path), MPI.COMM_SELF)
    fem = solve_cem(data, sigma, CURRENTS, CONTACTS, degree=degree)
    getdp = solve_getdp(mesh_path, sigma, CURRENTS, CONTACTS, degree=degree, output_dir=output_dir)
    points = getdp.sample_coordinates_m
    potential = _sample(fem.potential_v, points).ravel()
    electric = _sample(fem.electric_field_v_m, points)
    current = _sample(fem.current_density_a_m2, points)
    row = dict(
        degree=degree, cells=len(points),
        dolfinx_voltage_v=fem.electrode_potential_v[101], getdp_voltage_v=getdp.electrode_potential_v[101],
        voltage_relative_difference=_relative(getdp.electrode_potential_v[101], fem.electrode_potential_v[101]),
        potential_sample_relative_rms=_relative(getdp.potential_v, potential),
        electric_sample_relative_rms=_relative(getdp.electric_field_v_m, electric),
        current_sample_relative_rms=_relative(getdp.current_density_a_m2, current),
        getdp_bulk_power_w=getdp.bulk_power_w, getdp_contact_power_w=getdp.contact_power_w,
        getdp_supplied_power_w=getdp.supplied_power_w,
        bulk_power_relative_difference=_relative(getdp.bulk_power_w, fem.bulk_power_w),
        contact_power_relative_difference=_relative(getdp.contact_power_w, fem.contact_power_w),
        getdp_max_current_error_a=max(abs(getdp.electrode_current_a[t] - i) for t, i in CURRENTS.items()),
        getdp_relative_power_error=abs((getdp.bulk_power_w + getdp.contact_power_w) / getdp.supplied_power_w - 1),
        resistance_lower_ohm="", resistance_upper_ohm="",
    )
    if output_dir is not None:
        np.savez_compressed(Path(output_dir) / "comparison.npz", coordinates_m=points,
                            getdp_u_v=getdp.potential_v, dolfinx_u_v=potential,
                            getdp_e_v_m=getdp.electric_field_v_m, dolfinx_e_v_m=electric,
                            getdp_j_a_m2=getdp.current_density_a_m2, dolfinx_j_a_m2=current)
    return row


def _plot(root, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    fig.suptitle("Independent CEM verification · GetDP × DOLFINx · synthetic", weight="bold")
    labels = [f'{r["case"]} P{r["degree"]}' for r in rows]
    axes[0].semilogy(range(len(rows)), [r["current_sample_relative_rms"] for r in rows], "o", color="#7c3aed")
    axes[0].axhline(1e-8, color="#d97706", linestyle="--", label="test tolerance")
    axes[0].set(xticks=range(len(rows)), xticklabels=labels, ylabel="Relative J difference · sample RMS",
                title="Same geometry and polynomial degree")
    axes[0].tick_params(axis="x", labelrotation=55, labelsize=8)
    axes[0].grid(alpha=.2)
    axes[0].legend()
    patch = [r for r in rows if r["case"].startswith("finite_patch") and r["degree"] == 2]
    final = next(r for r in rows if r["case"] == "adaptive")
    axes[1].axhspan(final["resistance_lower_ohm"], final["resistance_upper_ohm"], color="#14b8a6", alpha=.2,
                    label="adaptive primal/complementary interval")
    axes[1].plot([r["cells"] for r in patch], [r["getdp_voltage_v"] / .001 for r in patch],
                 "o-", color="#7c3aed", label="GetDP P2 · independent global meshes")
    axes[1].plot(final["cells"], final["getdp_voltage_v"] / .001, "*", markersize=13,
                 color="#0891b2", label="GetDP P2 · exported adaptive mesh")
    axes[1].set(xlabel="Tetrahedra", ylabel="Resistance (Ω)", title="Refinement and energy interval", xscale="log")
    axes[1].grid(alpha=.2)
    axes[1].legend(fontsize=8)
    fig.savefig(root / "summary.png", dpi=160)
    plt.close(fig)


def write_benchmark(output, *, plot=False):
    """Archive nine comparisons, native GetDP results, field samples and hashes."""
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    mesh_dir = root / "meshes"
    mesh_dir.mkdir(exist_ok=True)
    inputs = {}
    for name in ("homogeneous", "layered", "finite_patch_0", "finite_patch_1", "finite_patch_2"):
        source = ROOT / f"data/synthetic_fem/{name}.msh"
        shutil.copyfile(source, mesh_dir / source.name)
        inputs[str(source.relative_to(ROOT))] = hashlib.sha256(source.read_bytes()).hexdigest()
    source = ROOT / "data/synthetic_fem_adaptive/mesh.xdmf"
    with XDMFFile(MPI.COMM_SELF, str(source), "r") as file:
        domain = file.read_mesh(name="mesh")
        domain.topology.create_connectivity(2, 3)
        adaptive = SimpleNamespace(mesh=domain, cell_tags=file.read_meshtags(domain, name="materials"),
                                   facet_tags=file.read_meshtags(domain, name="boundaries"))
    for path in (source, source.with_suffix(".h5")):
        inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    write_mesh(adaptive, mesh_dir / "adaptive.msh")
    tensor = np.array([[.5, .12, -.05], [.12, .8, .09], [-.05, .09, .4]])
    cases = [("homogeneous", "homogeneous", .5, 2), ("layered", "layered", {1: .2, 2: .8}, 2)]
    cases += [(f"finite_patch_{i}", f"finite_patch_{i}", .5, p) for i in range(3)
              for p in ((1, 2) if i == 0 else (2,))]
    cases += [("tensor", "finite_patch_0", tensor, p) for p in (1, 2)]
    cases += [("adaptive", "adaptive", .5, 2)]
    rows = []
    for name, mesh_name, sigma, degree in cases:
        row = compare_solvers(mesh_dir / f"{mesh_name}.msh", sigma, degree=degree,
                              output_dir=root / f"{name}_p{degree}")
        rows.append(dict(case=name, **row))
    # Recompute the interval from the reloaded mesh; do not copy a published bound.
    primal = solve_cem(adaptive, .5, CURRENTS, CONTACTS, degree=2)
    conservative = reconstruct_current(adaptive, primal, degree=2)
    rows[-1]["resistance_lower_ohm"] = primal.supplied_power_w / .001**2
    rows[-1]["resistance_upper_ohm"] = conservative.complementary_power_w / .001**2
    with (root / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    parameters = dict(synthetic=True, units="SI", current_convention="inward positive",
                      electrode_currents_a=CURRENTS, contact_impedance_ohm_m2=CONTACTS, reference_tag=102,
                      tensor_s_m=tensor.tolist(), sample_location="tetrahedral barycenters",
                      field_metric="unweighted sample RMS, not volume L2 or discretization error",
                      cases=[dict(case=name, mesh=f"meshes/{mesh_name}.msh", degree=degree,
                                  conductivity_s_m=sigma.tolist() if isinstance(sigma, np.ndarray) else sigma)
                             for name, mesh_name, sigma, degree in cases],
                      limits="Serial, linear tetrahedral geometry; no B/phase or experimental validation")
    (root / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    if plot:
        _plot(root, rows)
    executable = Path(find_getdp())
    manifest = dict(python=platform.python_version(), platform=platform.platform(),
                    versions={p: importlib.metadata.version(p) for p in ("fenics-dolfinx", "gmsh", "numpy", "scipy")},
                    getdp_executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
                    input_files_sha256=inputs,
                    code_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in (Path(__file__), *(Path(__file__).with_name(name) for name in
                                           ("getdp_cem.py", "getdp_cem.pro", "fem_tank.py", "fem_flux.py", "fem_benchmark.py")))},
                    files_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted(root.rglob("*")) if p.is_file() and p.name not in ("manifest.json", "README.md")})
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    write_benchmark(args.output, plot=args.plot)
