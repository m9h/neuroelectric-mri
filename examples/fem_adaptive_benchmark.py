"""Release the adaptive conservative-current benchmark and its tagged mesh."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

import numpy as np
from dolfinx import fem
from dolfinx.io import XDMFFile
from mpi4py import MPI

if __package__:
    from .fem_tank import build_tank_mesh
    from .fem_adaptive import solve_adaptive
    from .fem_benchmark import _sample
else:
    from fem_tank import build_tank_mesh
    from fem_adaptive import solve_adaptive
    from fem_benchmark import _sample


def _plot(root, history, fields):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), layout="constrained")
    fig.suptitle("Adaptive conservative current · DOLFINx / RT · synthetic", weight="bold")
    counts = [row["cells"] for row in history]
    lower = np.array([row["primal_power_w"] for row in history]) / 1e-6
    upper = np.array([row["complementary_power_w"] for row in history]) / 1e-6
    axes[0].fill_between(counts, lower, upper, color="#7c3aed", alpha=.18)
    axes[0].plot(counts, lower, "o-", color="#7c3aed", label="Primal lower")
    axes[0].plot(counts, upper, "o-", color="#0891b2", label="Complementary upper")
    axes[0].set(xlabel="Tetrahedra", ylabel="Resistance (Ω)", title="Energy interval · 1 mA")
    axes[0].grid(alpha=.2)
    axes[0].legend()
    x, z = fields["slice_x_m"] * 1e3, fields["slice_z_m"] * 1e3
    current = np.linalg.norm(fields["slice_j_a_m2"], axis=-1)
    artist = axes[1].pcolormesh(x, z, current.T, shading="nearest", cmap="viridis")
    fig.colorbar(artist, ax=axes[1], label="A/m²", shrink=.85)
    indicator = fields["slice_indicator_w_m3"]
    artist = axes[2].pcolormesh(x, z, indicator.T, shading="nearest", cmap="magma",
                               norm=LogNorm(vmin=max(float(indicator.min()), 1e-12), vmax=float(indicator.max())))
    fig.colorbar(artist, ax=axes[2], label="Indicator / cell volume (W/m³)", shrink=.85)
    axes[1].set_title("Equilibrated current")
    axes[2].set_title("Where refinement is needed")
    for ax in axes[1:]:
        ax.set(xlabel="x (mm)", ylabel="z (mm)", aspect="equal")
    fig.savefig(root / "summary.png", dpi=160)
    plt.close(fig)


def write_benchmark(output, *, max_refinements=2, plot=False):
    """Write fields, refinement history, mesh, stopping status and provenance."""
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    data = build_tank_mesh(mesh_size_m=.03, patch_fraction=.5, mesh_path=root / "initial.msh")
    result = solve_adaptive(data, .5, {101: .001, 102: -.001}, {101: .002, 102: .005},
                            max_refinements=max_refinements)
    data, solution, flux = result.mesh_data, result.solution, result.reconstruction
    domain = data.mesh
    domain.name, data.cell_tags.name, data.facet_tags.name = "mesh", "materials", "boundaries"
    with XDMFFile(domain.comm, str(root / "mesh.xdmf"), "w") as file:
        file.write_mesh(domain)
        file.write_meshtags(data.cell_tags, domain.geometry)
        file.write_meshtags(data.facet_tags, domain.geometry)
    count = domain.topology.index_map(3).size_local
    cell_ids = np.arange(count, dtype=np.int32)
    vertices = domain.geometry.x[domain.geometry.dofmaps[0]]
    centers = vertices.mean(axis=1)
    volumes = abs(np.linalg.det(vertices[:, 1:] - vertices[:, :1])) / 6
    density_space = fem.functionspace(domain, ("DG", 0))
    density = fem.Function(density_space)
    density.x.array[[density_space.dofmap.cell_dofs(c)[0] for c in cell_ids]] = flux.cell_indicator_squared_w / volumes
    domain.topology.create_cell_permutations()
    flux_space = flux.current_density_a_m2.function_space
    x, z = np.linspace(.0001, .1199, 101), np.linspace(.0001, .0599, 51)
    xx, zz = np.meshgrid(x, z, indexing="ij")
    points = np.column_stack([xx.ravel(), np.full(xx.size, .04), zz.ravel()])
    fields = dict(
        geometry_m=domain.geometry.x, geometry_cell_dofs=domain.geometry.dofmaps[0],
        cell_centers_m=centers, cell_volume_m3=volumes,
        cell_tag_indices=data.cell_tags.indices, cell_tags=data.cell_tags.values,
        potential_v=solution.potential_v.x.array,
        potential_cell_dofs=np.array([solution.potential_v.function_space.dofmap.cell_dofs(c) for c in cell_ids]),
        rt_coefficients=flux.current_density_a_m2.x.array,
        rt_cell_dofs=np.array([flux_space.dofmap.cell_dofs(c) for c in cell_ids]),
        cell_permutations=domain.topology.get_cell_permutation_info(),
        cell_center_j_a_m2=flux.current_density_a_m2.eval(centers, cell_ids),
        cell_center_raw_j_a_m2=solution.current_density_a_m2.eval(centers, cell_ids),
        cell_current_imbalance_a=flux.cell_current_imbalance_a,
        cell_indicator_squared_w=flux.cell_indicator_squared_w,
        slice_x_m=x, slice_z_m=z,
        slice_j_a_m2=_sample(flux.current_density_a_m2, points).reshape(*xx.shape, 3),
        slice_indicator_w_m3=_sample(density, points).reshape(xx.shape),
    )
    np.savez_compressed(root / "fields.npz", **fields)
    with (root / "history.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=result.history[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(result.history)
    parameters = dict(
        synthetic=True, lengths_m=[.12, .08, .06], conductivity_s_m=.5,
        electrode_currents_a={101: .001, 102: -.001}, contact_impedance_ohm_m2={101: .002, 102: .005},
        initial_mesh_size_m=.03, patch_fraction=.5, marking_fraction=.5,
        relative_gap_tolerance=.01, max_refinements=max_refinements,
        converged=result.converged, stop_reason="tolerance" if result.converged else "refinement_limit",
        achieved_relative_energy_gap=result.history[-1]["relative_energy_gap"],
        potential_family="Lagrange", potential_degree=2, current_family="RT", current_basix_degree=2,
        current_mapping="contravariant Piola; coefficients are basis moments, not point values",
        coefficient_order="DOLFINx local ordering; XDMF reload may renumber cells and DOFs; re-solve for restart",
        bound_scope="CEM energy / two-electrode resistance on the stated geometry, up to measured numerical residuals; not a pointwise J or MRI phase error bound",
        slice_y_m=.04,
    )
    (root / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    generated = ["initial.msh", "mesh.xdmf", "mesh.h5", "fields.npz", "history.csv", "parameters.json"]
    if plot:
        _plot(root, result.history, fields)
        generated.append("summary.png")
    repository = Path(__file__).resolve().parents[1]
    sources = [f"examples/{name}.py" for name in
               ("fem_tank", "fem_flux", "fem_adaptive", "fem_benchmark", "fem_adaptive_benchmark")]
    manifest = dict(synthetic=True, schema_version=1, python=platform.python_version(), platform=platform.platform(),
                    mpi_library=MPI.Get_library_version(),
                    versions={p: importlib.metadata.version(p) for p in
                              ("fenics-dolfinx", "fenics-basix", "fenics-ffcx", "fenics-ufl", "numpy", "scipy", "gmsh", "mpi4py")},
                    source_sha256={p: hashlib.sha256((repository / p).read_bytes()).hexdigest() for p in sources},
                    files_sha256={p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in generated})
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return result.history


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--max-refinements", type=int, default=2)
    args = parser.parse_args()
    history = write_benchmark(args.output, max_refinements=args.max_refinements, plot=args.plot)
    print(f"Wrote {len(history)} adaptive levels to {args.output}; converged={history[-1]['converged']}")
