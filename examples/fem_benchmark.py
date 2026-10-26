"""Export the complete-electrode FEM verification release."""

import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

import numpy as np
from dolfinx import geometry

if __package__:
    from .fem_tank import build_tank_mesh, solve_cem
else:
    from fem_tank import build_tank_mesh, solve_cem


def _sample(field, points):
    domain = field.function_space.mesh
    tree = geometry.bb_tree(domain, 3)
    candidates = geometry.compute_collisions_points(tree, points)
    cells = geometry.compute_colliding_cells(domain, candidates, points)
    if any(len(cells.links(i)) == 0 for i in range(len(points))):
        raise RuntimeError("A requested slice sample lies outside the mesh")
    return field.eval(points, np.array([cells.links(i)[0] for i in range(len(points))], dtype=np.int32))


def _row(case, data, solution, h, degree, analytic_voltage=None):
    voltage = solution.electrode_potential_v[101] - solution.electrode_potential_v[102]
    return dict(case=case, mesh_size_m=h, degree=degree,
                cells=data.mesh.topology.index_map(3).size_local,
                potential_dofs=solution.potential_v.function_space.dofmap.index_map.size_local,
                voltage_v=voltage,
                analytic_voltage_v="" if analytic_voltage is None else analytic_voltage,
                voltage_relative_error="" if analytic_voltage is None else abs(voltage / analytic_voltage - 1),
                relative_residual=solution.relative_residual,
                max_electrode_current_error_a=max(abs(solution.electrode_current_a[t] - i)
                                                  for t, i in solution.applied_current_a.items()),
                bulk_power_w=solution.bulk_power_w, contact_power_w=solution.contact_power_w,
                supplied_power_w=solution.supplied_power_w,
                relative_power_error=abs((solution.bulk_power_w + solution.contact_power_w)
                                         / solution.supplied_power_w - 1))


def _plot(output, fields, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x, z = fields["slice_x_m"] * 1e3, fields["slice_z_m"] * 1e3
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), layout="constrained")
    fig.suptitle("Complete-electrode tank · DOLFINx / Gmsh · synthetic", weight="bold")
    for ax, values, title, label, cmap in (
        (axes[0], fields["slice_potential_v"] * 1e3, "Finite electrodes · potential", "mV", "magma"),
        (axes[1], np.linalg.norm(fields["slice_j_a_m2"], axis=-1), "Current concentration", "A/m²", "viridis"),
    ):
        artist = ax.pcolormesh(x, z, values.T, shading="nearest", cmap=cmap)
        fig.colorbar(artist, ax=ax, label=label, shrink=.85)
        ax.set(xlabel="x (mm)", ylabel="z (mm)", title=title, aspect="equal")
    for degree, color in ((1, "#7c3aed"), (2, "#0891b2")):
        points = [r for r in rows if r["case"] == "finite_patch" and r["degree"] == degree]
        axes[2].plot([r["potential_dofs"] for r in points], [r["voltage_v"] * 1e3 for r in points],
                     "o-", color=color, label=f"CG{degree}")
    axes[2].set(xlabel="Potential degrees of freedom", ylabel="Electrode voltage (mV)",
                xscale="log", title="h / p refinement · 1 mA")
    axes[2].grid(alpha=.2)
    axes[2].legend()
    fig.savefig(output / "summary.png", dpi=160)
    plt.close(fig)


def write_benchmark(output, mesh_sizes_m=(.03, .02, .012), plot=False):
    """Write mesh, fields, analytical comparisons and h/p refinement evidence."""
    sizes = np.asarray(mesh_sizes_m, dtype=float)
    if (sizes.ndim != 1 or len(sizes) < 2 or not np.all(np.isfinite(sizes))
            or np.any(sizes <= 0) or np.any(np.diff(sizes) >= 0)):
        raise ValueError("Provide at least two strictly decreasing positive finite mesh sizes")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    contacts = {101: .002, 102: .005}
    currents = {101: .001, 102: -.001}
    lengths = [.12, .08, .06]
    area = lengths[1] * lengths[2]
    rows, generated = [], []
    for layered in (False, True):
        case = "layered" if layered else "homogeneous"
        mesh_file = f"{case}.msh"
        data = build_tank_mesh(layered=layered, mesh_path=output / mesh_file)
        generated.append(mesh_file)
        sigma = {1: .2, 2: .8} if layered else .5
        resistance_area = .06 / .2 + .06 / .8 if layered else .12 / .5
        exact = .001 * (resistance_area + sum(contacts.values())) / area
        for degree in (1, 2):
            out = solve_cem(data, sigma, currents, contacts, degree=degree)
            rows.append(_row(case, data, out, .02, degree, exact))
    for index, h in enumerate(sizes):
        mesh_file = f"finite_patch_{index}.msh"
        data = build_tank_mesh(mesh_size_m=float(h), patch_fraction=.5, mesh_path=output / mesh_file)
        generated.append(mesh_file)
        for degree in (1, 2):
            out = solve_cem(data, .5, currents, contacts, degree=degree)
            rows.append(_row("finite_patch", data, out, float(h), degree))
    # Archive coefficients and local block-dof maps for the finest CG2 solution.
    reverse = solve_cem(data, .5, {t: -i for t, i in currents.items()}, contacts, degree=2)
    off = solve_cem(data, .5, dict.fromkeys(currents, 0.), contacts, degree=2)
    p_space = out.potential_v.function_space
    j_space = out.current_density_a_m2.function_space
    n = data.mesh.topology.index_map(3).size_local
    x = np.linspace(.0001, .1199, 101)
    z = np.linspace(.0001, .0599, 51)
    xx, zz = np.meshgrid(x, z, indexing="ij")
    points = np.column_stack([xx.ravel(), np.full(xx.size, .04), zz.ravel()])
    fields = dict(
        potential_dof_coordinates_m=p_space.tabulate_dof_coordinates(),
        potential_cell_dofs=np.array([p_space.dofmap.cell_dofs(c) for c in range(n)]),
        potential_v=out.potential_v.x.array,
        potential_reversed_v=reverse.potential_v.x.array,
        potential_off_v=off.potential_v.x.array,
        vector_dof_coordinates_m=j_space.tabulate_dof_coordinates(),
        vector_cell_dofs=np.array([j_space.dofmap.cell_dofs(c) for c in range(n)]),
        electric_field_v_m=out.electric_field_v_m.x.array.reshape(-1, 3),
        current_density_a_m2=out.current_density_a_m2.x.array.reshape(-1, 3),
        cell_indices=data.cell_tags.indices, cell_tags=data.cell_tags.values,
        mesh_geometry_m=data.mesh.geometry.x, mesh_geometry_dofmap=data.mesh.geometry.dofmaps[0],
        slice_x_m=x, slice_z_m=z,
        slice_potential_v=_sample(out.potential_v, points).reshape(xx.shape),
        slice_j_a_m2=_sample(out.current_density_a_m2, points).reshape(*xx.shape, 3),
    )
    np.savez_compressed(output / "fields.npz", **fields)
    with (output / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    parameters = dict(lengths_m=lengths, mesh_sizes_m=sizes.tolist(),
                      patch_fraction=.5, electrode_edge_mesh_size_factor=.125,
                      electrode_edge_grading_distance="mesh_size_m", electrode_currents_a=currents,
                      contact_impedance_ohm_m2=contacts, reference_electrode=102,
                      conductivity_s_m={"homogeneous": .5, "layered": {1: .2, 2: .8}, "finite_patch": .5},
                      archived_field_mesh=f"finite_patch_{len(sizes) - 1}.msh", archived_field_degree=2,
                      archived_cell_order="DOLFINx; use archived dof maps, not Gmsh element ordering",
                      current_sign="positive into tank", contact_units="ohm m^2",
                      slice_y_m=.04, current_field="DG1 gradient reconstruction; not H(div) equilibrated",
                      convergence_reference="No exact finite-patch solution; report successive h/p voltage changes")
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    generated += ["fields.npz", "metrics.csv", "parameters.json"]
    if plot:
        _plot(output, fields, rows)
        generated.append("summary.png")
    root = Path(__file__).resolve().parents[1]
    sources = ["examples/fem_tank.py", "examples/fem_benchmark.py"]
    manifest = dict(synthetic=True, solver="DOLFINx CEM / SciPy serial", schema_version=1,
                    python=platform.python_version(), platform=platform.platform(),
                    versions={p: importlib.metadata.version(p) for p in
                              ("fenics-dolfinx", "fenics-basix", "fenics-ffcx", "fenics-ufl",
                               "numpy", "scipy", "gmsh", "mpi4py")},
                    source_sha256={p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in sources},
                    files_sha256={p: hashlib.sha256((output / p).read_bytes()).hexdigest() for p in generated})
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    rows = write_benchmark(args.output, plot=args.plot)
    print(f"Wrote {len(rows)} FEM cases/refinement levels to {args.output}")
