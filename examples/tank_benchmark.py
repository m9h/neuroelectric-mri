"""Reproduce scalar conductive-tank fields, nulls and convergence in SI units.

Analytical references: series resistance, and the harmonic polynomial
u = V0 * (X^4 - 6*X^2*Y^2 + Y^4), with X=x/L and Y=y/L.
No magnetic-field prediction is made without a specified complete circuit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import inspect
import json
from pathlib import Path
import platform

import numpy as np

if __package__:
    from .conductive_tank import solve_plate_current, solve_tank
else:
    from conductive_tank import solve_plate_current, solve_tank


def _write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _diagnostics(result):
    currents = result.boundary_current_a
    throughput = .5 * sum(abs(i) for i in currents.values())
    cell_balance = sum(np.diff(f, axis=a) for a, f in enumerate(result.face_current_a))
    return dict(
        boundary_net_current_a=float(sum(currents.values())),
        boundary_relative_imbalance=abs(sum(currents.values())) / max(throughput, 1e-30),
        max_cell_current_imbalance_a=float(np.max(abs(cell_balance))),
        boundary_power_w=result.boundary_power_w,
        dissipated_power_w=result.dissipated_power_w,
        relative_power_discrepancy=abs(result.boundary_power_w - result.dissipated_power_w)
            / max(abs(result.dissipated_power_w), 1e-30),
    )


def _harmonic_reference(x, y, z, length=.1, voltage=.01, sigma=.5):
    # z is retained to make the 3D coordinate contract explicit.
    x, y = x / length, y / length
    u = voltage * (x**4 - 6 * x*x*y*y + y**4)
    j = -sigma * voltage / length * np.stack(
        [4*x**3 - 12*x*y*y, 4*y**3 - 12*x*x*y, np.zeros_like(z)], axis=-1)
    return u, j


def _convergence(sizes):
    rows = []
    previous_error, previous_n = None, None
    for n in sizes:
        centers = (np.arange(n) + .5) * .1 / n
        coordinates = np.meshgrid(centers, centers, centers, indexing="ij")
        exact_u, exact_j = _harmonic_reference(*coordinates)
        boundaries = {}
        for axis, name in enumerate("xyz"):
            for side, endpoint in (("-", 0.), ("+", .1)):
                axes = [centers] * 3
                axes[axis] = np.array([endpoint])
                u, _ = _harmonic_reference(*np.meshgrid(*axes, indexing="ij"))
                boundaries[name + side] = np.squeeze(u, axis=axis)
        result = solve_tank(np.full((n, n, n), .5), (.1 / n,) * 3, boundaries)
        error = float(np.linalg.norm(result.potential_v - exact_u) / np.linalg.norm(exact_u))
        order = (None if previous_error is None else
                 float(np.log(previous_error / error) / np.log(n / previous_n)))
        rows.append(dict(cells_per_axis=n, n_cells=n**3, cell_size_m=.1 / n,
            potential_relative_l2_error=error, potential_observed_order=order,
            current_density_relative_l2_error=float(np.linalg.norm(result.current_density_a_m2 - exact_j)
                                                     / np.linalg.norm(exact_j)),
            **_diagnostics(result)))
        previous_error, previous_n = error, n
    return rows


def _plot_summary(output, arrays, convergence, spacing):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sigma = arrays["inclusion_conductivity_s_m"]
    z = sigma.shape[2] // 2
    x = (np.arange(sigma.shape[0]) + .5) * spacing[0] * 1000
    y = (np.arange(sigma.shape[1]) + .5) * spacing[1] * 1000
    u = arrays["inclusion_potential_v"][:, :, z].T * 1000
    j = arrays["inclusion_current_density_a_m2"][:, :, z, :]
    extent = [0, sigma.shape[0] * spacing[0] * 1000, 0, sigma.shape[1] * spacing[1] * 1000]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), layout="constrained")
    im = axes[0].imshow(u, origin="lower", extent=extent, cmap="viridis", aspect="equal")
    fig.colorbar(im, ax=axes[0], label="Potential (mV)", shrink=.8)
    axes[0].set_title("Potential: fixed total current")
    im = axes[1].imshow(np.linalg.norm(j, axis=-1).T, origin="lower", extent=extent,
                        cmap="magma", aspect="equal")
    axes[1].quiver(x, y, j[..., 0].T, j[..., 1].T, color="white", angles="xy")
    fig.colorbar(im, ax=axes[1], label="|J| (A/m²)", shrink=.8)
    axes[1].set_title("Current avoids the inclusion")
    for ax in axes[:2]:
        ax.contour(x, y, sigma[:, :, z].T, levels=[.3], colors="cyan", linewidths=1)
        ax.set(xlabel="x (mm)", ylabel="y (mm)")
    n = np.array([r["cells_per_axis"] for r in convergence])
    error = np.array([r["potential_relative_l2_error"] for r in convergence])
    axes[2].loglog(n, error, "o-", color="#6e45bc", label="Potential error")
    axes[2].loglog(n, error[-1] * (n[-1] / n)**2, "--", color="#e59629", label="Second-order slope")
    axes[2].set(xlabel="Cells per axis", ylabel="Relative L2 error",
                title="Independent harmonic reference")
    axes[2].grid(True, which="both", alpha=.2)
    axes[2].legend(frameon=False)
    fig.suptitle("Synthetic conductive tank • scalar conductivity • no scanner measurements")
    fig.savefig(output / "summary.png", dpi=180)
    plt.close(fig)


def write_benchmark(output: Path, convergence_sizes=(4, 8, 16, 32), plot: bool = False):
    """Write synthetic fields, analytical convergence and a hash manifest.

    Parameters
    ----------
    output : Path
        Output folder. Generated named files are replaced on repeat runs.
    convergence_sizes : sequence of integers
        At least two strictly increasing positive cube sizes, >=2 cells/axis.
    plot : bool
        Also write a static scientific summary (requires matplotlib).

    Returns
    -------
    dict
        The release manifest, also written as manifest.json.
    """
    sizes = tuple(convergence_sizes)
    if (len(sizes) < 2 or any(not isinstance(n, (int, np.integer)) or n < 2 for n in sizes)
            or any(b <= a for a, b in zip(sizes, sizes[1:]))):
        raise ValueError("convergence_sizes must be strictly increasing integers >=2; at least two")
    output = Path(output)
    shape, lengths, current = (8, 6, 4), np.array([.12, .08, .06]), .001
    spacing = lengths / shape
    homogeneous = np.full(shape, .5)
    layered = np.full(shape, .2)
    layered[4:] = .8
    inclusion = homogeneous.copy()
    inclusion[2:6, 1:5, 1:3] = .1
    cases = [("homogeneous", homogeneous, current), ("layered", layered, current),
             ("inclusion", inclusion, current), ("scaled_inclusion", 3*inclusion, current),
             ("reverse_polarity", inclusion, -current), ("current_off", inclusion, 0.)]
    coordinates = np.stack(np.meshgrid(*[(np.arange(n) + .5) * h
                                         for n, h in zip(shape, spacing)], indexing="ij"), axis=-1)
    arrays = {"cell_centers_m": coordinates, "spacing_m": spacing}
    metrics, case_parameters = [], {}
    for name, sigma, target_current in cases:
        result = solve_plate_current(sigma, spacing, target_current)
        for key in ("conductivity_s_m", "potential_v", "current_density_a_m2",
                    "electric_field_v_m", "divergence_a_m3"):
            arrays[name + "_" + key] = getattr(result, key)
        for axis, face_current in zip("xyz", result.face_current_a):
            arrays[f"{name}_face_current_{axis}_a"] = face_current
        voltage = float(result.boundary_potential_v["x-"].flat[0])
        row = dict(case=name, specified_current_a=target_current,
                   recovered_current_a=result.boundary_current_a["x+"],
                   required_plate_voltage_v=voltage, **_diagnostics(result))
        metrics.append(row)
        case_parameters[name] = dict(target_current_a=target_current, drive_axis=0,
            dirichlet_voltages_v={"x-": voltage, "x+": 0.},
            insulated_faces=["y-", "y+", "z-", "z+"])
    # Independent closed-form slab references, saved alongside numerical fields.
    x = coordinates[..., 0]
    jx = current / (lengths[1] * lengths[2])
    arrays["homogeneous_reference_potential_v"] = jx * (lengths[0] - x) / .5
    arrays["layered_reference_potential_v"] = jx * (
        np.maximum(.06 - x, 0) / .2 + (.12 - np.maximum(x, .06)) / .8)
    convergence = _convergence(sizes)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "fields.npz", **arrays)
    _write_csv(output / "metrics.csv", metrics)
    _write_csv(output / "convergence.csv", convergence)
    parameters = dict(
        shape=shape, lengths_m=lengths.tolist(), spacing_m=spacing.tolist(),
        coordinate_frame="right-handed Cartesian, axes x/y/z", origin_m=[0, 0, 0],
        potential_gauge="x+ electrode at 0 V", boundary_current_sign="positive outward",
        face_current_sign="positive toward increasing coordinate on every face",
        conductivity_units="S/m", current_density_units="A/m^2", electric_field_units="V/m",
        cases=case_parameters, convergence_sizes=[int(n) for n in sizes],
        convergence_reference=dict(domain_lengths_m=[.1]*3, conductivity_s_m=.5,
            potential="0.01*(X^4-6*X^2*Y^2+Y^4) V; X=x/0.1, Y=y/0.1",
            boundary="exact potential sampled at all six boundary-face centers"),
    )
    (output / "parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    files = ["fields.npz", "metrics.csv", "convergence.csv", "parameters.json"]
    if plot:
        _plot_summary(output, arrays, convergence, spacing)
        files.append("summary.png")
    manifest = dict(
        schema="neuroelectric-mri.synthetic-tank/0.1", data_kind="synthetic_numerical",
        scanner_data_available=False, measured_current_available=False,
        complete_circuit_field_available=False, conductivity_inversion_available=False,
        python_version=platform.python_version(),
        dependencies={name: importlib.metadata.version(name) for name in ("numpy", "scipy")},
        source_sha256={"generator": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       "solver": hashlib.sha256(Path(inspect.getfile(solve_tank)).read_bytes()).hexdigest()},
        files_sha256={name: hashlib.sha256((output / name).read_bytes()).hexdigest() for name in files},
        arrays={key: dict(shape=value.shape, dtype=str(value.dtype)) for key, value in arrays.items()},
        limitations=["Scalar conductivity on orthogonal Cartesian cells; no anisotropy",
                     "Ideal full-face electrodes, no contact impedance or lead fields",
                     "Cell J/E are reconstructions; face fluxes define conservation and discrete power",
                     "No independent FEM/SimNIBS comparison or acquired-data validation",
                     "Fixed-current global conductivity scaling is not identifiable from J alone"],
    )
    if plot:
        manifest["dependencies"]["matplotlib"] = importlib.metadata.version("matplotlib")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {len(cases)} synthetic tank controls and {len(sizes)} convergence grids to {output}")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grids", type=int, nargs="+", default=[4, 8, 16, 32])
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    write_benchmark(args.output, args.grids, args.plot)
