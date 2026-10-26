"""Synthetic current-phase versus concentration benchmark using mrs-jax.

The singlet basis and oracle phase correction test software separation only.
They do not validate myo-inositol quantification or an acquired MREIT result.
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
import subprocess

import numpy as np

from mrs_jax.current_encoding import CurrentEncoding
from mrs_jax.mrsi_sim import TissueModel, make_lorentzian_basis, simulate_mrsi
from current_phantom import loop_field


def fit_concentrations(fid, basis, times, t2star, b0, phase=None):
    """Real linear least squares with known basis, decay, B0 and optional phase."""
    if phase is not None:
        fid = fid * np.exp(-1j * phase)
    output = np.empty((*t2star.shape, basis.shape[0]))
    for index in np.ndindex(t2star.shape):
        design = (basis.T * np.exp(-times[:, None] / (t2star[index] + 1e-10))
                  * np.exp(2j * np.pi * b0[index] * times[:, None]))
        a = np.concatenate([design.real, design.imag], axis=0)
        y = np.concatenate([fid[index].real, fid[index].imag])
        output[index] = np.linalg.lstsq(a, y, rcond=None)[0]
    return output


def run_benchmark(noise_sd=0.0, seed=20261025):
    """Generate six controls with independent phase and concentration parameters."""
    if not np.isfinite(noise_sd) or noise_sd < 0:
        raise ValueError("noise_sd must be finite and nonnegative")
    rng = np.random.default_rng(seed)
    shape, n_spectral, dt = (6, 6, 1), 256, .0005
    times = np.arange(n_spectral) * dt
    axis = np.linspace(-.025, .025, shape[0])
    positions = np.stack(np.meshgrid(axis, axis, [0.], indexing="ij"), axis=-1)
    vertices = [(-.05, -.05, 0), (.05, -.05, 0), (.05, .05, 0), (-.05, .05, 0)]
    field = np.array([loop_field(vertices, p, 1)[2]
                      for p in positions.reshape(-1, 3)]).reshape(shape)
    windows = [[.008, .048], [.060, .090]]
    sensitivity = [[0, float(times[-1]), 1]]
    current = np.array([.002, -.001])
    enc = CurrentEncoding(times, windows, sensitivity)
    basis_dict = make_lorentzian_basis({"mIns": 3.56, "Cr": 3.03, "NAA": 2.01},
                                      n_spectral, dt, linewidth_hz=3)
    names = list(basis_dict)
    basis = np.array(list(basis_dict.values()))
    baseline = np.broadcast_to([7., 8., 12.], (*shape, 3)).copy()
    roi = np.zeros(shape, dtype=bool)
    roi[2:4, 2:4, :] = True
    t2 = np.full(shape, .08)
    b0 = np.broadcast_to(np.linspace(-4, 4, shape[0])[:, None, None], shape).copy()
    cases = [
        ("null", 0., 0., 0.),
        ("phase_only", 1., 0., 0.),
        ("concentration_only", 0., .35, 0.),
        ("combined", 1., .35, 0.),
        ("reverse_polarity", -1., 0., 0.),
        ("delayed_current", 1., 0., .003),
    ]
    arrays = dict(positions_m=positions, field_per_amp_t=field, sample_times_s=times,
                  basis_fids=basis, roi_mask=roi, t2star_s=t2, b0_hz=b0)
    rows = []
    for name, polarity, delta, delay in cases:
        conc = baseline.copy()
        conc[..., 0] += delta * roi
        phase = enc.phase(field, current * polarity, phase_sign=1, current_delay_s=delay)
        model = TissueModel(np.full(shape, 2), {n: conc[..., i] for i, n in enumerate(names)},
                            np.full(shape, 1.3), t2, b0, (10, 10, 10))
        kspace = simulate_mrsi(model, basis_dict, n_spectral, dt, phase_rad=phase)
        if noise_sd:
            kspace += noise_sd * (rng.normal(size=kspace.shape) + 1j * rng.normal(size=kspace.shape))
        image = np.fft.ifftn(kspace, axes=(0, 1, 2))
        blind = fit_concentrations(image, basis, times, t2, b0)
        aware = fit_concentrations(image, basis, times, t2, b0, phase)
        # An incorrect zero-delay assumption exposes timing-model mismatch.
        assumed_phase = enc.phase(field, current * polarity, phase_sign=1)
        wrong_timing = fit_concentrations(image, basis, times, t2, b0, assumed_phase)
        rows.append(dict(
            case=name, current_scale=polarity, current_delay_s=delay,
            true_roi_mins_change_mM=delta,
            unaware_roi_mins_change_mM=float(np.mean(blind[..., 0][roi]) - 7),
            oracle_roi_mins_change_mM=float(np.mean(aware[..., 0][roi]) - 7),
            zero_delay_roi_mins_change_mM=float(np.mean(wrong_timing[..., 0][roi]) - 7),
            unaware_max_mins_error_mM=float(np.max(abs(blind[..., 0] - conc[..., 0]))),
            oracle_max_mins_error_mM=float(np.max(abs(aware[..., 0] - conc[..., 0]))),
        ))
        arrays.update({f"{name}_kspace": kspace, f"{name}_phase_rad": phase,
                       f"{name}_concentrations_mM": conc, f"{name}_oracle_fit_mM": aware,
                       f"{name}_unaware_fit_mM": blind})
    config = dict(
        shape=shape, n_spectral=n_spectral, dwell_time_s=dt, metabolite_order=names,
        loop_vertices_m=vertices, loop_closure="last vertex to first", b0_unit_vector=[0, 0, 1],
        frame="right-handed synthetic scanner Cartesian, coordinates in metres",
        current_windows_s=windows, base_current_a=current.tolist(), sensitivity_intervals=sensitivity,
        phase_sign=1, phase_convention="+gamma integral_0^t q(u) B_parallel(u) du",
        reference_mapping="Negate neuroelectric-mri Mx+iMy analytical phase for this MRS convention",
        noise_sd_per_real_imag_kspace_component=noise_sd, seed=seed,
        basis="Illustrative single Lorentzian per metabolite; NOT a coupled-spin mI basis",
        basis_centre_freq_hz=123.25e6, basis_shifts_ppm={"mIns": 3.56, "Cr": 3.03, "NAA": 2.01},
        basis_linewidth_hz=3., concentration_units="nominal mM under unit-amplitude synthetic basis",
    )
    return rows, arrays, config


def write_benchmark(output, noise_sd=0.0, seed=20261025):
    """Write complex samples, ground truth, metrics and source-hash provenance."""
    rows, arrays, config = run_benchmark(noise_sd, seed)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "signals.npz", **arrays)
    with (output / "metrics.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "parameters.json").write_text(json.dumps(config, indent=2) + "\n")
    sources = {"generator": Path(__file__), "loop_field": Path(inspect.getfile(loop_field)),
               "current_encoding": Path(inspect.getfile(CurrentEncoding)),
               "mrsi_sim": Path(inspect.getfile(simulate_mrsi))}
    mrs_root = sources["current_encoding"].parents[2]
    try:
        revision = subprocess.check_output(["git", "-C", str(mrs_root), "rev-parse", "HEAD"],
                                           text=True, stderr=subprocess.DEVNULL).strip()
        dirty = bool(subprocess.check_output(["git", "-C", str(mrs_root), "status", "--porcelain"],
                                             text=True))
    except (FileNotFoundError, subprocess.CalledProcessError):
        revision, dirty = None, None
        # pip VCS installations retain the immutable source revision even
        # though their installed modules are no longer in a Git checkout.
        try:
            direct = importlib.metadata.distribution("mrs-jax").read_text("direct_url.json")
            if direct:
                revision = json.loads(direct).get("vcs_info", {}).get("commit_id")
        except importlib.metadata.PackageNotFoundError:
            pass
    manifest = dict(
        schema="neuroelectric-mri.synthetic-mrs-current/0.1", data_kind="synthetic",
        scanner_data_available=False, measured_current_available=False,
        mrs_jax_revision=revision, mrs_jax_worktree_dirty=dirty,
        python_version=platform.python_version(),
        dependencies={n: importlib.metadata.version(n) for n in ("numpy", "scipy")},
        arrays={k: {"shape": v.shape, "dtype": str(v.dtype)} for k, v in arrays.items()},
        source_sha256={k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in sources.items()},
        files_sha256={n: hashlib.sha256((output / n).read_bytes()).hexdigest()
                      for n in ("signals.npz", "metrics.csv", "parameters.json")},
        limitations=["Oracle phase uses known simulation truth; no conductivity inference",
                     "Matched forward/fitting basis is a software check, not independent mI validation",
                     "No coils, motion, lipid baseline, hardware response or scanner sequence",
                     "Phase-encoded MRSI assumes the same excitation-relative waveform for every encode",
                     "Synthetic waveform amplitudes are not hardware operating recommendations"],
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for row in rows:
        print(f"{row['case']:20s} true={row['true_roi_mins_change_mM']:+.4f} "
              f"unaware={row['unaware_roi_mins_change_mM']:+.4f} "
              f"oracle={row['oracle_roi_mins_change_mM']:+.4f} nominal mM")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--noise-sd", type=float, default=0.)
    parser.add_argument("--seed", type=int, default=20261025)
    args = parser.parse_args()
    write_benchmark(args.output, args.noise_sd, args.seed)
