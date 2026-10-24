"""Ideal field-to-phase budgets; no anatomical or scanner realism is implied.

SI units inside the functions. Phase convention: dM/dt = M x (gamma B),
signal = Mx + i My, so a positive Bz gives negative signal phase.
"""

from __future__ import annotations

import argparse
import cmath
import csv
import math
from pathlib import Path

GAMMA_RAD_S_T = 2 * math.pi * 42.57747892e6
MU0_APPROX = 4 * math.pi * 1e-7


def phase_radians(field_t: float, duration_s: float) -> float:
    """Uniform static longitudinal perturbation, with no refocusing."""
    if not math.isfinite(field_t) or not math.isfinite(duration_s) or duration_s < 0:
        raise ValueError("Field must be finite and duration finite/nonnegative")
    return -GAMMA_RAD_S_T * field_t * duration_s


def encoded_phase(intervals) -> float:
    """Exact piecewise-constant integral of -gamma * sensitivity * Bz.

    Each interval is (duration_s, field_t, sensitivity); sensitivity is +1/-1
    for ideal refocusing. These are toggling-frame weights, not RF simulations.
    """
    contributions = []
    for duration_s, field_t, sensitivity in intervals:
        if sensitivity not in (-1, 1):
            raise ValueError("Sensitivity must be +1 or -1")
        contributions.append(sensitivity * phase_radians(field_t, duration_s))
    return math.fsum(contributions)


def normalized_signal(phases) -> complex:
    """Equal-weight, nonrelaxing isochromats, all initially transverse."""
    phases = tuple(phases)
    if not phases or any(not math.isfinite(p) for p in phases):
        raise ValueError("At least one finite phase is required")
    return sum(cmath.exp(1j * p) for p in phases) / len(phases)


def symmetric_magnitude_loss(phase: float) -> float:
    """Fractional magnitude loss for equal +phase and -phase populations.

    Stable at small phase; restricted to the positive-cosine branch.
    """
    if not math.isfinite(phase) or abs(phase) > math.pi / 2:
        raise ValueError("This small-phase example requires abs(phase) <= pi/2")
    return 2 * math.sin(phase / 2) ** 2


def loop_axis_field_tesla(current_a: float, radius_m: float, z_m: float = 0.0) -> float:
    """Biot-Savart field on the axis of one ideal thin circular wire loop.

    Vacuum permeability is approximated as 4*pi*1e-7 H/m. This is an external
    closed-loop phantom, not a neuronal source or conductive-tissue solution.
    """
    if not all(map(math.isfinite, (current_a, radius_m, z_m))) or radius_m <= 0:
        raise ValueError("Inputs must be finite and radius strictly positive")
    return MU0_APPROX * current_a * radius_m**2 / (2 * (radius_m**2 + z_m**2) ** 1.5)


def ideal_polarity_pairs(phase: float, sigma_phase_rad: float = 0.01,
                         target_snr: float = 5.0) -> int | None:
    """Pairs needed for ideal SNR with independent equal-variance phase noise.

    Difference of +I/-I observations has signal 2*phase and noise sqrt(2)*sigma.
    Returns None for zero signal. Does not model power, multiple comparisons,
    correlated noise, nuisance bias, phase wrapping or available scan duration.
    """
    if not all(map(math.isfinite, (phase, sigma_phase_rad, target_snr))):
        raise ValueError("Inputs must be finite")
    if sigma_phase_rad <= 0 or target_snr <= 0:
        raise ValueError("Noise and target SNR must be positive")
    if phase == 0:
        return None
    return max(1, math.ceil((target_snr * sigma_phase_rad / (math.sqrt(2) * abs(phase))) ** 2))


def budget_rows():
    for duration_ms in (10, 40):
        for field_nt in (0.01, 0.07, 0.1, 0.2, 0.6, 1.0):
            phase = phase_radians(field_nt * 1e-9, duration_ms * 1e-3)
            yield {
                "field_nT": field_nt,
                "duration_ms": duration_ms,
                "phase_rad": phase,
                "phase_deg": math.degrees(phase),
                "symmetric_fractional_magnitude_loss": symmetric_magnitude_loss(phase),
                "ideal_polarity_pairs_sigma_0p01rad_snr5": ideal_polarity_pairs(phase),
            }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, help="Write deterministic derived examples")
    args = parser.parse_args()
    rows = list(budget_rows())
    if args.csv:
        with args.csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} analytical example rows to {args.csv}")
        return
    print("Idealized uniform-field phase budgets; no detection claim")
    print(" B (nT)  time (ms)    phase (rad)   |phase| (deg)   +/-field magnitude loss")
    for row in rows:
        print(f"{row['field_nT']:7.2f} {row['duration_ms']:10d} "
              f"{row['phase_rad']:14.7g} {abs(row['phase_deg']):15.7g} "
              f"{row['symmetric_fractional_magnitude_loss']:25.7g}")
    field = 0.1e-9
    steady = encoded_phase([(0.02, field, 1), (0.02, field, -1)])
    bipolar = encoded_phase([(0.02, field, 1), (0.02, -field, -1)])
    print(f"\nStatic field refocused at TE=40 ms: {steady:.6g} rad")
    print(f"Field reversing at refocusing pulse: {bipolar:.6g} rad")
    print(f"1 microamp, 1 cm loop, center field: {loop_axis_field_tesla(1e-6, 0.01)*1e9:.6g} nT")
    print("Pair counts in CSV assume independent 0.01-rad phase noise and SNR=5.")


if __name__ == "__main__":
    main()
