"""Analytical closed-wire-loop and timing references in SI units.

No tissue solver, scanner driver, hardware interface or measured data.
Ideal thin wires in free space; piecewise-constant currents and ideal RF toggling.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from phase_budget import GAMMA_RAD_S_T, MU0_APPROX


def vector(value):
    value = tuple(value)
    if len(value) != 3 or not all(math.isfinite(x) for x in value):
        raise ValueError("Expected three finite coordinates")
    return value


def dot(a, b):
    return math.fsum(x * y for x, y in zip(a, b))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def segment_field(start_m, end_m, observation_m, current_a):
    """Exact finite straight-wire field, current flowing start -> end.

    A segment is a mathematical component; loop_field supplies a closed return
    path. Observation on the ideal wire is singular and rejected.
    """
    start, end, point = map(vector, (start_m, end_m, observation_m))
    if not math.isfinite(current_a):
        raise ValueError("Current must be finite")
    length_vector = subtract(end, start)
    length = math.sqrt(dot(length_vector, length_vector))
    if length == 0:
        raise ValueError("Wire segment must have nonzero length")
    direction = tuple(x / length for x in length_vector)
    relative = subtract(point, start)
    along_start = dot(relative, direction)
    along_end = along_start - length
    perpendicular = tuple(x - along_start * u for x, u in zip(relative, direction))
    rho2 = dot(perpendicular, perpendicular)
    if rho2 == 0:
        if 0 <= along_start <= length:
            raise ValueError("Observation lies on an ideal thin wire")
        return (0.0, 0.0, 0.0)
    angular_factor = (along_start / math.sqrt(rho2 + along_start**2)
                      - along_end / math.sqrt(rho2 + along_end**2))
    factor = MU0_APPROX * current_a * angular_factor / (4 * math.pi * rho2)
    return tuple(factor * x for x in cross(direction, perpendicular))


def loop_field(vertices_m, observation_m, current_a):
    """Sum fields of a closed polygon; final vertex connects to the first.

    Supply distinct consecutive vertices without repeating the first at end.
    """
    vertices = tuple(map(vector, vertices_m))
    if len(vertices) < 3:
        raise ValueError("At least three vertices required")
    fields = [segment_field(a, b, observation_m, current_a)
              for a, b in zip(vertices, vertices[1:] + vertices[:1])]
    return tuple(math.fsum(b[i] for b in fields) for i in range(3))


def intervals(values):
    """Validate sorted, nonoverlapping (start_s, end_s, value) intervals.

    Negative times are allowed for waveform shifts relative to sequence zero.
    Gaps have zero value. Endpoints have zero measure in the phase integral.
    """
    result = tuple(tuple(v) for v in values)
    previous_end = -math.inf
    for interval in result:
        if len(interval) != 3 or not all(math.isfinite(x) for x in interval):
            raise ValueError("Intervals require three finite values")
        start, end, _ = interval
        if start >= end or start < previous_end:
            raise ValueError("Intervals must be ordered, positive-length and nonoverlapping")
        previous_end = end
    return result


def waveform_phase(field_per_amp_t, current_intervals, sensitivity_intervals,
                   current_delay_s=0.0):
    """Exact -gamma * integral B_parallel(t)*q(t) dt.

    Current intervals hold amperes. Sensitivity q is -1, 0 or +1. Positive
    current_delay shifts the current later relative to the sequence, without
    shifting its RF/ADC windows. Field must already be projected along B0.
    """
    if not all(map(math.isfinite, (field_per_amp_t, current_delay_s))):
        raise ValueError("Field coefficient and delay must be finite")
    currents, sensitivity = intervals(current_intervals), intervals(sensitivity_intervals)
    if any(q not in (-1, 0, 1) for _, _, q in sensitivity):
        raise ValueError("Sensitivity must be -1, 0 or +1")
    terms = []
    for start, end, current in currents:
        start, end = start + current_delay_s, end + current_delay_s
        for q_start, q_end, q in sensitivity:
            overlap = max(0.0, min(end, q_end) - max(start, q_start))
            terms.append(overlap * current * q)
    return -GAMMA_RAD_S_T * field_per_amp_t * math.fsum(terms)


def write_reference(output):
    """Write deterministic synthetic field/phase data and a hash manifest."""
    output.mkdir(parents=True, exist_ok=True)
    half_width, current = 0.05, 1e-4
    vertices = [(-half_width, -half_width, 0), (half_width, -half_width, 0),
                (half_width, half_width, 0), (-half_width, half_width, 0)]
    waveform = [(0, 0.02, current), (0.02, 0.04, -current)]
    sensitivity = [(0, 0.02, 1), (0.02, 0.04, -1)]
    delay = 0.0001
    rows = []
    for z in (-0.02, 0.0, 0.02):
        for y in (-0.03, -0.015, 0.0, 0.015, 0.03):
            for x in (-0.03, -0.015, 0.0, 0.015, 0.03):
                field = loop_field(vertices, (x, y, z), current)
                coefficient = field[2] / current
                phase = waveform_phase(coefficient, waveform, sensitivity)
                delayed = waveform_phase(coefficient, waveform, sensitivity, delay)
                rows.append(dict(x_m=x, y_m=y, z_m=z, bx_T=field[0], by_T=field[1],
                                 bz_T=field[2], phase_rad=phase,
                                 delayed_phase_rad=delayed, delay_bias_rad=delayed-phase))
    with (output / "field_phase.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with (output / "encoding.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["start_s", "end_s", "synthetic_current_A", "ideal_sensitivity"])
        writer.writerows((*w, q[2]) for w, q in zip(waveform, sensitivity))
    files = ["encoding.csv", "field_phase.csv"]
    manifest = dict(
        schema="neuroelectric-mri.synthetic-loop/0.1", data_kind="synthetic_analytical",
        description="Closed square loop, free-space thin-wire field; ideal refocused phase",
        geometry=dict(vertices_m=vertices, closure="last vertex to first",
                      frame="right-handed synthetic scanner Cartesian", b0_unit_vector=[0, 0, 1]),
        phase_convention="-gamma integral q(t) B_parallel(t) dt; signal Mx+iMy",
        time_origin="ideal transverse excitation at t=0; acquisition endpoint t=0.04 s",
        current_delay_s=delay, gamma_rad_s_T=GAMMA_RAD_S_T, mu0_H_m=MU0_APPROX,
        measured_current_available=False, scanner_data_available=False,
        limitations=["No conductive tissue or conductivity inversion", "No RF/gradient/ADC simulation",
                     "No noise, voxel averaging, relaxation, hardware or scanner validation"],
        files_sha256={f: hashlib.sha256((output / f).read_bytes()).hexdigest() for f in files},
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {len(rows)} synthetic points and metadata to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    write_reference(parser.parse_args().output)
