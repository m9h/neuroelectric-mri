"""Independent field limits, symmetry and sequence-timing regressions."""

import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
from current_phantom import loop_field, segment_field, waveform_phase, write_reference
from phase_budget import GAMMA_RAD_S_T, MU0_APPROX, loop_axis_field_tesla


class CurrentPhantomTests(unittest.TestCase):
    def setUp(self):
        self.a = 0.05
        self.loop = [(-self.a, -self.a, 0), (self.a, -self.a, 0),
                     (self.a, self.a, 0), (-self.a, self.a, 0)]

    def test_square_center_closed_form(self):
        b = loop_field(self.loop, (0, 0, 0), 1)
        expected = MU0_APPROX * math.sqrt(2) / (math.pi * self.a)
        self.assertAlmostEqual(b[2] / expected, 1, places=13)
        self.assertEqual(b[:2], (0, 0))

    def test_long_wire_limit_and_direction(self):
        b = segment_field((-100, 0, 0), (100, 0, 0), (0, 0.02, 0), 1)
        expected = MU0_APPROX / (2 * math.pi * 0.02)
        self.assertLess(abs(b[2] / expected - 1), 3e-8)
        self.assertEqual(b[:2], (0, 0))

    def test_polygon_converges_to_circular_loop(self):
        errors = []
        exact = loop_axis_field_tesla(1, self.a, 0.02)
        for count in (32, 128):
            vertices = [(self.a * math.cos(2 * math.pi * k / count),
                         self.a * math.sin(2 * math.pi * k / count), 0)
                        for k in range(count)]
            errors.append(abs(loop_field(vertices, (0, 0, 0.02), 1)[2] / exact - 1))
        self.assertLess(errors[1], errors[0] / 10)
        self.assertLess(errors[1], 2e-4)

    def test_polarity_and_path_reversal(self):
        p = (0.01, -0.02, 0.03)
        b = loop_field(self.loop, p, 0.001)
        for opposite in (loop_field(self.loop, p, -0.001),
                         loop_field(list(reversed(self.loop)), p, 0.001)):
            for x, y in zip(b, opposite):
                self.assertAlmostEqual(x, -y, delta=1e-20)

    def test_rotation_creates_longitudinal_null(self):
        rotated = [(z, y, -x) for x, y, z in self.loop]
        b = loop_field(rotated, (0, 0, 0), 1)
        self.assertGreater(abs(b[0]), 0)
        self.assertEqual(b[2], 0)

    def test_refocusing_static_and_matched_bipolar_current(self):
        q = [(0, 0.02, 1), (0.02, 0.04, -1)]
        self.assertEqual(waveform_phase(1e-6, [(0, 0.04, 0.001)], q), 0)
        phase = waveform_phase(1e-6, [(0, 0.02, 0.001), (0.02, 0.04, -0.001)], q)
        self.assertAlmostEqual(phase / (-GAMMA_RAD_S_T * 1e-9 * 0.04), 1)

    def test_delay_bias_has_boundary_and_refocusing_contributions(self):
        q = [(0, 0.02, 1), (0.02, 0.04, -1)]
        i = [(0, 0.02, 1), (0.02, 0.04, -1)]
        baseline = waveform_phase(1e-9, i, q)
        for delay in (-1e-4, 1e-4):
            shifted = waveform_phase(1e-9, i, q, delay)
            # Loss |delay| at acquisition boundary plus 2|delay| at RF sign change.
            expected = -GAMMA_RAD_S_T * 1e-9 * (0.04 - 3 * abs(delay))
            self.assertAlmostEqual(shifted / expected, 1, places=13)
            self.assertGreater(shifted, baseline)

    def test_gap_and_zero_sensitivity(self):
        self.assertEqual(waveform_phase(1, [(0, 1, 1)], [(2, 3, 1)]), 0)
        self.assertEqual(waveform_phase(1, [(0, 1, 1)], [(0, 1, 0)]), 0)

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):
            segment_field((0, 0, 0), (1, 0, 0), (0.5, 0, 0), 1)
        with self.assertRaises(ValueError):
            segment_field((0, 0, 0), (0, 0, 0), (1, 1, 1), 1)
        with self.assertRaises(ValueError):
            waveform_phase(1, [(0, 2, 1), (1, 3, 1)], [(0, 3, 1)])
        with self.assertRaises(ValueError):
            waveform_phase(1, [(0, 1, 1)], [(0, 1, 2)])
        with self.assertRaises(ValueError):
            waveform_phase(float("nan"), [], [])

    def test_fixture_is_deterministic_and_explicitly_synthetic(self):
        with tempfile.TemporaryDirectory() as folder:
            a, b = Path(folder) / "a", Path(folder) / "b"
            write_reference(a)
            write_reference(b)
            for name in ("encoding.csv", "field_phase.csv", "manifest.json"):
                self.assertEqual((a / name).read_bytes(), (b / name).read_bytes())
            manifest = json.loads((a / "manifest.json").read_text())
            self.assertEqual(manifest["data_kind"], "synthetic_analytical")
            self.assertFalse(manifest["scanner_data_available"])


if __name__ == "__main__":
    unittest.main()
