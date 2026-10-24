"""Analytical oracles and physical invariants for the illustrative calculations."""

import math
import unittest

from examples.phase_budget import (
    encoded_phase, ideal_polarity_pairs, loop_axis_field_tesla,
    normalized_signal, phase_radians, symmetric_magnitude_loss,
)


class PhaseBudgetTests(unittest.TestCase):
    def test_reference_phase_in_si_units(self):
        # 0.1 nT for 40 ms -> 0.06131157 degrees in absolute value.
        self.assertAlmostEqual(phase_radians(0.1e-9, 0.04), -0.00107008875987, places=12)
        self.assertAlmostEqual(abs(math.degrees(phase_radians(0.1e-9, 0.04))),
                               0.06131156964, places=10)

    def test_uniform_field_changes_phase_without_dephasing(self):
        phase = phase_radians(1e-9, 0.04)
        signal = normalized_signal([phase] * 100)
        self.assertAlmostEqual(abs(signal), 1.0, places=14)
        self.assertAlmostEqual(math.atan2(signal.imag, signal.real), phase, places=14)

    def test_equal_opposite_fields_cancel_phase_not_magnitude_loss(self):
        phase = phase_radians(0.1e-9, 0.04)
        signal = normalized_signal([phase, -phase])
        self.assertEqual(signal.imag, 0.0)
        self.assertAlmostEqual(symmetric_magnitude_loss(phase), 1 - abs(signal), places=15)
        # Small-phase asymptotic prediction is independent of the exact formula.
        self.assertAlmostEqual(symmetric_magnitude_loss(phase) / (phase**2 / 2), 1, places=6)

    def test_static_field_refocuses_but_bipolar_field_accumulates(self):
        field = 2e-9
        self.assertEqual(encoded_phase([(0.02, field, 1), (0.02, field, -1)]), 0)
        self.assertAlmostEqual(encoded_phase([(0.02, field, 1), (0.02, -field, -1)]),
                               phase_radians(field, 0.04))
        self.assertEqual(encoded_phase([(0.02, field, 1), (0.02, -field, 1)]), 0)

    def test_polarity_reversal_cancels_common_background(self):
        neural = phase_radians(0.2e-9, 0.04)
        background = 0.5
        difference = (background + neural) - (background - neural)
        self.assertAlmostEqual(difference, 2 * neural, places=14)

    def test_loop_center_and_axial_symmetry(self):
        # Circular loop: at center B = mu0 I/(2R) = 6.2831853e-11 T.
        self.assertAlmostEqual(loop_axis_field_tesla(1e-6, 0.01) / 1e-11,
                               6.283185307179586, places=12)
        self.assertEqual(loop_axis_field_tesla(1e-6, 0.01, 0.03),
                         loop_axis_field_tesla(1e-6, 0.01, -0.03))
        self.assertEqual(loop_axis_field_tesla(-1e-6, 0.01),
                         -loop_axis_field_tesla(1e-6, 0.01))

    def test_pair_count_meets_ideal_noise_budget(self):
        phase, sigma, snr = 0.001, 0.01, 5
        count = ideal_polarity_pairs(phase, sigma, snr)
        achieved = math.sqrt(2 * count) * phase / sigma
        previous = math.sqrt(2 * (count - 1)) * phase / sigma
        self.assertGreaterEqual(achieved, snr)
        self.assertLess(previous, snr)
        self.assertIsNone(ideal_polarity_pairs(0))

    def test_small_phase_loss_does_not_round_to_zero(self):
        self.assertGreater(symmetric_magnitude_loss(1e-12), 0)

    def test_invalid_physical_inputs(self):
        cases = [lambda: phase_radians(1, -1), lambda: phase_radians(math.nan, 1),
                 lambda: normalized_signal([]), lambda: loop_axis_field_tesla(1, 0),
                 lambda: symmetric_magnitude_loss(math.pi),
                 lambda: ideal_polarity_pairs(1, -1),
                 lambda: encoded_phase([(1, 1, 0)])]
        for call in cases:
            with self.subTest(call=call), self.assertRaises(ValueError):
                call()


if __name__ == "__main__":
    unittest.main()
