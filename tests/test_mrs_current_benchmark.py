"""Optional cross-repository checks; install requirements-simulation.txt first."""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "examples"))
AVAILABLE = importlib.util.find_spec("mrs_jax") is not None
if AVAILABLE:
    import numpy as np
    from mrs_jax.current_encoding import CurrentEncoding
    from current_phantom import waveform_phase
    from mrs_current_benchmark import run_benchmark, write_benchmark


@unittest.skipUnless(AVAILABLE, "Install the pinned mrs-jax simulation requirements")
class SpectroscopyBenchmarkTests(unittest.TestCase):
    def test_independent_scalar_reference_sign_and_delay(self):
        times = np.array([0, .007, .02, .03, .04])
        current = [(0, .02, .001), (.02, .04, -.001)]
        sensitivity = [(0, .02, 1), (.02, .04, -1)]
        enc = CurrentEncoding(times, [v[:2] for v in current], sensitivity)
        for delay in (-.001, 0., .001):
            actual = enc.phase(2e-6, [v[2] for v in current], phase_sign=1,
                               current_delay_s=delay)
            expected = []
            for t in times:
                clipped = [(a, min(b, t), q) for a, b, q in sensitivity if a < t]
                expected.append(-waveform_phase(2e-6, current, clipped, delay))
            np.testing.assert_allclose(actual, expected, atol=1e-15)

    def test_null_concentration_phase_and_combined_recovery(self):
        rows, arrays, _ = run_benchmark()
        result = {r["case"]: r for r in rows}
        for row in rows:
            self.assertLess(row["oracle_max_mins_error_mM"], 2e-5)
            self.assertAlmostEqual(row["oracle_roi_mins_change_mM"],
                                   row["true_roi_mins_change_mM"], places=5)
        self.assertLess(abs(result["null"]["unaware_roi_mins_change_mM"]), 2e-5)
        self.assertGreater(result["phase_only"]["unaware_max_mins_error_mM"], 1e-3)
        self.assertGreater(abs(result["delayed_current"]["zero_delay_roi_mins_change_mM"]), 1e-4)
        np.testing.assert_array_equal(arrays["phase_only_concentrations_mM"],
                                      arrays["null_concentrations_mM"])
        np.testing.assert_allclose(arrays["phase_only_phase_rad"],
                                    -arrays["reverse_polarity_phase_rad"])

    def test_seeded_noise_and_release_hashes(self):
        a = run_benchmark(noise_sd=.01, seed=19)[1]
        b = run_benchmark(noise_sd=.01, seed=19)[1]
        np.testing.assert_array_equal(a["combined_kspace"], b["combined_kspace"])
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            write_benchmark(path)
            manifest = json.loads((path / "manifest.json").read_text())
            self.assertEqual(manifest["data_kind"], "synthetic")
            self.assertFalse(manifest["scanner_data_available"])
            for name, expected in manifest["files_sha256"].items():
                self.assertEqual(hashlib.sha256((path / name).read_bytes()).hexdigest(), expected)


if __name__ == "__main__":
    unittest.main()
