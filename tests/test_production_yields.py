import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "compute_sequential_yields", ROOT / "scripts" / "compute_sequential_yields.py"
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

FIGURE_SPEC = importlib.util.spec_from_file_location(
    "export_detachment_probabilities", ROOT / "scripts" / "export_detachment_probabilities.py"
)
assert FIGURE_SPEC is not None and FIGURE_SPEC.loader is not None
FIGURE_MODULE = importlib.util.module_from_spec(FIGURE_SPEC)
FIGURE_SPEC.loader.exec_module(FIGURE_MODULE)


class ProductionYieldSourceTests(unittest.TestCase):
    def _write(self, directory: str, payload: dict) -> Path:
        path = Path(directory) / "grid.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_legacy_grid_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self._write(
                directory,
                {
                    "source_kind": "legacy_matlab_workspace",
                    "data_policy": {},
                },
            )
            with self.assertRaisesRegex(ValueError, "directly from Octopus"):
                MODULE.load_mean_loss_grid(path)

    def test_direct_tddft_grid_is_accepted_without_scaling(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = {
                "source_kind": "octopus_rt_tddft_density",
                "data_policy": {
                    "normalization": "none",
                    "scale_factor": 1.0,
                    "value_transform": "none",
                    "velocity_specific_adjustments": "none",
                },
                "heights_bohr": [2.0, 3.0],
                "velocities_au": [0.1, 0.2],
                "mean_loss": [[0.2, 0.3], [0.1, 0.2]],
                "runs": [
                    {
                        "production_manifest_sha256": "0" * 64,
                        "isolated_manifest_sha256": "1" * 64,
                        "provenance": {
                            "source_files": [{"path": "density.ncdf", "sha256": "2" * 64}]
                        },
                    }
                ],
            }
            model, loaded = MODULE.load_mean_loss_grid(self._write(directory, payload))
            self.assertEqual(loaded["source_kind"], "octopus_rt_tddft_density")
            self.assertAlmostEqual(model.probability({"surface_height": 2.5, "v_parallel": 0.1}), 0.15)

    def test_detachment_export_preserves_single_electron_probability(self):
        self.assertEqual(FIGURE_MODULE.pdet_from_mean_loss(0.25), 0.25)
        self.assertEqual(FIGURE_MODULE.pdet_from_mean_loss(1.0), 1.0)
        with self.assertRaises(ValueError):
            FIGURE_MODULE.pdet_from_mean_loss(1.25)

    def test_sequential_cli_uses_equation10_and_explicit_capture_inputs(self):
        # Complete I/O integration with arithmetic fixtures, not research data.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            payload = {
                "source_kind": "octopus_rt_tddft_density",
                "data_policy": {"normalization": "none", "scale_factor": 1.0,
                                "value_transform": "none", "velocity_specific_adjustments": "none"},
                "heights_bohr": [1.2, 3.5], "velocities_au": [0.2],
                "mean_loss": [[0.7], [0.1]],
                "runs": [{"production_manifest_sha256": "0"*64,
                          "isolated_manifest_sha256": "1"*64,
                          "provenance": {"source_files": [{"path": "fixture", "sha256": "2"*64}]}}],
            }
            grid = self._write(directory, payload)
            config = root/"capture.yaml"
            config.write_text("parameters: {gamma_bohr_inverse: 0.7}\ngamma_source: arithmetic fixture\n", encoding="utf-8")
            events = root/"events.csv"
            events.write_text("trajectory_id,event_index,event_time,surface_height,v_parallel,energy_defect_au\n"
                              "test,0,0,1.2,0.2,0.01\ntest,1,1,3.5,0.2,0.01\n", encoding="utf-8")
            output = root/"yields.json"
            MODULE.main([str(events), "--config", str(config), "--tddft-grid", str(grid), "--output", str(output)])
            result = json.loads(output.read_text(encoding="utf-8"))
            capture = 0.5/math.cosh(math.pi*(.01+.2**2/2)/(2*.7*.2))**2
            expected = (1-capture)*capture+capture*(1-.1)
            summary = result["summaries"][0]
            self.assertEqual(set(summary["final_states"]), {"F_minus", "F_neutral"})
            self.assertAlmostEqual(summary["final_f_minus_fraction"], expected, places=15)


if __name__ == "__main__":
    unittest.main()
