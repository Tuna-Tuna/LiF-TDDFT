import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.legacy_run_process import replay


class LegacyReplayTests(unittest.TestCase):
    def test_two_state_replay_and_domain_gate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            arrays = {
                "vcal": np.asarray([[0.2, 0.4]]),
                "zcal": np.asarray([[2.0, 3.0]]),
                "result_hsemt": np.ones((2, 2)),
                "result_esemt": np.zeros((2, 2)),
                "result_loss": np.asarray([[0.4, 0.5], [0.1, 0.2]]),
                "Zhr": np.asarray([[2.5, 0.0], [1.5, 0.0]]),
                "zPcapture": np.asarray([[0.5, 0.0], [0.6, 0.0]]),
                "zPloss": np.asarray([[0.2, 0.0], [0.3, 0.0]]),
                "Pfin": np.asarray([[0.0, 0.5], [0.0, 0.6]]),
                "last_non_zero": np.asarray([[0.5], [0.6]]),
                "Zturn": np.asarray([[2.5, 1.5]]),
                "expe": np.asarray([[0.2, 0.5], [0.4, 0.6]]),
            }
            records = []
            for name, value in arrays.items():
                file_name = f"{name}.csv"
                np.savetxt(root / file_name, value, delimiter=",")
                records.append({"name": name, "file": file_name})
            (root / "manifest.json").write_text(
                json.dumps({"source_mat": "fixture", "arrays": records}), encoding="utf-8"
            )
            result = replay(root)
            self.assertAlmostEqual(result["velocity_summaries"][0]["replay_absolute_error"], 0.0)
            self.assertEqual(result["velocity_summaries"][1]["out_of_domain_events"], 1)


if __name__ == "__main__":
    unittest.main()
