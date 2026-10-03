import tempfile
import unittest
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.octopus_results import (
    align_isolated_population,
    discover_density_frames,
    extract_paired_plateau,
    integrate_density_on_rectilinear_grid,
    load_octopus_density,
    parse_td_coordinates,
)


class OctopusResultTests(unittest.TestCase):
    def test_parse_coordinates_reads_position_block_only(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "coordinates"
            path.write_text(
                "# iter time positions velocities forces\n"
                "0 0.0 1.0 2.0 3.0 9.0 9.0 9.0\n"
                "100 1.0D+00 0.5 2.0 3.0 8.0 8.0 8.0\n",
                encoding="utf-8",
            )
            series = parse_td_coordinates(path, natoms=1, projectile_index=0)
            np.testing.assert_array_equal(series.iterations, [0, 100])
            np.testing.assert_allclose(series.times, [0.0, 1.0])
            np.testing.assert_allclose(series.projectile_positions[:, 0], [1.0, 0.5])

    def test_discover_density_frames_is_strict(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "td.0000000"
            second = root / "td.0000100"
            first.mkdir()
            second.mkdir()
            (first / "density.ncdf").touch()
            (second / "density.ncdf").touch()
            frames = discover_density_frames(root)
            self.assertEqual(sorted(frames), [0, 100])

    def test_density_integral_uses_native_grid_volume(self):
        axis = np.linspace(-1.0, 1.0, 21)
        density = np.ones((axis.size, axis.size, axis.size))
        population = integrate_density_on_rectilinear_grid(
            density,
            (axis, axis, axis),
            np.zeros(3),
            0.5,
            transition_half_width=0.1,
        )
        self.assertAlmostEqual(population, 4.0 * np.pi * 0.5**3 / 3.0, delta=0.08)

    def test_load_octopus_density_netcdf(self):
        from netCDF4 import Dataset

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "density.ncdf"
            axis = np.linspace(-1.0, 1.0, 5)
            with Dataset(path, "w") as dataset:
                for name in ("x", "y", "z"):
                    dataset.createDimension(name, axis.size)
                    coordinate = dataset.createVariable(name, "f8", (name,))
                    coordinate[:] = axis
                density = dataset.createVariable("density", "f8", ("x", "y", "z"))
                density[...] = 1.0
            values, axes = load_octopus_density(path)
            self.assertEqual(values.shape, (5, 5, 5))
            for loaded in axes:
                np.testing.assert_allclose(loaded, axis)

    def test_reference_alignment_selects_recorded_samples_only(self):
        source_x = np.asarray([-2.0, -1.0, 0.0])
        source_population = np.asarray([8.0, 7.8, 7.5])
        aligned = align_isolated_population(
            np.asarray([0.0, -2.0, -1.0]), source_x, source_population
        )
        np.testing.assert_array_equal(aligned, [7.5, 8.0, 7.8])
        for position in (-1.5, -0.5, 0.1):
            with self.subTest(position=position), self.assertRaisesRegex(ValueError, "recorded"):
                align_isolated_population(np.asarray([position]), source_x, source_population)

    def test_reference_alignment_rejects_missing_and_ambiguous_samples(self):
        with self.assertRaisesRegex(ValueError, "finite"):
            align_isolated_population(np.asarray([0.0]), np.asarray([0.0]), np.asarray([np.nan]))
        with self.assertRaisesRegex(ValueError, "exactly one"):
            align_isolated_population(np.asarray([0.0]), np.asarray([-1e-11, 1e-11]), np.asarray([7.0, 8.0]))
        with self.assertRaisesRegex(ValueError, "exactly one"):
            align_isolated_population(np.asarray([0.0]), np.asarray([0.0, 1e-11]), np.asarray([7.0, 8.0]))

    def test_paired_plateau_comes_from_density_frames(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stages = {name: root / name for name in ("production", "isolated")}
            for stage in stages.values():
                (stage / "td.general").mkdir(parents=True)
                coordinate_lines = [
                    f"{iteration} {float(iteration)} {0.4 - 0.1 * iteration} 0.0 0.0"
                    for iteration in range(5)
                ]
                (stage / "td.general" / "coordinates").write_text(
                    "\n".join(coordinate_lines) + "\n", encoding="utf-8"
                )
                for iteration in range(5):
                    frame = stage / f"td.{iteration:07d}"
                    frame.mkdir()
                    (frame / "density.ncdf").touch()

            axis = np.asarray([-1.0, 0.0, 1.0])

            def loader(path):
                iteration = int(Path(path).parent.name.split(".")[1])
                is_production = Path(path).parents[1].name == "production"
                value = 1.0 - 0.01 * iteration if is_production else 1.0
                return np.full((3, 3, 3), value), (axis, axis, axis)

            mean_loss, provenance = extract_paired_plateau(
                stages["production"],
                stages["isolated"],
                radius=2.0,
                transition_half_width=0.01,
                tail_fraction=0.4,
                production_natoms=1,
                isolated_natoms=1,
                density_loader=loader,
            )
            self.assertAlmostEqual(mean_loss, 0.81)
            self.assertEqual(provenance["frame_count"], 5)
            self.assertEqual(len(provenance["source_files"]), 12)


if __name__ == "__main__":
    unittest.main()
