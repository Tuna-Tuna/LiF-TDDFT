import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from lif_tddft.analysis.interpolation import DetachmentProbabilityInterpolator, OutOfDomainError
from lif_tddft.analysis.data_contract import (
    assert_no_scaling_policy,
    assert_values_identical,
    validate_probability_values,
)
from lif_tddft.analysis.moving_flux import closure_residual, outward_flux
from lif_tddft.analysis.orthogonalization import cross_overlap, occupied_overlap_descriptors
from lif_tddft.analysis.pauli_exchange_energy import allowed_label
from lif_tddft.analysis.population import (
    moving_sphere_population,
    paired_projectile_mean_loss,
    surface_corrected_local_deficit,
)
from lif_tddft.campaign import expand_campaign, independent_histories, load_config, ramp_velocities
from lif_tddft.clusters import build_one_active_site_cluster, topology_check
from lif_tddft.models.capture_demkov import DemkovCapture
from lif_tddft.models.demkov_parameters import load_demkov_parameters
from lif_tddft.models.capture_energy import (
    dynamic_image_interaction,
    electrostatic_energy_defect,
    mott_littleton_polarization,
)
from lif_tddft.models.charge_state import propagate_ordered_events, transition_matrix
from lif_tddft.models.detachment_rt_tddft import RTDetachment, detachment_probability
from lif_tddft.octopus_input import render_input
from lif_tddft.restart import validate_restart
from lif_tddft.trajectory.events import identify_events
from lif_tddft.trajectory.integrate import integrate_grazing_trajectory
from lif_tddft.uncertainty import combine
from lif_tddft.units import angstrom_to_bohr, bohr_to_angstrom


class RevisionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = load_config(ROOT / "config" / "campaign_revision.yaml")
        # Explicit test fixtures; these are not published production parameters.
        cls.config['heights_bohr'] = [1.2, 2.0, 3.5, 10.0]
        cls.config['numerics']['cap_strength_magnitude_by_velocity'] = {
            str(v): 0.2 for v in cls.config['relative_velocities_au']}
        cls.config['pseudopotentials'].update(Li='test_Li.psf', F='test_F.psf')
        cls.specs = expand_campaign(cls.config)

    def test_units_round_trip(self):
        self.assertAlmostEqual(bohr_to_angstrom(angstrom_to_bohr(3.5)), 3.5, places=12)

    def test_campaign_is_independent(self):
        self.assertEqual(len(self.specs), 48)
        self.assertTrue(independent_histories(self.specs))
        self.assertEqual(ramp_velocities(0.10, 0.05, 0.025), [0.05, 0.075, 0.1])

    def test_input_contract(self):
        spec = self.specs[0]
        text = render_input(ROOT / "octopus" / "templates" / "single_flyby.inp.j2", spec, self.config)
        self.assertIn("MoveIons = yes", text)
        self.assertIn("IonsConstantVelocity = yes", text)
        self.assertIn("TDDynamics = ehrenfest", text)
        self.assertIn("| no", text)
        self.assertIn(" current", text)
        ramp = render_input(
            ROOT / "octopus" / "templates" / "single_flyby.inp.j2",
            spec,
            self.config,
            relative_velocity_au=0.05,
            propagation_time_expression="0.04*fs",
        )
        self.assertIn("TDPropagationTime = 0.04*fs", ramp)

    def test_isolated_projectile_template(self):
        spec = next(item for item in self.specs if item.variant == "isolated_projectile")
        text = render_input(
            ROOT / "octopus" / "templates" / "isolated_projectile.inp.j2",
            spec,
            self.config,
        )
        self.assertIn("Isolated moving F-", text)
        self.assertNotIn('"Li"', text)

    def test_ground_state_restart_does_not_require_td_coordinates(self):
        with tempfile.TemporaryDirectory() as directory:
            restart = Path(directory) / "restart"
            restart.mkdir()
            (restart / "states").write_text("ok", encoding="utf-8")
            self.assertTrue(validate_restart(directory, require_coordinates=False)["valid"])
            self.assertFalse(validate_restart(directory, require_coordinates=True)["valid"])

    def test_population_constant_density(self):
        axis = np.linspace(-1, 1, 21)
        points = np.asarray(np.meshgrid(axis, axis, axis, indexing="ij")).reshape(3, -1).T
        voxel = (axis[1] - axis[0]) ** 3
        value = moving_sphere_population(np.ones(points.shape[0]), points, np.zeros(3), 0.5, voxel, transition_half_width=0.1)
        self.assertAlmostEqual(value, 4 * np.pi * 0.5**3 / 3, delta=0.08)

    def test_paired_mean_loss_and_surface_correction(self):
        interacting = np.asarray([8.0, 7.7, 7.5])
        isolated = np.asarray([8.0, 7.9, 7.8])
        loss = paired_projectile_mean_loss(interacting, isolated)
        np.testing.assert_allclose(loss, [0.0, 0.2, 0.3])
        corrected = surface_corrected_local_deficit(loss, np.asarray([0.2, 0.25, 0.3]))
        np.testing.assert_allclose(corrected, [0.0, 0.25, 0.4])

    def test_flux_velocity_term_and_closure(self):
        normals = np.asarray([[1.0, 0.0, 0.0]])
        flux = outward_flux(np.zeros((1, 3)), np.ones(1), normals, np.ones(1), np.asarray([0.2, 0, 0]))
        self.assertAlmostEqual(flux, -0.2)
        result = closure_residual(np.arange(5.0), -2 * np.arange(5.0), 2 * np.ones(5))
        self.assertLess(result["relative_rms"], 1e-12)

    def test_overlap_descriptor(self):
        a = np.asarray([[1.0, 0.0]])
        b = np.asarray([[0.5, np.sqrt(0.75)]])
        s = cross_overlap(a, b, np.ones(2))
        self.assertAlmostEqual(occupied_overlap_descriptors(s)["Omega_occ"], 0.25)

    def test_bounded_interpolation(self):
        interp = DetachmentProbabilityInterpolator(np.asarray([3.0, 4.0]), np.asarray([0.2, 0.4]), np.asarray([[0.2, 0.4], [0.1, 0.2]]))
        self.assertAlmostEqual(interp(3.5, 0.2), 0.15)
        with self.assertRaises(OutOfDomainError):
            interp(3.5, 0.3)
        with self.assertRaises(OutOfDomainError):
            interp(2.0, 0.3)
        with self.assertRaisesRegex(ValueError, "clipping or rescaling is forbidden"):
            DetachmentProbabilityInterpolator(
                np.asarray([3.0, 4.0]),
                np.asarray([0.2, 0.4]),
                np.asarray([[0.2, 1.2], [0.1, 0.2]]),
            )

    def test_topology_gate(self):
        small = build_one_active_site_cluster("small")
        large = build_one_active_site_cluster("large")
        result = topology_check(small, large)
        self.assertTrue(result["passed"])
        self.assertEqual(len(small.atoms), 18)
        self.assertEqual(len(large.atoms), 50)
        self.assertEqual(sum(atom.species == "F" and atom.coordinates_bohr[2] == 0 for atom in large.atoms), 1)

    def test_events_are_ordered(self):
        times = np.arange(5.0)
        position = np.asarray([[-2, 0, 3], [-1, 0, 3], [0, 0, 3], [1, 0, 3], [2, 0, 3]], dtype=float)
        velocity = np.tile(np.asarray([1.0, 0.0, 0.0]), (5, 1))
        events = identify_events("t1", times, position, velocity, {"site": np.zeros(3)})
        self.assertEqual(events[0]["event_time"], 2.0)
        self.assertEqual(events[0]["event_index"], 0)

    def test_grazing_trajectory_keeps_parallel_speed_uniform(self):
        result = integrate_grazing_trajectory(
            0.0, 3.5, 0.3, -0.01, 10.0, 0.1, 5, lambda x, z: 0.2
        )
        np.testing.assert_allclose(np.diff(result["x"]), 0.03)
        np.testing.assert_allclose(result["v_parallel"], 0.3)
        self.assertGreater(result["v_normal"][-1], result["v_normal"][0])

    def test_transition_is_normalized(self):
        matrix = transition_matrix(0.2, 0.3)
        np.testing.assert_allclose(matrix.sum(axis=0), np.ones(2))
        self.assertEqual(matrix.shape, (2, 2))
        with self.assertRaisesRegex(ValueError, "clipping is forbidden"):
            transition_matrix(0.2, 1.01)

    def test_probability_data_cannot_be_scaled(self):
        policy = {
            "normalization": "none",
            "scale_factor": 1.0,
            "value_transform": "none",
            "velocity_specific_adjustments": "none",
        }
        assert_no_scaling_policy(policy)
        values = validate_probability_values(np.asarray([0.1, 0.4]), label="test")
        assert_values_identical(values, values.copy(), label="test")
        with self.assertRaisesRegex(ValueError, "scale_factor"):
            assert_no_scaling_policy(dict(policy, scale_factor=1.1))

    def test_single_electron_probability_domain(self):
        for value in (0.0, 0.25, 1.0):
            self.assertEqual(detachment_probability(value), value)
        for value in (-0.01, 1.01, np.nan, np.inf):
            with self.assertRaises(ValueError):
                detachment_probability(value)

    def test_detachment_probabilities_preserve_calculated_nodes(self):
        values = np.asarray([[0.6, 1.0], [0.2, 0.4]])
        model = RTDetachment.from_mean_loss_grid(np.asarray([2., 2.5]), np.asarray([.4, .5]), values)
        self.assertEqual(model.probability(dict(surface_height=2., v_parallel=.5)), 1.0)
        self.assertAlmostEqual(model.probability(dict(surface_height=2.25, v_parallel=.4)), 0.4)
        with self.assertRaises(OutOfDomainError):
            model.probability(dict(surface_height=1.9, v_parallel=.4))
        with self.assertRaises(ValueError):
            RTDetachment.from_mean_loss_grid(np.asarray([2., 2.5]), np.asarray([.4, .5]), values+0.1)

    def test_detachment_grid_does_not_fill_missing_nodes(self):
        for values in ([[0.1, np.nan], [0.3, 0.4]], [[0.1, 0.2], [0.3, np.inf]]):
            with self.assertRaisesRegex(ValueError, "finite"):
                RTDetachment.from_mean_loss_grid(np.asarray([2., 3.]), np.asarray([.1, .2]), np.asarray(values))

    def test_demkov_matches_capture_formula(self):
        gamma = 0.7  # Arithmetic test input, not a manuscript parameter.
        event = {"energy_defect_au": 0.2, "v_parallel": 0.1, "surface_height": 3.5}
        expected = 0.5 / np.cosh(
            np.pi * (0.2 + 0.1**2 / 2) / (2 * gamma * 0.1)
        ) ** 2
        model = DemkovCapture(gamma)
        self.assertAlmostEqual(model.probability(event), expected, places=15)
        changed_height = dict(event, surface_height=9.0)
        self.assertEqual(model.probability(event), model.probability(changed_height))
        stable_tail = model._sech_squared(350.0)
        expected_tail = 4.0 * np.exp(-700.0) / (1.0 + np.exp(-700.0)) ** 2
        self.assertEqual(stable_tail, expected_tail)
        self.assertGreater(stable_tail, 0.0)

    def test_demkov_requires_explicit_gamma_and_source(self):
        with self.assertRaisesRegex(ValueError,'gamma_c and gamma_source are required'):
            load_demkov_parameters(ROOT/'config/demkov_parameters.yaml')
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'capture.yaml'
            path.write_text('parameters: {gamma_bohr_inverse: 0.7}\ngamma_source: arithmetic test fixture\n',encoding='utf-8')
            parameters=load_demkov_parameters(path)
            self.assertEqual(parameters.gamma_bohr_inverse,0.7)
            self.assertEqual(parameters.source,'arithmetic test fixture')

    def test_capture_energy_equations(self):
        sites = np.asarray([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
        charges = np.asarray([1.0, -1.0])
        value = electrostatic_energy_defect(np.asarray([0.0, 0.0, 2.0]), sites, charges)
        self.assertAlmostEqual(value, 0.0)
        self.assertAlmostEqual(
            mott_littleton_polarization(np.zeros(3), sites, np.ones(2)), 0.0
        )
        image = dynamic_image_interaction(
            0.3, 3.0, np.geomspace(1.0e-8, 5.0, 20000), np.ones(20000)
        )
        self.assertAlmostEqual(image, -1.0 / (4.0 * 3.0), places=4)

    def test_two_state_recursion_matches_equation10(self):
        events = [dict(event_index=i, event_time=i, p_det=d, p_cap=c)
                  for i,(d,c) in enumerate([(0.2,0.3),(0.7,0.1),(1.0,0.4),(0.0,0.2)])]
        result = propagate_ordered_events(events, lambda e:e['p_det'], lambda e:e['p_cap'])
        fraction = 0.0
        for event,row in zip(events,result['history']):
            fraction = (1-fraction)*event['p_cap']+fraction*(1-event['p_det'])
            self.assertAlmostEqual(row['state'][0], fraction, places=15)
            self.assertEqual(len(row['state']), 2)
        self.assertEqual(set(result['states']), {'F_minus','F_neutral'})

    def test_naming_and_uncertainty_gates(self):
        self.assertEqual(allowed_label({}), "frozen-fragment orthogonalization descriptor")
        result = combine({"uncertainty_tddft_numerical": 0.03, "uncertainty_demkov_model": 0.04})
        self.assertAlmostEqual(result["quadrature_total"], 0.05)


if __name__ == "__main__":
    unittest.main()
