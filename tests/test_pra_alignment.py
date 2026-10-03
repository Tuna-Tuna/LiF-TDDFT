"""Independent numerical and workflow checks against the final PRA definitions."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from lif_tddft.analysis.spherical_weight import spherical_weight
from lif_tddft.analysis.population import moving_sphere_population
from lif_tddft.analysis.moving_flux import smooth_outward_flux, closure_residual
from lif_tddft.analysis.interpolation import DetachmentProbabilityInterpolator, OutOfDomainError
from lif_tddft.campaign import load_config, expand_campaign
from lif_tddft.octopus_input import render_input
from lif_tddft.units import angstrom_to_bohr

ROOT=Path(__file__).resolve().parents[1]


class PRAAlignmentTests(unittest.TestCase):
    def config(self):
        config=load_config(ROOT/'config/campaign_revision.yaml')
        config['heights_bohr']=[1.2,3.5,10.0]
        config['numerics']['cap_strength_magnitude_by_velocity']={str(v):0.2 for v in config['relative_velocities_au']}
        config['pseudopotentials'].update(Li='test_Li.psf',F='test_F.psf')
        return config

    def test_missing_parameters_are_not_estimated(self):
        with self.assertRaisesRegex(ValueError, 'heights_bohr is required'):
            expand_campaign(load_config(ROOT/'config/campaign_revision.yaml'))
        config=self.config()
        spec=expand_campaign(config)[0]
        template=ROOT/'octopus/templates/single_flyby.inp.j2'
        config['numerics']['cap_strength_magnitude_by_velocity']['0.1']=None
        with self.assertRaisesRegex(ValueError, 'CAP magnitude is required'):
            render_input(template,spec,config)
        config=self.config()
        config['pseudopotentials']['F']=None
        with self.assertRaisesRegex(ValueError,'pseudopotential filename is required'):
            render_input(template,spec,config)

    def test_generated_stages_and_physical_settings(self):
        module_spec=importlib.util.spec_from_file_location('generator',ROOT/'octopus/generate_inp_files.py')
        module=importlib.util.module_from_spec(module_spec)
        module_spec.loader.exec_module(module)
        config=self.config()
        spec=expand_campaign(config)[0]
        with tempfile.TemporaryDirectory() as directory:
            inp=Path(directory)/'inp'
            inp.write_text('test',encoding='utf-8')
            manifest=module.build_manifest(spec,config,inp)
        self.assertEqual(manifest['stages'][-1]['name'],'production_constant_velocity')
        self.assertTrue(all(stage['ions_constant_velocity'] for stage in manifest['stages'][1:]))
        text=render_input(ROOT/'octopus/templates/single_flyby.inp.j2',spec,config)
        self.assertIn('IonsConstantVelocity = yes',text)
        self.assertIn('dt = 0.0004*fs',text)
        self.assertIn(str(angstrom_to_bohr(9.0)),text)
        self.assertIn('"F"    | 0*Ha |  0*Ha | 0+zcoor',text)
        with self.assertRaisesRegex(ValueError,'constant velocity'):
            render_input(ROOT/'octopus/templates/single_flyby.inp.j2',spec,config,constant_velocity=False)

    def test_smooth_weight_and_analytical_gradient(self):
        points=np.array([[0,0,0],[3.2,0,0],[3.5,0,0],[3.8,0,0],[4,0,0]])
        w,g=spherical_weight(points,np.zeros(3),3.5,0.3)
        np.testing.assert_allclose(w,[1,1,.5,0,0],atol=1e-14)
        eps=1e-6
        plus=points.copy();minus=points.copy();plus[:,0]+=eps;minus[:,0]-=eps
        derivative=(spherical_weight(plus,np.zeros(3),3.5,.3)[0]-spherical_weight(minus,np.zeros(3),3.5,.3)[0])/(2*eps)
        self.assertAlmostEqual(g[2,0],derivative[2],places=8)

    def test_generated_plan_passes_runner_dry_run(self):
        modules=[]
        for name,path in [('generate_pra','octopus/generate_inp_files.py'),
                          ('run_pra','python/run_octopus_workflow.py')]:
            spec=importlib.util.spec_from_file_location(name,ROOT/path)
            module=importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            modules.append(module)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            config=self.config()
            config['relative_velocities_au']=[0.1]
            path=root/'campaign.json'
            path.write_text(json.dumps(config),encoding='utf-8')
            rows=modules[0].generate(path,root/'runs')
            self.assertEqual(len(rows),6)
            plan=Path(rows[0]['manifest'])
            modules[1].run_plan(plan,path,'not-launched',[],dry_run=True)
            status=json.loads((plan.parent/'run_status.json').read_text(encoding='utf-8'))
            self.assertEqual(status['state'],'validated')
            self.assertTrue(all(stage['state']=='validated' for stage in status['stages']))

    def test_smooth_population_matches_radial_integral(self):
        axis=np.linspace(-.7,.7,71)
        points=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1).reshape(-1,3)
        volume=(axis[1]-axis[0])**3
        actual=moving_sphere_population(np.ones(len(points)),points,np.zeros(3),.5,volume,transition_half_width=.1)
        expected=4*np.pi*(.5**3/3+.5*.1**2*(1-8/np.pi**2))
        self.assertAlmostEqual(actual,expected,delta=1e-4)
        # A rigidly co-moving density has zero transport through its moving weight.
        velocity=np.array([.2,0,0])
        flux=smooth_outward_flux(np.tile(velocity,(len(points),1)),np.ones(len(points)),points,np.zeros(3),velocity,.5,.1,volume)
        self.assertEqual(flux,0.0)
        check=closure_residual(np.arange(5.),-3*np.arange(5.),2*np.ones(5),absorption_rates=np.ones(5))
        self.assertEqual(check['relative_rms'],0.0)

    def test_pchip_is_shape_preserving_and_height_only(self):
        h=np.array([1.2,2,3.5,10.])
        values=np.array([[1.],[.8],[.1],[0.]])
        interp=DetachmentProbabilityInterpolator(h,np.array([.2]),values)
        sampled=np.array([interp(x,.2) for x in np.linspace(1.2,10,101)])
        self.assertTrue(np.all(np.diff(sampled)<=1e-14))
        self.assertTrue(np.all((sampled>=0)&(sampled<=1)))
        for height,velocity in [(1.19,.2),(10.01,.2),(2.,.3)]:
            with self.assertRaises(OutOfDomainError): interp(height,velocity)
        for x,y in zip(h,values[:,0]): self.assertEqual(interp(x,.2),y)

if __name__=='__main__': unittest.main()
