"""Self-contained guards for the optional historical transport diagnostic."""
import copy
import math
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parent/'fork_compat'))
import transport
import observer
import source_observer_pilot as pilot
import source_observer_ensemble as ensemble


def row():
    state = dict(position_Mpc=[0,0,0],direction=[0,0,1],flight_Mpc=.25,beta=1,
                 field_calls=5,scale_calls=1,zero_scale_calls=0,
                 direction_norm_error=0,time_error_Mpc=0,
                 analytic_position_error_Mpc=.01,analytic_direction_error=.02)
    return dict(case='free',kind=0,mode=0,particle=0,accuracy=10,outer_steps=1,
                difference_position_Mpc=.001,difference_direction=.002,difference_beta=0,
                fork=copy.deepcopy(state),modern=copy.deepcopy(state))


class TransportDiagnosticTests(unittest.TestCase):
    def test_collection_deadline_kills_child_and_keeps_no_successful_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); (root/'tmp').mkdir()
            script=root/'wait.py'; script.write_text('import time; time.sleep(60)\n')
            collection=ensemble.Collection(root,Path(sys.executable),script)
            timer=threading.Timer(.1,collection.stop); timer.start()
            try:
                with self.assertRaisesRegex(RuntimeError,'Failed stream'):
                    collection.run(dict(name='timeout',count=1),'fork')
            finally:
                timer.cancel(); collection.stop()
            self.assertFalse((root/'timeout-fork.json').exists())
            self.assertTrue(all(process.poll() is not None for process in collection.processes))

    def test_final_design_and_fixture_reservations_are_disjoint(self):
        final,fixture=ensemble.design_for('final'),ensemble.design_for('fixture')
        self.assertEqual(final['total_histories'],48000000)
        self.assertEqual(final['family_size'],39)
        self.assertEqual([case['count'] for case in final['cases']],[8000000]*3)
        seeds=lambda d:{case[key] for case in d['cases'] for key in ('fork_seed','modern_seed')}
        self.assertFalse(seeds(final)&seeds(fixture))
        self.assertEqual(seeds(final),set(range(2004101,2004107)))

    def test_ensemble_intervals_keep_boundary_uncertainty_and_correct_two_sample_factor(self):
        design=ensemble.design_for('final')
        rows=self.ensemble_rows(design)
        result=ensemble.confidence_report(rows,design)
        self.assertEqual(len(result['differences']),33)
        self.assertEqual(len(result['hit_fractions']),6)
        self.assertTrue(result['compatible_at_reported_resolution'])
        difference=result['differences'][0]
        fraction=result['hit_fractions'][0]
        self.assertGreater(difference['half_width'],0)
        self.assertLess(difference['half_width'],.001)
        self.assertAlmostEqual(difference['half_width']/fraction['half_width'],math.sqrt(2))
        self.assertGreater(fraction['interval'][1],1) # No clipping at all-hit boundary.
        rows[1]['hits']=4000000; rows[1]['sums'][:5]=[4000000]*5
        self.assertFalse(ensemble.confidence_report(rows,design)['compatible_at_reported_resolution'])

    @staticmethod
    def ensemble_rows(design):
        return [dict(case=case['name'],implementation=impl,histories=case['count'],
                     failures=0,emitted_weight=case['count'],hits=case['count'],
                     sums=[case['count']]*5+[0]*6,past_end_histories=0)
                for case in design['cases'] for impl in ('fork','modern')]

    def test_ensemble_cannot_accept_partial_nonfinite_duplicate_or_low_precision_data(self):
        design=ensemble.design_for('final'); rows=self.ensemble_rows(design)
        for bad in (rows[:-1],rows+[rows[0]]):
            with self.assertRaisesRegex(ValueError,'Incomplete or duplicate'):
                ensemble.confidence_report(bad,design)
        broken=copy.deepcopy(rows); broken[0]['sums'][5]=float('nan')
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            ensemble.confidence_report(broken,design)
        broken=copy.deepcopy(rows); broken[0]['sums'][2]=-1
        with self.assertRaisesRegex(ValueError,'arrival counts'):
            ensemble.confidence_report(broken,design)
        fixture=ensemble.design_for('fixture')
        self.assertNotIn('differences',ensemble.confidence_report(self.ensemble_rows(fixture),fixture))
        fixture['mode']='final'
        with self.assertRaisesRegex(ValueError,'Insufficient'):
            ensemble.confidence_report(self.ensemble_rows(fixture),fixture)

    def test_cost_pilot_freezes_separate_seeds_counts_and_geometry(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            path,design=pilot.inputs(Path(a))
            self.assertEqual(design,pilot.inputs(Path(b))[1])
            self.assertEqual(design['total_histories'],6000)
            self.assertEqual([case['count'] for case in design['cases']],[1000]*3)
            seeds=[case[key] for case in design['cases'] for key in ('fork_seed','modern_seed')]
            self.assertEqual(len(set(seeds)),6)
            self.assertEqual(design['observer_radius_Mpc'],.01)
            self.assertEqual(design['observer_distance_Mpc'],.1)
            self.assertEqual(design['input_sha256'],transport.kernels.digest(path))

    def test_cost_pilot_rejects_partial_duplicate_nonfinite_and_missing_timings(self):
        design=dict(cases=[dict(name='tiny',count=2)],purpose='Cost only')
        rows=[]
        for implementation in ('fork','modern'):
            for i in range(2):
                rows.append(dict(record='history',case='tiny',implementation=implementation,
                    primary=i,weight=1,direction=[0,0,1],historical_steps=1,
                    flight=dict(hit=i==0,past_end=False,end_error_Mpc=0,position_Mpc=[0,0,.1])))
            rows.append(dict(record='timing',case='tiny',implementation=implementation,compute_seconds=.1))
        result=pilot.summarize(rows,design)
        self.assertEqual(result['histories'],4)
        self.assertEqual(result['groups'][0]['hits'],1)
        self.assertEqual(result['groups'][1]['emitted_weight'],2)
        self.assertNotIn('accepted',result)
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            pilot.summarize(rows[1:],design)
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            pilot.summarize(rows+[rows[0]],design)
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            pilot.summarize(rows[:-1],design)
        broken=copy.deepcopy(rows); broken[0]['weight']=float('nan')
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            pilot.summarize(broken,design)

    def test_observer_design_is_fixed_and_does_not_sample_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            cases,design=observer.inputs(Path(directory))
            self.assertEqual(len(cases.read_text().splitlines()),210)
            self.assertEqual(design['points'],174)
            self.assertEqual(design['flights'],36)
            self.assertFalse(design['source_sampling'])
            self.assertFalse(design['fields'])
            self.assertEqual(transport.kernels.digest(cases),design['sha256'])

    def test_observer_hash_mismatch_precedes_reference_extraction(self):
        with tempfile.TemporaryDirectory() as directory, \
             patch.object(observer.kernels,'digest',return_value='changed'), \
             patch.object(observer.transport,'extract') as extract:
            with self.assertRaisesRegex(ValueError,'pinned snapshot'):
                observer.extract(Path(directory),Path(directory))
            extract.assert_not_called()

    def test_observer_report_retains_disagreements_without_a_verdict(self):
        point=dict(case='surface',kind=0,fork_hit=False,modern_hit=True,
                   fork_exceeds_remaining=True,modern_exceeds_remaining=False)
        flight=dict(case='cap',kind=1,fork=dict(hit=True,past_end=True),
                    modern=dict(hit=True,past_end=False),variable_rate=dict(hit=True,past_end=False),
                    modern_paths_equal=False,difference_position_Mpc=.1,
                    difference_flight_Mpc=.1,difference_rate_path_time_Mpc=.01)
        result=observer.summarize([point,flight])
        self.assertEqual(result['membership_differences'],[point])
        self.assertEqual(result['ballistic_end_time_overshoots'],[flight])
        self.assertEqual(result['modern_rate_path_differences'],[flight])
        self.assertEqual(result['modern_steps_exceeding_remaining'],0)
        self.assertNotIn('accepted',result)

    def test_shared_compiler_helper_preserves_original_default_sources(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            build=root/'build'
            recipe=build/'tests/CMakeFiles/grid_magnetic_field_tests.dir/link.txt'
            recipe.parent.mkdir(parents=True)
            recipe.write_text('c++ -O2 test.o -o tests ../mcray/libmcray.a -lgsl\n')
            output=root/'out'; output.mkdir()
            with patch.object(transport.kernels.subprocess,'check_output',return_value='compiler'), \
                 patch.object(transport.kernels.subprocess,'run') as compile:
                transport.kernels.compile_probe(build,output,[])
                command=compile.call_args.args[0]
                self.assertIn(str(output/'reference.cpp'),command)
                self.assertIn(str(transport.HERE/'compare.cpp'),command)
                self.assertIn(str(output/'fork-kernel-probe'),command)
                transport.kernels.compile_probe(build,output,[],[output/'custom.cpp'],'custom')
                command=compile.call_args.args[0]
                self.assertIn(str(output/'custom.cpp'),command)
                self.assertNotIn(str(output/'reference.cpp'),command)

    def test_frozen_design_covers_species_boundaries_and_refinements(self):
        with tempfile.TemporaryDirectory() as directory:
            cases,grid,design = transport.inputs(Path(directory))
            rows = [line.split() for line in cases.read_text().splitlines()]
            self.assertEqual(len(rows),756)
            self.assertEqual(len({tuple(line[:6]) for line in rows}),756)
            self.assertEqual({line[3] for line in rows},{'0','1','2'})
            self.assertEqual({line[4] for line in rows},{'10','20','40'})
            self.assertEqual({line[5] for line in rows},{'1','4','16'})
            grids = [line for line in rows if line[1]=='3']
            self.assertEqual({line[2] for line in grids},{'0','1','2','3'})
            self.assertEqual(len({line[0] for line in grids}),6)
            self.assertFalse(design['ensembles'])
            self.assertEqual(len(grid.read_text().splitlines()),3*4*5)
            for name,digest in design['sha256'].items():
                self.assertEqual(transport.kernels.digest(Path(directory)/name),digest)

    def test_inputs_do_not_depend_on_output_path(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            self.assertEqual(transport.inputs(Path(a))[2],transport.inputs(Path(b))[2])

    def test_reference_mismatch_fails_before_source_extraction(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(transport.kernels,'digest',return_value='changed'), \
                 patch.object(transport.kernels,'extract') as extract:
                with self.assertRaisesRegex(ValueError,'pinned snapshot'):
                    transport.extract(Path(directory),Path(directory))
                extract.assert_not_called()

    def test_summarize_keeps_every_refinement_and_reports_no_verdict(self):
        rows = [row(), row() | dict(accuracy=20,difference_position_Mpc=.0005)]
        summary = transport.summarize(rows)
        self.assertEqual(summary['cases'],2)
        self.assertEqual(len(summary['convergence']),2)
        self.assertEqual(summary['max_position_difference_Mpc'],.001)
        self.assertNotIn('accepted',summary)
        self.assertEqual(summary['convergence'][0]['fork_max_analytic_position_error_Mpc'],.01)

    def test_duplicate_or_nonfinite_trajectories_cannot_form_a_report(self):
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            transport.summarize([row(),row()])
        broken = row()
        broken['modern']['position_Mpc'][1] = float('nan')
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            transport.summarize([broken])

    def test_constant_field_orientations_have_separate_summary_groups(self):
        axial=row() | dict(case='axial',kind=1,difference_position_Mpc=.001)
        oblique=row() | dict(case='oblique',kind=1,difference_position_Mpc=.002)
        summary=transport.summarize([axial,oblique])
        groups={group['case']:group for group in summary['convergence']}
        self.assertEqual(set(groups),{'axial','oblique'})
        self.assertEqual(groups['axial']['max_position_difference_Mpc'],.001)
        self.assertEqual(groups['oblique']['max_position_difference_Mpc'],.002)
        self.assertEqual(summary['max_position_difference_Mpc'],.002)


if __name__=='__main__':
    unittest.main()
