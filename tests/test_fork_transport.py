"""Self-contained guards for the optional historical transport diagnostic."""
import copy
from decimal import Decimal
import math
import signal
from pathlib import Path
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parent/'fork_compat'))
import transport
import observer
import source_observer_pilot as pilot
import source_observer_ensemble as ensemble
import vmf as vmf_probe
import coupled
import field_source_pilot as field_pilot
import eev_field_ensemble as eev
import run_eev_field_server as server
import subprocess
import tarfile
import io
import os


def row():
    state = dict(position_Mpc=[0,0,0],direction=[0,0,1],flight_Mpc=.25,beta=1,
                 field_calls=5,scale_calls=1,zero_scale_calls=0,
                 direction_norm_error=0,time_error_Mpc=0,
                 analytic_position_error_Mpc=.01,analytic_direction_error=.02)
    return dict(case='free',kind=0,mode=0,particle=0,accuracy=10,outer_steps=1,
                difference_position_Mpc=.001,difference_direction=.002,difference_beta=0,
                fork=copy.deepcopy(state),modern=copy.deepcopy(state))


class TransportDiagnosticTests(unittest.TestCase):
    def test_eev_final_family_and_disjoint_reservations(self):
        final,fixture=eev.design_for('final'),eev.design_for('fixture')
        self.assertEqual(final['total_histories'],162000000)
        self.assertEqual(fixture['total_histories'],18000)
        self.assertEqual(final['accuracy'],40)
        seeds=lambda d:{c[k] for c in d['cases'] for k in ('fork_seed','modern_seed')}
        self.assertEqual(seeds(final),set(range(2010101,2010119)))
        self.assertFalse(seeds(final)&seeds(fixture))
        rows=[dict(case=c['name'],implementation=impl,histories=c['count'],failures=0,
                   emitted_weight=c['count'],hits=c['count'],past_end_histories=0,
                   max_engine_moves=1,sums=[c['count']]*5+[0]*9)
              for c in final['cases'] for impl in ('fork','modern')]
        report=eev.confidence_report(rows,final)
        self.assertEqual(len(report['differences']),126)
        self.assertEqual(len(report['hit_fractions']),18)
        self.assertAlmostEqual(report['differences'][0]['half_width'],.000980855225113378)
        self.assertTrue(report['compatible_at_reported_resolution'])
        # A final-direction-only discrepancy must not be lost among source moments.
        rows[1]['sums'][8]=rows[1]['histories']
        self.assertFalse(eev.confidence_report(rows,final)['compatible_at_reported_resolution'])
        with self.assertRaises(ValueError): eev.confidence_report(rows[:-1],final)
        broken=copy.deepcopy(rows); broken[0]['sums'][11]=float('nan')
        with self.assertRaises(ValueError): eev.confidence_report(broken,final)
        altered=copy.deepcopy(final); altered['alpha']=.1
        with self.assertRaisesRegex(ValueError,'Frozen'): eev.confidence_report(rows,altered)

    def test_final_ledger_is_exclusive(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'ledger.json'
            eev.reserve(path,{'revision':'first'})
            with self.assertRaises(FileExistsError): eev.reserve(path,{'revision':'second'})
            self.assertIn('first',path.read_text())

    def test_server_fixture_workflow_has_no_final_collection(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); archive=root/'input.tar.gz'
            with tarfile.open(archive,'w:gz') as out:
                member=tarfile.TarInfo('snapshot/file'); member.size=1
                out.addfile(member,io.BytesIO(b'x'))
            revision='a'*40; calls=[]
            def execute(supervisor,name,command,cwd,env):
                calls.append((name,command))
                if name=='clone':
                    candidate=Path(command[-1])/'tests/fork_compat'
                    candidate.mkdir(parents=True)
                    (candidate/'run_eev_field_server.py').write_bytes(Path(server.__file__).read_bytes())
            args=SimpleNamespace(data_dir=root/'results',revision=revision,repository='candidate.bundle',
                workers=2,fixture_only=True,reference_archive=archive,gtest_archive=archive)
            with patch.object(server.platform,'platform',return_value='workflow-fixture'),patch.object(server.Supervisor,'execute',execute),patch.object(server,'GTEST_SHA256',server.sha256(archive)),patch.object(server.subprocess,'run',return_value=SimpleNamespace(returncode=0)),patch.object(server.subprocess,'check_output',return_value=revision+'\n'):
                self.assertEqual(server.run(args),0)
            self.assertEqual([name for name,_ in calls],['clone','checkout','versions','configure','build','ctest','python-guards','fixture'])
            self.assertFalse((args.data_dir/'eev-field-ensemble-v1.reservation.json').exists())
            self.assertEqual(server.LIMIT,28800)

    def test_server_archive_rejects_escape_and_links(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for index,name in enumerate(('../outside','root/link')):
                archive=root/(str(index)+'.tar.gz')
                with tarfile.open(archive,'w:gz') as out:
                    member=tarfile.TarInfo(name)
                    if index: member.type=tarfile.SYMTYPE; member.linkname='/etc/passwd'
                    out.addfile(member)
                with self.assertRaises(ValueError): server.extract_archive(archive,root/('extract'+str(index)))

    @unittest.skipUnless(sys.platform.startswith('linux'),'Linux process-tree watchdog')
    def test_server_deadline_kills_worker_in_separate_session(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); pidfile=root/'worker.pid'
            script=root/'parent.py'
            script.write_text(
                "import subprocess, sys, threading, time\n"
                "def launch():\n"
                "    worker = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'], start_new_session=True)\n"
                "    with open(sys.argv[1], 'w') as output:\n"
                "        output.write(str(worker.pid))\n"
                "    time.sleep(60)\n"
                "thread = threading.Thread(target=launch)\n"
                "thread.start()\n"
                "thread.join()\n"
                "time.sleep(60)\n"
            )
            supervisor = server.Supervisor(root, seconds=2)
            with self.assertRaises(subprocess.TimeoutExpired):
                supervisor.execute('deadline', [sys.executable, str(script), str(pidfile)], root, os.environ.copy())
            pid = int(pidfile.read_text())
            status = Path('/proc') / str(pid) / 'stat'
            # SIGKILL delivery is asynchronous for a non-child PID. It can
            # briefly remain visible before exiting or becoming a zombie.
            deadline = time.monotonic() + 2
            try:
                while True:
                    try:
                        state = status.read_text().rsplit(')', 1)[1].split()[0]
                    except (FileNotFoundError, ProcessLookupError):
                        break
                    if state == 'Z':
                        break
                    self.assertLess(time.monotonic(), deadline, 'Worker survived the watchdog')
                    time.sleep(.01)
            except BaseException:
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                raise

    def test_field_pilot_reservations_are_separate_and_match_revised_model(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            cases,grid,design=field_pilot.inputs(Path(a),'pilot')
            _,_,fixture=field_pilot.inputs(Path(b),'fixture')
            self.assertEqual(design['particles'],['electron','photon'])
            self.assertEqual(design['energy_eV'],1e15)
            self.assertEqual(design['total_histories'],12000)
            self.assertEqual(fixture['total_histories'],24)
            self.assertEqual({c['species'] for c in design['cases']},{0,3})
            self.assertEqual({c['accuracy'] for c in design['cases']},{20})
            seeds=lambda d:{c[k] for c in d['cases'] for k in ('fork_seed','modern_seed')}
            self.assertEqual(seeds(design),set(range(2005101,2005113)))
            self.assertFalse(seeds(design)&seeds(fixture))
            self.assertFalse(seeds(design)&set(range(2003101,2003107)))
            self.assertFalse(seeds(design)&set(range(2004101,2004107)))
            self.assertEqual(field_pilot.kernels.digest(grid),design['grid_sha256'])
            self.assertEqual(len(cases.read_text().splitlines()),6)

    def test_eev_pilot_preserves_old_design_and_reserves_new_streams(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b,tempfile.TemporaryDirectory() as c:
            cases,grid,design=field_pilot.inputs(Path(a),'pilot','eev')
            _,_,fixture=field_pilot.inputs(Path(b),'fixture','eev')
            old_cases,old_grid,old=field_pilot.inputs(Path(c),'pilot')
            self.assertEqual(design['energy_eV'],1e18)
            self.assertEqual(design['particles'],['electron','photon','proton'])
            self.assertEqual(design['total_histories'],18000)
            self.assertEqual(fixture['total_histories'],36)
            self.assertEqual(len(cases.read_text().splitlines()),9)
            self.assertEqual(grid.read_bytes(),old_grid.read_bytes())
            seeds=lambda d:{v[k] for v in d['cases'] for k in ('fork_seed','modern_seed')}
            self.assertEqual(seeds(design),set(range(2007101,2007119)))
            self.assertEqual(seeds(fixture),set(range(2008101,2008119)))
            self.assertFalse(seeds(design)&(seeds(old)|seeds(fixture)))
            self.assertEqual(old['energy_eV'],1e15)
            self.assertEqual(old['total_histories'],12000)
            self.assertEqual(len(old_cases.read_text().splitlines()),6)

    def test_field_pilot_rejects_incomplete_and_invalid_final_states(self):
        design=dict(cases=[dict(name='electron',count=2)],purpose='Cost only',mode='fixture')
        rows=[]
        for impl in ('fork','modern'):
            for i in range(2):
                rows.append(dict(record='history',case='electron',implementation=impl,primary=i,
                    weight=1,direction=[0,0,1],historical_steps=2,
                    flight=dict(hit=False,past_end=False,end_error_Mpc=0,position_Mpc=[0,0,.2],
                                direction=[0,0,1],moves=2,field_calls=8,zero_scale_calls=0)))
            rows.append(dict(record='timing',case='electron',implementation=impl,compute_seconds=1))
        report=field_pilot.summarize(rows,design)
        self.assertEqual(report['histories'],4)
        self.assertEqual(report['groups'][0]['max_engine_moves'],2)
        self.assertNotIn('compatible_at_reported_resolution',report)
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            field_pilot.summarize(rows[:-1],design)
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            field_pilot.summarize(rows+[rows[0]],design)
        for bad in ([0,0], [0,0,float('nan')]):
            broken=copy.deepcopy(rows); broken[0]['flight']['direction']=bad
            with self.assertRaises(ValueError): field_pilot.summarize(broken,design)
        broken=copy.deepcopy(rows); broken[0]['weight']=.5
        with self.assertRaises(ValueError): field_pilot.summarize(broken,design)

    def test_requested_eev_species_matrix_keeps_energy_and_refinements(self):
        with tempfile.TemporaryDirectory() as directory:
            cases,_,design=coupled.inputs(Path(directory),('electron','photon','proton'))
            self.assertEqual(design['species'],['electron','photon','proton'])
            self.assertEqual(design['energy_eV'],1e18)
            rows=[line.split() for line in cases.read_text().splitlines()]
            self.assertEqual(len(rows),729)
            self.assertEqual({int(row[3]) for row in rows},{0,2,3})
            coverage={tuple(row[:4]):set() for row in rows}
            for row in rows: coverage[tuple(row[:4])].add(int(row[4]))
            self.assertTrue(all(values=={10,20,40} for values in coverage.values()))
        for species in ((),('electron','electron'),('muon',)):
            with self.assertRaises(ValueError): coupled.inputs(Path('.'),species)

    def test_coupled_matrix_reuses_existing_fields_and_freezes_all_refinements(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            cases,grid,design=coupled.inputs(Path(a))
            other,other_grid,other_design=coupled.inputs(Path(b))
            self.assertEqual(design,other_design)
            self.assertEqual(design['cases'],729)
            self.assertEqual(design['initial_directions'],27)
            self.assertEqual(cases.read_bytes(),other.read_bytes())
            self.assertEqual(grid.read_bytes(),other_grid.read_bytes())
            self.assertEqual(transport.kernels.digest(grid),design['reused_grid_sha256'])
            self.assertFalse(design['source_sampling'])
            rows=[line.split() for line in cases.read_text().splitlines()]
            self.assertEqual(len({tuple(row[:5]) for row in rows}),729)
            self.assertEqual({row[4] for row in rows},{'10','20','40'})
            self.assertEqual(sum(row[1].startswith('tangent') for row in rows),216)

    def test_coupled_report_preserves_disagreements_and_decimal_refinement(self):
        rows=[]
        for accuracy in (10,20,40):
            state=dict(hit=True,past_end=False,zero_scale_calls=0,
                       position_Mpc=[Decimal('0.1')+Decimal(accuracy)*Decimal('1e-23'),0,0],
                       direction=[0,0,1],flight_Mpc=Decimal('.1'))
            rows.append(dict(profile='cone',ray='tangent_1',field=2,particle=0,accuracy=accuracy,
                fork=copy.deepcopy(state),modern=copy.deepcopy(state),same_hit=accuracy!=10,
                difference_position_Mpc=Decimal('1e-20'),difference_direction=0,difference_time_Mpc=0))
        rows[0]['modern']['hit']=False
        result=coupled.summarize(rows)
        self.assertEqual(len(result['hit_differences']),1)
        self.assertEqual(len(result['groups']),3)
        self.assertGreater(result['refinements'][0]['fork']['position_change_Mpc'],0)
        self.assertNotIn('accepted',result)
        with self.assertRaisesRegex(ValueError,'Duplicate'):
            coupled.summarize(rows+[rows[0]])
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            coupled.summarize(rows[:-1])
        broken=copy.deepcopy(rows); broken[0]['modern']['flight_Mpc']=float('nan')
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            coupled.summarize(broken)

    def test_vmf_point_design_separates_ordinary_and_numerical_controls(self):
        with tempfile.TemporaryDirectory() as a,tempfile.TemporaryDirectory() as b:
            points,count=vmf_probe.inputs(Path(a))
            other,other_count=vmf_probe.inputs(Path(b))
            self.assertEqual(count,87); self.assertEqual(other_count,count)
            self.assertEqual(points.read_bytes(),other.read_bytes())
            rows=[line.split() for line in points.read_text().splitlines()]
            self.assertEqual(sum(row[0]=='ordinary' for row in rows),70)
            self.assertEqual(sum(row[0]=='cancellation' for row in rows),12)
            self.assertEqual(sum(row[0]=='range' for row in rows),5)
            self.assertIn(-100.,{float(row[3]) for row in rows})
            self.assertIn(1e-14,{float(row[3]) for row in rows})

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
