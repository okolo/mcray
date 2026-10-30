"""Self-contained guards for the optional historical transport diagnostic."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parent/'fork_compat'))
import transport


def row():
    state = dict(position_Mpc=[0,0,0],direction=[0,0,1],flight_Mpc=.25,beta=1,
                 field_calls=5,scale_calls=1,zero_scale_calls=0,
                 direction_norm_error=0,time_error_Mpc=0,
                 analytic_position_error_Mpc=.01,analytic_direction_error=.02)
    return dict(case='free',kind=0,mode=0,particle=0,accuracy=10,outer_steps=1,
                difference_position_Mpc=.001,difference_direction=.002,difference_beta=0,
                fork=copy.deepcopy(state),modern=copy.deepcopy(state))


class TransportDiagnosticTests(unittest.TestCase):
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


if __name__=='__main__':
    unittest.main()
