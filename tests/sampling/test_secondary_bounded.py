import copy
import tempfile
from pathlib import Path
import unittest
from analyze_bounded_secondary import analyze_run
import secondary_statistics as stats


class BoundedFamilies(unittest.TestCase):
    def fixture(self):
        row=dict(p=[.5,.5,0.],error=[0.,0.,0.],mean=.5,mean_error=0.,target=.1,possible=[True,True,False])
        counts={'old':{'tail':[500,500,0]},'cdf':{'tail':[500,500,0]},'current':{'tail':[500,500,0]}}
        return dict(references={'ICS':{'tail':row}},counts={'ICS':counts},planned_n={'ICS':1000},complete=True,
                    controls=['old','cdf'],sampler_labels={'old':'Rejection','cdf':'CDF','current':'Bounded'})

    def test_global_level_is_split_and_module_default_restored(self):
        with tempfile.TemporaryDirectory() as directory:
            result=analyze_run(self.fixture(),Path(directory)/'report.json')
        self.assertEqual(result['alpha_allocation'],{'old':.005,'cdf':.005})
        self.assertEqual(result['family_size'],18)
        self.assertEqual(stats.ALPHA,.01)
        self.assertFalse(any(t['significant'] for t in result['tests']))

    def test_deviation_in_either_control_survives_combination(self):
        data=self.fixture();data['counts']['ICS']['cdf']['tail']=[0,0,1000]
        with tempfile.TemporaryDirectory() as directory:
            result=analyze_run(data,Path(directory)/'report.json')
        self.assertEqual(result['status'],'deviation')
        self.assertTrue(any(t['significant'] for t in result['tests'] if t['comparison_family']=='cdf'))

    def test_incomplete_runs_cannot_establish_equivalence(self):
        data=self.fixture();data['complete']=False
        with tempfile.TemporaryDirectory() as directory:
            result=analyze_run(data,Path(directory)/'report.json')
        self.assertEqual(result['status'],'inconclusive')
        self.assertFalse(any(t['significant'] for t in result['tests']))

    def test_single_control_retains_original_level_and_family(self):
        data=self.fixture();data['controls']=['old']
        with tempfile.TemporaryDirectory() as directory:
            result=analyze_run(data,Path(directory)/'report.json')
        self.assertEqual(result['alpha_allocation'],{'old':.01})
        self.assertEqual(result['family_size'],9)

    def test_unknown_and_duplicate_control_families_are_rejected(self):
        for controls in ([],['cdf','cdf'],['unknown']):
            data=copy.deepcopy(self.fixture());data['controls']=controls
            with tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):analyze_run(data,Path(directory)/'report.json')
