import unittest
from secondary_statistics import holm, interval, composite_p, analyze
class ExactStatistics(unittest.TestCase):
    def test_holm_preserves_original_order_and_family(self):
        self.assertEqual(holm([.5,.001,.02]),[.5,.003,.04])
    def test_zero_counts_have_positive_upper_limit(self):
        lo,hi=interval(0,1000,.001)
        self.assertEqual(lo,0)
        self.assertGreater(hi,.007)
    def test_reference_uncertainty_is_conservative(self):
        self.assertEqual(composite_p(1,1000,0,.01),1)
        self.assertLess(composite_p(50,1000,.009,.011),1e-10)
    def test_pair_trials_are_events_not_particles(self):
        row=dict(p=[.5,0,.5],error=[0,0,0],mean=1.,mean_error=0.,possible=[1,0,1],target=.1)
        counts={'case':{m:{'all':[500,0,500]} for m in ('old','current')}}
        result=analyze({'case':{'all':row}},counts,{'case':1000},True)
        self.assertEqual(result['family_size'],9)
        self.assertEqual(result['spectra'][0]['N'],1000)
        self.assertEqual(result['spectra'][0]['mean'],1)
    def test_incomplete_run_cannot_pass_or_claim_significance(self):
        row=dict(p=[.9,.1,0],error=[0,0,0],mean=.1,mean_error=0.,possible=[1,1,0],target=.1)
        counts={'case':{m:{'all':[500,500,0]} for m in ('old','current')}}
        result=analyze({'case':{'all':row}},counts,{'case':1000},False)
        self.assertEqual(result['status'],'inconclusive')
        self.assertFalse(any(t['significant'] for t in result['tests']))

class AcceptanceCriteria(unittest.TestCase):
    def test_detects_large_tail_bias_and_direct_difference(self):
        row=dict(p=[.999,.001,0],error=[1e-12,1e-12,1e-12],mean=.001,mean_error=3e-12,possible=[1,1,0],target=.1)
        counts={'case':{'old':{'tail':[998000,2000,0]},'current':{'tail':[999000,1000,0]}}}
        result=analyze({'case':{'tail':row}},counts,{'case':1000000},True)
        self.assertEqual(result['status'],'deviation')
        self.assertTrue(any(t['significant'] and t['comparison']=='old_current' for t in result['tests']))
        self.assertEqual(result['spectra'][0]['status'],'deviation')
    def test_unobserved_possible_tail_cannot_establish_equality(self):
        row=dict(p=[1-1e-9,1e-9,0],error=[1e-14]*3,mean=1e-9,mean_error=3e-14,possible=[1,1,0],target=.1)
        counts={'case':{m:{'tail':[1000,0,0]} for m in ('old','current')}}
        result=analyze({'case':{'tail':row}},counts,{'case':1000},True)
        self.assertEqual(result['status'],'inconclusive')
        self.assertGreater(result['spectra'][0]['mean_ci'][1],0)

    def test_wrong_ensemble_size_cannot_claim_significance(self):
        row=dict(p=[.9,.1,0],error=[0,0,0],mean=.1,mean_error=0.,possible=[1,1,0],target=.1)
        counts={'case':{m:{'all':[500,500,0]} for m in ('old','current')}}
        result=analyze({'case':{'all':row}},counts,{'case':10000},True)
        self.assertFalse(result['complete'])
        self.assertEqual(result['status'],'inconclusive')
        self.assertFalse(any(t['significant'] for t in result['tests']))

if __name__=='__main__':unittest.main()
