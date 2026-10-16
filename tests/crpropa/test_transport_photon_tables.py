#!/usr/bin/env python3
"""Preserve the pinned generator's conditional-yield edge semantics."""
import unittest
import numpy as np
from prepare_transport_photons import emission_ratio


class PhotonTablesTest(unittest.TestCase):
    def test_official_ratios_and_zero_support(self):
        ratio=emission_ratio(np.array([0.,1.,3.,1.]),np.array([0.,2.,2.,0.]))
        np.testing.assert_array_equal(ratio[:3],[0.,.5,1.5])
        self.assertTrue(np.isfinite(ratio[3]));self.assertGreater(ratio[3],1.)
        # CRPropa compares a uniform deviate to this ratio, saturating above one.
        self.assertEqual(min(ratio[2],1.),1.)


if __name__=='__main__':unittest.main()
