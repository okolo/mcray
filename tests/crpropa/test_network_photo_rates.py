#!/usr/bin/env python3
"""Independent checks of the validation-only analytic photon integral."""
import unittest
import numpy as np
from scipy.integrate import quad
from diagnose_network_photo_rates import PhotonIntegral, interpolate


class PhotonIntegralTest(unittest.TestCase):
    def test_against_adaptive_quadrature(self):
        energy = np.array([1., 2., 4.])
        for values in ([2., 4., 1.], [0., 4., 0.], [4., 0., 3.], [1., 2., 4.]):
            density = np.array(values)
            integral = PhotonIntegral(energy, density)
            for x in (.5, 1., 1.1, 1.9, 2., 2.1, 3.9, 4., 5.):
                with self.subTest(density=values, lower=x):
                    expected = quad(lambda v: float(interpolate(np.asarray(v),energy,density))/v**2,
                        max(x, 1.), 4., points=[2.] if x<2 else None, epsabs=1e-12)[0] if x<4 else 0.
                    np.testing.assert_allclose(integral(np.array([x]))[0], expected, rtol=1e-11, atol=1e-13)

    def test_steep_tail_remains_finite(self):
        energy = np.geomspace(1e-7, 1., 4097)
        density = np.exp(-energy/1e-3)*energy
        result = PhotonIntegral(energy,density)(energy)
        self.assertTrue(np.all(np.isfinite(result)))
        self.assertTrue(np.all(result>=0))
        self.assertTrue(np.all(np.diff(result)<=0))
        self.assertEqual(result[-1],0.)


if __name__ == '__main__':
    unittest.main()
