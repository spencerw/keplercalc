from unittest import TestCase
import numpy as np
import pynbody as pb
import os
from keplercalc import hello, kep2cart, cart2kep, orb_params

class TestKeplercalc(TestCase):
    def test_hello1(self):
        output = hello('Spencer')
        self.assertEqual(output, 'Hello Spencer!')

    def test_cart2kep2cart(self):
        tol = 1e-10

        # Earth orbiting Sun, slightly inclined so angles are defined
        m1, m2 = 1, 1e-20

        a = 1
        e = 0.05
        inc = 0.1
        asc_node = np.pi
        omega = np.pi
        M = np.pi

        X, Y, Z, vx, vy, vz = kep2cart(a, e, inc, asc_node, omega, M, m1, m2)

        self.assertTrue(np.fabs(X - -1.05) < tol)
        self.assertTrue(np.fabs(Y - 3.782338790704024e-16) < tol)
        self.assertTrue(np.fabs(Z - -2.5048146051777413e-17) < tol)
        self.assertTrue(np.fabs(vx - -3.490253699036788e-16) < tol)
        self.assertTrue(np.fabs(vy - -0.9464377445249709) < tol)
        self.assertTrue(np.fabs(vz - 0.09496052074620637) < tol)

        a, e, inc, asc_node, omega, M = cart2kep(X, Y, Z, vx, vy, vz, m1, m2)

        self.assertTrue(np.fabs(a - 1) < tol)
        self.assertTrue(np.fabs(e - 0.05) < tol)
        self.assertTrue(np.fabs(inc - 0.1) < tol)
        self.assertTrue(np.fabs(asc_node - np.pi) < tol)
        self.assertTrue(np.fabs(omega - np.pi) < tol)
        self.assertTrue(np.fabs(M - np.pi) < tol)

        # Now try converting back to cartesian
        X1, Y1, Z1, vx1, vy1, vz1 = kep2cart(a, e, inc, asc_node, omega, M, m1, m2)

        self.assertTrue(np.fabs(X1 - X) < tol)
        self.assertTrue(np.fabs(Y1 - Y) < tol)
        self.assertTrue(np.fabs(Z1 - Z) < tol)
        self.assertTrue(np.fabs(vx1 - vx) < tol)
        self.assertTrue(np.fabs(vy1 - vy) < tol)
        self.assertTrue(np.fabs(vz1 - vz) < tol)

    def test_snap(self):
        tol = 1e-8

        path = os.path.dirname(os.path.abspath(__file__))
        snap = pb.load(path + '/test.ic')
        pl = orb_params(snap)

        self.assertTrue(np.mean(pl['a']) - 2.96233536 < tol)
        self.assertTrue(np.mean(pl['e']) - 0.06809352 < tol)
        self.assertTrue(np.mean(pl['inc']) - 0.005003 < tol)
        self.assertTrue(np.mean(pl['asc_node']) - 3.15458748 < tol)
        self.assertTrue(np.mean(pl['omega']) - 3.13253456)
        self.assertTrue(np.mean(pl['M']) - 3.13972053 < tol)