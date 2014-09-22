__author__ = 'taowu'

import unittest
from fmin import fmin
from scipy.optimize import minimize_scalar

class TestFmin(unittest.TestCase):
    def test_fmin(self):
        """
        Test the fmin function which finds the x that gives the minimal output for f(x)
        """
        from math import pi
        fx = lambda x:(x - pi) ** 2 + 2 ** 0.5
        self.assertAlmostEqual(fmin(-10.0, 10.0, fx, 0), pi)
        # compare result of fmin and minimize_scalar
        self.assertAlmostEqual(fmin(-10.0, 10.0, fx, 0), minimize_scalar(fx, bounds=(-10.0, 10.0), method='bounded').x)

