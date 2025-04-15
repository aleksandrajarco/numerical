import argparse
import unittest
from unittest.mock import patch, mock_open
import tempfile
import os
from math import cos
import builtins

# Assume this imports everything from the module
import sys

sys.modules["matplotlib.pyplot"] = __import__(
    "matplotlib.pyplot"
)  # Avoid errors if matplotlib not in test env
from monteCarlo import MonteCarlo, plot_results


class TestMonteCarlo(unittest.TestCase):
    def setUp(self):
        # Default mocked args
        self.mock_args = argparse.Namespace(
            xrange=[0, 1], yrange=[0, 1], nsteps=1000, stepsize=100
        )

        self.mc = MonteCarlo(self.mock_args)

    def test_monte_carlo_integration_temp_file_exists(self):
        result_file = self.mc.monte_carlo_integration()
        self.assertTrue(result_file.startswith("/tmp"))
        self.assertTrue(os.path.exists(result_file))

    def test_monte_carlo_integration_results(self):
        result_file = self.mc.monte_carlo_integration()
        with open(result_file, "r") as f:
            lines = f.readlines()
            stepsize = self.mock_args.stepsize
            nsteps = self.mock_args.nsteps
            xrange_min = self.mock_args.xrange[0]
            xrange_max = self.mock_args.xrange[1]
            expected_len = (
                self.mock_args.nsteps // stepsize - 1
            )  # nsteps / stepsize - 1
            self.assertEqual(len(lines), expected_len)
            for line in lines:
                n, integral = map(float, line.strip().split())
                self.assertGreaterEqual(n, stepsize)
                self.assertLessEqual(n, nsteps)
                self.assertGreaterEqual(integral, xrange_min)
                self.assertLessEqual(integral, xrange_max)
                self.assertLessEqual(integral, 1)

    def test_plot_results_plot_exists(self):
        plot_file = self.mc.monte_carlo_integration()
        plot_results(plot_file, "integrationResultPlot.png")
        self.assertTrue(os.path.exists("integrationResultPlot.png"))
        self.assertGreater(
            os.path.getsize("integrationResultPlot.png"), 0, "Result file is empty"
        )
        os.remove(plot_file)

    def tearDown(self):
        pass
