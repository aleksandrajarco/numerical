'''
Created on Jan 18, 2014

@author: ola
'''

from math import cos
import argparse
from random import uniform
import os
import tempfile
import matplotlib.pyplot as plt
from typing import List


def parse_args() -> argparse.Namespace:
    """
    Parse arguments from the command line.

    Returns:
        argparse.Namespace: Parsed command-line arguments with the following attributes:
            - xrange (List[float]): Boundary of the integration range as [x1, x2].
            - yrange (List[float]): Minimum and maximum of the function in the [xa, xb] range as [ymin, ymax].
            - nsteps (int): Total number of steps for Monte Carlo integration.
            - stepsize (int): Step size for iteration.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-x", dest="xrange", required=True, help="Boundary of integration range", type=float, nargs=2)
    parser.add_argument("-y", dest="yrange", required=True,
                        help="Minimum and maximum of the function in [xa, xb] range", type=float, nargs=2)
    parser.add_argument("-n", dest="nsteps", required=False, help="Number of steps", type=int, default=50000)
    parser.add_argument("-s", dest="stepsize", required=False, help="Step size", type=int, default=1000)

    args = parser.parse_args()

    # Validate xrange and yrange
    x1, x2 = args.xrange
    if x1 >= x2:
        parser.error("xrange must have a smaller value first (x1 < x2).")
    y1, y2 = args.yrange
    if y1 >= y2:
        parser.error("yrange must have a smaller value first (ymin < ymax).")
    if args.stepsize <= 0 or args.nsteps <= 0:
        parser.error("Step size and number of steps must be positive integers.")

    return args


def plot_results(file_name: str, output_image: str = "MCIntegration.png") -> None:
    """
    Generate and save a plot from a results file.

    Args:
        file_name (str): Name of the input file containing Monte Carlo integration results.
        output_image (str): Name of the output image file to save the plot. Default is 'MCIntegration.png'.
    """
    x_vals: List[float] = []
    y_vals: List[float] = []

    with open(file_name, 'r') as file:
        for line in file:
            n, integral = map(float, line.strip().split())
            x_vals.append(n)
            y_vals.append(integral)

    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_vals, label="Monte Carlo Integration", color='blue', linestyle='-', marker='o')
    plt.xlabel("Number of Steps (n)")
    plt.ylabel("Integral Value")
    plt.title("Monte Carlo Integration")
    plt.legend()
    plt.grid(True)
    plt.savefig(output_image)
    plt.close()


class MonteCarlo:
    """
    Class to perform and visualize Monte Carlo integration.

    Attributes:
        args (argparse.Namespace): The parsed command-line arguments.
        data_file (tempfile.NamedTemporaryFile): Temporary file to store Monte Carlo integration results.
    """

    def __init__(self, args=None):
        """
        Initialize the MonteCarlo class by parsing command-line arguments
        and creating a temporary file to store results.
        """
        self.args = args if args is not None else parse_args()
        self.data_file: tempfile.NamedTemporaryFile = tempfile.NamedTemporaryFile(delete=False, mode='w')

    def monte_carlo_integration(self) -> str:
        """
        Perform Monte Carlo integration based on user input and save results to a file.

        Args:
            args (argparse.Namespace): The parsed command-line arguments specifying the integration parameters.

        Returns:
            str: The name of the file where the Monte Carlo integration results are stored.
        """
        nsteps: int = self.args.nsteps
        stepsize: int = self.args.stepsize
        x1, x2 = self.args.xrange
        ymin, ymax = self.args.yrange

        with self.data_file:
            for n in range(stepsize, nsteps, stepsize):  # Iterate through steps
                p: float = 0.0
                for _ in range(n):
                    x = uniform(x1, x2)  # Random x within [x1, x2]
                    y = uniform(ymin, ymax)  # Random y within [ymin, ymax]
                    if abs(cos(x)) >= abs(y):  # Check if point under the curve
                        p += 1
                r = p / n
                integral: float = r * (x2 - x1) * (ymax - ymin)  # Calculate integral
                self.data_file.write(f"{n} {integral}\n")  # Write to temporary file
            self.data_file.flush()  # Ensure the data is written to the file
        return self.data_file.name


if __name__ == '__main__':
    args = parse_args()  # Parse command-line arguments
    MC = MonteCarlo()  # Create an instance of MonteCarlo
    MC.monte_carlo_integration()  # Perform Monte Carlo integration
    plot_results(MC.data_file.name, 'integrationResultPlot.png')  # Generate the plot
