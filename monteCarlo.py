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


def parse_args():
    """Parse arguments from command line."""
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-x", dest="xrange", required=True, help="boundary of integration range", type=float, nargs=2)
    parser.add_argument("-y", dest="yrange", required=True, help="minimum and maximum of function in [xa, xb] range", type=float, nargs=2)
    parser.add_argument("-n", dest="nsteps", required=False, help="number of steps", type=int, default=50000)
    parser.add_argument("-s", dest="stepsize", required=False, help="step size", type=int, default=1000)

    return parser.parse_args()


class MonteCarlo():
    def __init__(self):
        self.dataFile=tempfile.NamedTemporaryFile(delete=False, mode='w')
        self.PlotInput=tempfile.NamedTemporaryFile(delete=False, mode='w')

    def create_plot(self):
        """Generate plot using matplotlib."""
        # Read the data from the temporary data file
        x_vals = []
        y_vals = []
        with open(self.dataFile.name, 'r') as file:
            for line in file:
                n, integral = map(float, line.strip().split())
                x_vals.append(n)
                y_vals.append(integral)

        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(x_vals, y_vals, label="Monte Carlo Integration", color='blue', linestyle='-', marker='o')
        plt.xlabel("Number of Steps (n)")
        plt.ylabel("Integral Value")
        plt.title("Monte Carlo Integration")
        plt.legend()
        plt.grid(True)

        # Save the plot to a file
        plt.savefig("MCIntegration.png")  # Save as PNG
        plt.close()  # Close the plot

    def monte_carlo_integration(self, args):
        """Perform Monte Carlo integration based on user input."""
        nsteps = args.nsteps
        stepsize = args.stepsize
        x1, x2 = args.xrange
        ymin, ymax = args.yrange

        for n in range(stepsize, nsteps, stepsize):
            p = 0.0
            for i in range(n):
                x = uniform(x1, x2)
                y = uniform(ymin, ymax)
                if abs(cos(x)) >= abs(y):
                    p += 1
            r = p / n
            integral = r * (x2 - x1) * (ymax - ymin)
            self.dataFile.write(f"{n} {integral}\n")
            self.dataFile.flush()  # Ensure the data is written to the file
        self.dataFile.close()

    def apply(self, args):
        """Run Monte Carlo integration and plot results."""
        self.monte_carlo_integration(args)
        self.create_plot()
        os.unlink(self.PlotInput.name)
        os.unlink(self.dataFile.name)

if __name__ == '__main__':
    args = parse_args()
    MC = MonteCarlo()
    MC.apply(args)
    
