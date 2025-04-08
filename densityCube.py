'''
Created on Jan 25, 2014

Author: ola

Script reads Gaussian cube file (cubefile.txt) and performs the following tasks:
1) Reads the coordinates of a molecule and the electron density distribution values.
2) Evaluates first and second derivatives (partial) at each point using numerical differentiation.
3) Integrates the density distribution (via numerical integration) and calculates the number of electrons. Compares the result with the exact number of electrons and calculates the error.
4) Evaluates eigenvalues of Hessian matrices for each point in the cube file.
5) Finds the approximate positions of atoms by analyzing the electronic density distributions (stationary points, local maxima of electron density).
'''

import argparse

import numpy as np
from numpy import array
import sys
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process Gaussian cube file")
    parser.add_argument("-cube", dest="cube_file", required=True, help="Path to the cube file")
    parser.add_argument("-d_range", dest="d_range", required=True, nargs=2, type=float, help="Density range")
    return parser.parse_args()

class DensityCube(object):
    '''
    Class for handling Gaussian cube file operations and calculations
    '''

    def __init__(self, cube_file_path, d_range):
        self.cube_file_path = cube_file_path
        self.d_range = d_range
        self.grid_data = []
        self.atom_info = []
        self.grid_shape = []
        self.header = []
        self.title = ""
        self.xVoxel = None
        self.yVoxel = None
        self.zVoxel = None
        self._read_cube_file()
        self._parse_header()

    def _read_cube_file(self):
        with open(self.cube_file_path, 'r') as f:
            self.header = [next(f) for _ in range(6)]  # title + atom + voxel info
            atom_count = int(self.header[2].split()[0])
            self.atoms = [next(f) for _ in range(atom_count)]
            lines = f.readlines()
            self.grid_data = []

            # Parse the grid data (density values) and store it as a list of floats
            for line in lines:
                for grid in line.split():
                    self.grid_data.append(float(grid))

        self.grid_data = np.array(self.grid_data)  # Convert to NumPy array

    def _parse_header(self):
        '''Parses the header section of the cube file.'''
        # Parse the title lines
        self.title = self.header[:2]

        # Parse atom info (3rd line)
        atom_info = self.header[2].split()
        self.atom_count = int(atom_info[0])  # Number of atoms
        self.atom_x_origin = float(atom_info[1])
        self.atom_y_origin = float(atom_info[2])
        self.atom_z_origin = float(atom_info[3])

        # Parse voxel info (next 3 lines)
        self.xVoxelLine = self.header[3].split()
        self.xVoxel, self.xVoxelWidth = int(self.xVoxelLine[0]), float(self.xVoxelLine[1])

        self.yVoxelLine = self.header[4].split()
        self.yVoxel, self.yVoxelWidth = int(self.yVoxelLine[0]), float(self.yVoxelLine[2])

        self.zVoxelLine = self.header[5].split()  # Correct the index for zVoxel
        self.zVoxel, self.zVoxelWidth = int(self.zVoxelLine[0]), float(self.zVoxelLine[3])

        # create grid_shape array
        self.grid_shape = np.array([self.xVoxel, self.yVoxel, self.zVoxel])

        # reshape
        self.grid_data = np.array(self.grid_data).reshape(self.grid_shape)

        # Output parsed information for verification (optional)
        print(f"Title: {self.title}")
        print(f"Atom Info: {self.atom_count}, ({self.atom_x_origin}, {self.atom_y_origin}, {self.atom_z_origin})")
        print(f"yVoxel Info: {self.yVoxel}")
        print(f"zVoxel Info: {self.zVoxel}")

    # First derivative with respect to x
    def dfdx(self, z, y, x):
        if x < self.zVoxel - 1:
            return (self.grid_data[z, y, x + 1] - self.grid_data[z, y, x]) / self.xVoxelWidth
        else:
            return None

    # First derivative with respect to y
    def dfdy(self, z, y, x):
        if y < self.yVoxel - 1:
            return (self.grid_data[z, y + 1, x] - self.grid_data[z, y, x]) / self.yVoxelWidth
        else:
            return None

    # First derivative with respect to z
    def dfdz(self, z, y, x):
        if z < self.xVoxel - 1:
            return (self.grid_data[z + 1, y, x] - self.grid_data[z, y, x]) / self.zVoxelWidth
        else:
            return None

    # Second derivative with respect to x (d^2f/dx^2)
    def d2f_dx2(self, z, y, x):
        if x < self.zVoxel - 1 and x > 0:
            return ((self.grid_data[z, y, x + 1] - 2 * self.grid_data[z, y, x] + self.grid_data[
                        z, y, x - 1]) / (self.xVoxelWidth ** 2))
        else:
            return None

    # Second derivative with respect to y (d^2f/dy^2)
    def d2f_dy2(self, z, y, x):
        if y < self.yVoxel - 1 and y > 0:
            return ((self.grid_data[z, y + 1, x] - 2 * self.grid_data[z, y, x] + self.grid_data[
                        z, y - 1, x]) / (self.yVoxelWidth ** 2))
        else:
            return None

    # Second derivative with respect to z (d^2f/dz^2)
    def d2f_dz2(self, z, y, x):
        if z < self.xVoxel - 1 and z > 0:
            return ((self.grid_data[z + 1, y, x] - 2 * self.grid_data[z, y, x] + self.grid_data[
                        z - 1, y, x]) / (self.zVoxelWidth ** 2))
        else:
            return None

    #Mixed second derivative with respect to x and y (d ^ 2f / dxdy)
    def d2f_dxdy(self, z, y, x):
        if x < self.zVoxel - 1 and y < self.yVoxel - 1:
            return (self.grid_data[z, y + 1, x + 1] - self.grid_data[z, y - 1, x + 1] - self.grid_data[
                z, y + 1, x - 1] + self.grid_data[z, y - 1, x - 1]) / (2 * self.xVoxelWidth * self.yVoxelWidth)
        else:
            return None

    # Mixed second derivative with respect to x and z (d^2f/dxdz)
    def d2f_dxdz(self, z, y, x):
        if x < self.zVoxel - 1 and z < self.xVoxel - 1:
            return (self.grid_data[z + 1, y, x + 1] - self.grid_data[z - 1, y, x + 1] - self.grid_data[
                z + 1, y, x - 1] + self.grid_data[z - 1, y, x - 1]) / (2 * self.xVoxelWidth * self.zVoxelWidth)
        else:
            return None

    # Mixed second derivative with respect to y and z (d^2f/dydz)
    def d2f_dydz(self, z, y, x):
        if y < self.yVoxel - 1 and z < self.xVoxel - 1:
            return (self.grid_data[z + 1, y + 1, x] - self.grid_data[z - 1, y + 1, x] - self.grid_data[
                z + 1, y - 1, x] + self.grid_data[z - 1, y - 1, x]) / (2 * self.yVoxelWidth * self.zVoxelWidth)
        else:
            return None
    def calculate_hessian(self, z, y, x):
        '''Calculate the Hessian matrix at a given point.'''
        # Second-order derivatives in each direction
        d2f_dx2 = self.d2f_dx2(z, y, x)
        d2f_dy2 = self.d2f_dy2(z, y, x)
        d2f_dz2 = self.d2f_dz2(z, y, x)

        # Mixed partial derivatives
        d2f_dxdy = self.d2f_dxdy(z, y, x)
        d2f_dxdz = self.d2f_dxdz(z, y, x)
        d2f_dydz = self.d2f_dydz(z, y, x)

        # Constructing the Hessian matrix
        hessian_matrix = np.array([[d2f_dx2, d2f_dxdy, d2f_dxdz],
                                  [d2f_dxdy, d2f_dy2, d2f_dydz],
                                  [d2f_dxdz, d2f_dydz, d2f_dz2]])

        return hessian_matrix

    def calculate_hessian_for_all_points(self):
        """Calculate Hessian for all points in the grid using numpy."""
        # Create grid indices
        x, y, z = np.indices(self.grid_shape)  # Generate 3D grid of indices

        # Flatten the indices and grid data for vectorized processing
        x, y, z = x.flatten(), y.flatten(), z.flatten()
        grid_data = self.grid_data.flatten()

        # Compute Hessian matrices for all points in parallel
        hessians = np.array([self.calculate_hessian(x_i, y_i, z_i) for x_i, y_i, z_i  in zip(x, y, z)])

        return hessians

    def write_coords_and_density(self, output_file='density_coords.csv'):

        with open(output_file, 'w') as f:
            f.write("density,x,y,z\n")  # Write headers to the file

            # Generate a grid of coordinates using meshgrid
            x_coords = self.atom_x_origin + self.xVoxelWidth * np.arange(self.xVoxel)
            y_coords = self.atom_y_origin + self.yVoxelWidth * np.arange(self.yVoxel)
            z_coords = self.atom_z_origin + self.zVoxelWidth * np.arange(self.zVoxel)

            # Use np.meshgrid to create a grid of coordinates for all x, y, z
            X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords)

            # Flatten the grid to iterate over all points
            for x, y, z, density in zip(X.flatten(), Y.flatten(), Z.flatten(), self.grid_data.flatten()):
                f.write(f"{density},{x},{y},{z}\n")

    def integrate_density(self):
        '''Integrate the electron density over the cube grid.'''
        total_density = 0.0

        # Calculate the volume of a single voxel (assuming orthogonal axes)
        voxel_volume = self.xVoxelWidth * self.yVoxelWidth * self.zVoxelWidth
        print(f"voxel_width:{self.yVoxelWidth}")

        mask = (self.grid_data >= self.d_range[0]) & (self.grid_data <= self.d_range[1])
        total_density = np.sum(self.grid_data[mask]) * voxel_volume

        print(f"Integrated electron density (number of electrons): {total_density}")


if __name__ == "__main__":
    '''Main function to execute the Cube file processing'''
    args = parse_arguments()
    cube = DensityCube(args.cube_file, args.d_range)
    #print(cube.header)
    #cube.write_coords_and_density()
    #cube.integrate_density()
    cube.calculate_hessian_for_all_points()


