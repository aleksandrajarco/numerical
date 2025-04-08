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

        """for line in lines:
                self.grid_data.extend([float(val) for val in line.split()]) """


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
        print(f"xVoxel Info: {self.xVoxel}")
        print(f"self.xVoxelWidth: {self.xVoxelWidth}")
        print(f"yVoxel Info: {self.yVoxel}")
        print(f"zVoxel Info: {self.zVoxel}")

    def d2f_dz2(self, z, y, x):  # 1
        if z < self.xVoxel - 1 and z > 0:
            return (((a[z + 1][y][x]) - 2 * float(a[z][y][x]) + (a[z - 1][y][x])) / (self.delta_z * self.delta_z))
        else:
            return None

    def d2f_dy2(self, z, y, x, delta_y):  # 2
        if y < self.zVoxel - 1 and y > 0:
            return (((a[z][y + 1][x]) - 2 * float(a[z][y][x]) + (a[z][y - 1][x])) / (self.delta_y * self.delta_y))
        else:
            return None

    def d2f_dx2(self, z, y, x):  # 3
        if x < self.zVoxel - 1 and x > 0:
            return ((a[z][y][x + 1] - 2 * float(a[z][y][x]) + (a[z][y][x - 1])) / (self.delta_x * self.delta_x))
        else:
            return None

    def d2f_dzdy(self, z, y, x):  # 4
        if y < self.zVoxel - 1 and y > 0 and z < self.xVoxel - 1 and z > 0:
            return ((a[z + 1][y + 1][x] - float(a[z + 1][y - 1][x]) - float(a[z - 1][y + 1][x]) + float(
                a[z - 1][y - 1][x])) / (4 * self.delta_z * self.delta_y))
        else:
            return None

    def d2f_dzdx(self, z, y, x):  # 5
        if x < self.yVoxel - 1 and x > 0 and z < self.xVoxel - 1 and z > 0:
            return (((a[z + 1][y][x + 1] - float(a[z + 1][y][x - 1]) - float(a[z - 1][y][x + 1]) + float(
                a[z - 1][y][x + 1])) / (4 * self.delta_z * self.delta_y)))
        else:
            return None

    def d2f_dydx(self, z, y, x):  # 6
        if y < self.zVoxel - 1 and y > 0 and x < self.zVoxel - 1 and x > 0:
            return (((a[z][y + 1][x + 1] - float(a[z][y + 1][x - 1]) - float(a[z][y - 1][x + 1]) + float(
                a[z][y - 1][x + 1])) / (4 * self.delta_y * self.delta_x)))
        else:
            return None

    def d2f_dydz(self, z, y, x):  # 7
        if y < self.zVoxel - 1 and y > 0 and z < self.xVoxel - 1 and z > 0:
            return (((a[z + 1][y + 1][x] - float(a[z - 1][y + 1][x]) - float(a[z + 1][y - 1][x]) + float(
                a[z + 1][y - 1][x])) / (4 * self.delta_y * self.delta_z)))
        else:
            return None

    def d2f_dxdz(self, z, y, x):  # 8
        if x < self.yVxel - 1 and x > 0 and z < self.xVoxel - 1 and z > 0:
            return (((a[z + 1][y][x + 1] - float(a[z - 1][y][x + 1]) - float(a[z + 1][y][x - 1]) + float(
                a[z + 1][y][x - 1])) / (4 * self.delta_x * self.delta_z)))
        else:
            return None

    def d2f_dxdy(self, z, y, x):  # 9
        if x < self.yVoxel - 1 and x > 0 and y < self.zVoxel - 1 and y > 0:
            return ((a[z][y + 1][x + 1] - float(a[z][y - 1][x + 1]) - float(a[z][y - 1][x - 1]) + float(
                a[z][y - 1][x - 1])) / (4 * self.delta_x * self.delta_y))
        else:
            return None
        
    def hessian(self, z, y, x):
        '''Calculate the Hessian matrix at a given point.'''
        return (self.d2f_dx2(z, y, x), self.d2f_dxdy(z, y, x), self.d2f_dzdx(z, y, x),
                self.d2f_dxdy(z, y, x), self.d2f_dy2(z, y, x), self.d2f_dzdy(z, y, x),
                self.d2f_dxdz(z, y, x), self.d2f_dydz(z, y, x), self.d2f_dz2(z, y, x))

    def write_coords_and_density(self, output_file='density_coords.csv'):

        with open('output_file.csv', 'w') as f:
            f.write("density,x,y,z\n")  # Write headers to the file

            # Iterate over the grid dimensions and calculate the coordinates
            for z in range(int(self.zVoxel)):
                for y in range(int(self.yVoxel)):
                    for x in range(int(self.xVoxel)):
                        # Calculate the 3D coordinates based on the voxel information
                        coord_x = self.atom_x_origin + self.xVoxelWidth * x
                        coord_y = self.atom_y_origin + self.yVoxelWidth * y
                        coord_z = self.atom_z_origin + self.zVoxelWidth * z

                        # Extract the density value at (z, y, x)
                        density = self.grid_data[x,y,z]

                        # Write the density and the coordinates (x, y, z) to the file
                        f.write(f"{density}, {coord_x}, {coord_y}, {coord_z}\n")

    def integrate_density(self):
        '''Integrate the electron density over the cube grid.'''
        total_density = 0.0

        # Calculate the volume of a single voxel (assuming orthogonal axes)
        voxel_volume = self.xVoxelWidth * self.yVoxelWidth * self.zVoxelWidth
        print(f"voxel_width:{self.yVoxelWidth}")
        for x in range(self.xVoxel):
            for y in range(self.yVoxel):
                for z in range(self.zVoxel):
                    density = self.grid_data[x, y, z]

                    # Check if the density value is within the specified range
                    if self.d_range[0] <= density <= self.d_range[1]:
                        total_density += density * voxel_volume

        print(f"Integrated electron density (number of electrons): {total_density}")


if __name__ == "__main__":
    '''Main function to execute the Cube file processing'''
    args = parse_arguments()
    cube = DensityCube(args.cube_file, args.d_range)
    #print(cube.header)
   # cube.write_coords_and_density()
    cube.integrate_density()


