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
from numpy import array
import sys
import re

def argp():
    '''Parse arguments given by the user in the command line.'''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-cube", dest="cubeFile", required=True, help="density cube file")
    parser.add_argument("-d_range", dest="d_range", required=True, help="density range", type=float, action='append')
    args = parser.parse_args()
    return args

class DensityCube(object):
    '''
    Class for handling Gaussian cube file operations and calculations
    '''
    
    def __init__(self):
        '''Initialize necessary variables and file handlers.'''
        self.densityCubeFile = open(argp().cubeFile)
        self.densityCoordFile = open('density_coord.csv', 'w')
        self.densityIntegralsFile = open('integrals.csv', 'w')
        self.densityDerivativesFile = open('derivatives', 'w')
        self.grid_list = []

    def Proccess_cube_file(self):
        '''Process the cube file to extract density, atom, and voxel information.'''
        self.title = ''
        finishTitle = 0
        lineNumber = 0

        # Extract title
        while finishTitle == 0:
            line = self.densityCubeFile.readline()
            if re.search('[a-zA-Z]+', line):
                self.title += line
                lineNumber += 1
                lastPos = self.densityCubeFile.tell()
            else:
                finishTitle = 1
        self.title = self.title[:-1]
        self.densityCubeFile.seek(lastPos)

        # Extract atom information
        self.atomsLine = self.densityCubeFile.readline()
        self.atomsNumber, self.atom_x_origin, self.atom_y_origin, self.atom_z_origin = self.atomsLine.split()
        self.atomsNumber = int(self.atomsNumber)

        # Extract voxel information
        self.xVoxelLine = self.densityCubeFile.readline().split()
        self.xVoxel, self.xVoxelWidth = float(self.xVoxelLine[0]), float(self.xVoxelLine[1])
        self.yVoxelLine = self.densityCubeFile.readline().split()
        self.yVoxel, self.yVoxelWidth = float(self.yVoxelLine[0]), float(self.yVoxelLine[1])
        self.zVoxelLine = self.densityCubeFile.readline().split()
        self.zVoxel, self.zVoxelWidth = float(self.zVoxelLine[0]), float(self.zVoxelLine[1])

        # Store atoms information in a dictionary
        self.atomsDict = {}
        for i in range(self.atomsNumber):
            lineNumber += 1
            atomLine = self.densityCubeFile.readline()
            atomNr, atomNr2, atom_x, atom_y, atom_z = atomLine.split()
            self.atomsDict[i] = [atomNr, atom_x, atom_y, atom_z]

        # Store grid values
        for line in self.densityCubeFile.readlines():
            for grid in line.split():
                self.grid_list.append(grid)

    def d2f_dz2(self, z, y, x):
        '''Calculate second derivative with respect to z.'''
        if 0 < z < self.xVoxel - 1:
            return (a[z+1][y][x] - 2 * float(a[z][y][x]) + a[z-1][y][x]) / (self.delta_z * self.delta_z)
        return None

    def d2f_dy2(self, z, y, x, delta_y):
        '''Calculate second derivative with respect to y.'''
        if 0 < y < self.zVoxel - 1:
            return (a[z][y+1][x] - 2 * float(a[z][y][x]) + a[z][y-1][x]) / (self.delta_y * self.delta_y)
        return None

    def d2f_dx2(self, z, y, x):
        '''Calculate second derivative with respect to x.'''
        if 0 < x < self.zVoxel - 1:
            return (a[z][y][x+1] - 2 * float(a[z][y][x]) + a[z][y][x-1]) / (self.delta_x * self.delta_x)
        return None

    def d2f_dzdy(self, z, y, x):
        '''Calculate mixed second derivative with respect to z and y.'''
        if 0 < y < self.zVoxel - 1 and 0 < z < self.xVoxel - 1:
            return (a[z+1][y+1][x] - a[z+1][y-1][x] - a[z-1][y+1][x] + a[z-1][y-1][x]) / (4 * self.delta_z * self.delta_y)
        return None

    def d2f_dzdx(self, z, y, x):
        '''Calculate mixed second derivative with respect to z and x.'''
        if 0 < x < self.yVoxel - 1 and 0 < z < self.xVoxel - 1:
            return (a[z+1][y][x+1] - a[z+1][y][x-1] - a[z-1][y][x+1] + a[z-1][y][x+1]) / (4 * self.delta_z * self.delta_y)
        return None

    def d2f_dydx(self, z, y, x):
        '''Calculate mixed second derivative with respect to y and x.'''
        if 0 < y < self.zVoxel - 1 and 0 < x < self.zVoxel - 1:
            return (a[z][y+1][x+1] - a[z][y+1][x-1] - a[z][y-1][x+1] + a[z][y-1][x-1]) / (4 * self.delta_y * self.delta_x)
        return None

    def d2f_dydz(self, z, y, x):
        '''Calculate mixed second derivative with respect to y and z.'''
        if 0 < y < self.zVoxel - 1 and 0 < z < self.xVoxel - 1:
            return (a[z+1][y+1][x] - a[z-1][y+1][x] - a[z+1][y-1][x] + a[z+1][y-1][x]) / (4 * self.delta_y * self.delta_z)
        return None

    def d2f_dxdz(self, z, y, x):
        '''Calculate mixed second derivative with respect to x and z.'''
        if 0 < x < self.yVoxel - 1 and 0 < z < self.xVoxel - 1:
            return (a[z+1][y][x+1] - a[z-1][y][x+1] - a[z+1][y][x-1] + a[z+1][y][x-1]) / (4 * self.delta_x * self.delta_z)
        return None

    def d2f_dxdy(self, z, y, x):
        '''Calculate mixed second derivative with respect to x and y.'''
        if 0 < x < self.yVoxel - 1 and 0 < y < self.zVoxel - 1:
            return (a[z][y+1][x+1] - a[z][y-1][x+1] - a[z][y-1][x-1] + a[z][y-1][x-1]) / (4 * self.delta_x * self.delta_y)
        return None

    def hessian(self, z, y, x):
        '''Calculate the Hessian matrix at a given point.'''
        return (self.d2f_dx2(z, y, x), self.d2f_dxdy(z, y, x), self.d2f_dzdx(z, y, x),
                self.d2f_dxdy(z, y, x), self.d2f_dy2(z, y, x), self.d2f_dzdy(z, y, x),
                self.d2f_dxdz(z, y, x), self.d2f_dydz(z, y, x), self.d2f_dz2(z, y, x))

    def createCoordsDensityFile(self):
        '''Create a file with coordinates and density values.'''
        self.coord_list = []  # List of coordinates for all grid points

        # Calculate coordinates for each grid point
        for i in range(int(self.xVoxel)):
            for k in range(int(self.yVoxel)):
                for j in range(int(self.zVoxel)):
                    self.coord_list.append([self.xCoord + self.xVoxelWidth * int(i),
                                           self.yCoord + self.yVoxelWidth * int(k),
                                           self.zCoord + self.zVoxelWidth * int(j)])

        # Write coordinates and corresponding densities to the file
        for i in range(len(self.coord_list)):
            print >> self.densityCoordFile, self.grid_list[i], self.coord_list[i]

    def IntegrateDensity(self):
        '''Integrate the density distribution over a given range.'''
        self.densityCubeFile.seek(0)
        for i in range(10):
            self.densityCubeFile.readline()  # Jump to 10th line (start of density values)

        self.densityList = self.densityCubeFile.read().split()
        densityArray = array(self.densityList, float).reshape(83, 78, 78)

        integral = 0
        for z in range(len(densityArray)):
            for y in range(len(densityArray[z])):
                for x in range(len(densityArray[z][y])):
                    if self.d_range[0] < densityArray[z][y][x] < self.d_range[1]:
                        integral += pow(self.densityCubeFile, 3)

        print >> self.densityIntegralsFile, integral


if __name__ == "__main__":
    '''Main function to execute the Cube file processing'''
    densityCube = DensityCube()
    densityCube.Proccess_cube_file()
    densityCube.createCoordsDensityFile()
    densityCube.IntegrateDensity()

