'''
Created on Jan 25, 2014

@author: ola

Script reads Gaussian cube file (cubefile.txt) and does these things: 
1) reads the coordinates of a molecule and the electron density distribution values
2) evaluates in each point first and second derivatives (partial) by utilizing numerical differentiation 
3) integrates the density distribution (by using numerical integration) and finds the number of electrons by this. Compare with the exact number of electons and error 
4) evaluate eigenvalues of Hessian matrices for each point in cube file 
5) Find the approximate positions of atoms by analyzing electronic density distributions (stationary points, local maximas of electron density - find them)

'''
import argparse
from numpy import array
import sys


def argp():
    ''' parse arguments given by user in command line'''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-cube", dest="cubeFile", required=True, help="density cube file")
    parser.add_argument("-d_range", dest="d_range", required=True, help="density range", type=float, action='append')
    args = parser.parse_args()
    return args

class DensityCube(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        #self.densityCubeFile=open(sys.argv[1])
        #print argp()
        #sys.exit()
        self.densityCubeFile=open(argp().cubeFile)
        self.densityCoordFile=open('density_coord.csv','w')
        self.densityIntegralsFile=open('integrals.csv','w')
        
    def readDensityCubeFile(self):
        
        self.d_range=argp().d_range
        #print self.d_range
        #sys.exit()
        coordsLine= self.densityCubeFile.readlines()[2].split()
        self.xCoord=float(coordsLine[1])
        self.yCoord=float(coordsLine[2])
        self.zCoord=float(coordsLine[3])
        #print self.xCoord, self.yCoord,self.zCoord
        self.densityCubeFile.seek(0)

        self.xVoxelLine=self.densityCubeFile.readlines()[3].split()
        self.xVoxel, self.xVoxelWidth=float(self.xVoxelLine[0]), float(self.xVoxelLine[1])
        self.densityCubeFile.seek(0)
        self.yVoxelLine=self.densityCubeFile.readlines()[4].split()
        self.yVoxel, self.yVoxelWidth=float(self.yVoxelLine[0]), float(self.yVoxelLine[2])
        self.densityCubeFile.seek(0)
        self.zVoxelLine=self.densityCubeFile.readlines()[5].split()
        self.zVoxel, self.zVoxelWidth=float(self.zVoxelLine[0]), float(self.zVoxelLine[3])
        
        self.densityCubeFile.seek(0)
        self.grid_list=[] 
        
        for line in self.densityCubeFile.readlines()[10:]:
            for grid in line.split():
                self.grid_list.append(grid)
        
        #self.grid_list=[grid for grid in  [line.split() for line  in self.densityCubeFile.readlines()[10:]]]
        #print self.grid_list[0],self.grid_list[-1]
        print self.xVoxel, self.yVoxel, self.zVoxel
        print self.xVoxelWidth, self.yVoxelWidth, self.zVoxelWidth
        
        
    def createCoordsDensityFile(self):
        
        self.coord_list=[] #lista wpolrzednych wszystkich gridow
        
        for i in range(int(self.xVoxel)):
            for k in range(int(self.yVoxel)):
                for j in range(int(self.zVoxel)):
                    #print self.xCoord,self.xVoxelWidth,int(i)
                    self.coord_list.append([self.xCoord+self.xVoxelWidth*int(i),self.yCoord+self.yVoxelWidth*int(k),self.zCoord+self.zVoxelWidth*int(j)])
            
        
        print len(self.coord_list), len(self.grid_list)
        for i in range(len(self.coord_list)):
            print >> self.densityCoordFile,self.grid_list[i],self.coord_list[i]
            

    
    def IntegrateDensity(self):
        #plik_inp=open('density.cube2')
        
        self.densityCubeFile.seek(0)
        for i in range(10):
            self.densityCubeFile.readline() #jump to 10th line (start of density values)
        
        #print self.densityCubeFile.readline()
        #sys.exit()
        self.densityList=self.densityCubeFile.read().split()
        
        #print a
        densityArray=array(self.densityList,float).reshape(83,78,78)
        
        
        #integrating inside some isodensity value
        #d_range=(0,1)# add me 
        print self.xVoxelWidth
        
        
        integral=0
        for z in range(0,len(densityArray)):
                for y in range(0,len(densityArray[z])):
                        for x in range(0,len(densityArray[z][y])):
                            if densityArray[z][y][x]>self.d_range[0] and densityArray[z][y][x]<self.d_range[1]:
                                integral+= pow(self.xVoxelWidth,3)/8*( #all voxels have the same width
                                                densityArray[z-1][y-1][x-1]+
                                                densityArray[z][y-1][x-1]+
                                                densityArray[z-1][y][x-1]+
                                                densityArray[z][y][x-1]+
                                                densityArray[z-1][y-1][x]+
                                                densityArray[z][y-1][x]+
                                                densityArray[z][y][x]+
                                                densityArray[z-1][y][x])
                                print >> self.densityIntegralsFile, integral
        
        


if __name__ == '__main__':
    DC=DensityCube()
    DC.readDensityCubeFile()
    #DC.createCoordsDensityFile()
    DC.IntegrateDensity()
                