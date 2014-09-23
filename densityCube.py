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
import re

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
        self.densityDerivativesFile=open('derivatives','w')
        self.grid_list=[]
    	

    def Proccess_cube_file(self):
    	self.title=''
    	finishTitle=0
    	lineNumber=0
    	
    	#find title
    	while finishTitle==0:
    	 	line=self.densityCubeFile.readline()
    	 	if re.search('[a-zA-Z]+',line):
    	 		self.title+=line
    	 		lineNumber+=1
    	 		lastPos=self.densityCubeFile.tell()
    	 	else:
    	 		finishTitle=1
        self.title=self.title[0:-1]
    	self.densityCubeFile.seek(lastPos)
    	
    	#find atom
    	self.atomsLine=self.densityCubeFile.readline()
    	self.atomsNumber,self.atom_x_origin,self.atom_y_origin,self.atom_z_origin=self.atomsLine.split()
    	self.atomsNumber=int(self.atomsNumber)
    	#find voxels
    	self.xVoxelLine=self.densityCubeFile.readline().split()
    	self.xVoxel, self.xVoxelWidth=float(self.xVoxelLine[0]), float(self.xVoxelLine[1])
    	self.yVoxelLine=self.densityCubeFile.readline().split()
    	self.xVoxel, self.yVoxelWidth=float(self.yVoxelLine[0]), float(self.yVoxelLine[1])
    	self.zVoxelLine=self.densityCubeFile.readline().split()    	
    	self.zVoxel, self.zVoxelWidth=float(self.zVoxelLine[0]), float(self.zVoxelLine[1])
    	
    	
    	self.atomsDict={}
    	for i in range(self.atomsNumber):
    		lineNumber+=1
    		atomLine=self.densityCubeFile.readline()
    		atomNr,atomNr2,atom_x,atom_y,atom_z=atomLine.split()
    		self.atomsDict[i]=[atomNr,atom_x,atom_y,atom_z]
    	
    	 		
    	#print self.densityCubeFile.readline()
    	for line in self.densityCubeFile.readlines():
             for grid in line.split():
                 self.grid_list.append(grid)
     #=======================================================================
    	#print self.grid_list[-1]
    	#print atomsDict
    		
    	#self.atomsInfo=
    	#print atomsNumber
     	#print self.atomsLine
    
    def d2f_dz2(self,z,y,x): #1
        if z<self.xVoxel-1 and z>0:
            return (((a[z+1][y][x])-2*float(a[z][y][x])+(a[z-1][y][x]))/(self.delta_z*self.delta_z))
        else:
            return None 
 
    def d2f_dy2(self,z,y,x,delta_y):#2
        if y<self.zVoxel-1 and y>0:        
            return  (((a[z][y+1][x])-2*float(a[z][y][x])+(a[z][y-1][x]))/(self.delta_y*self.delta_y))
        else:
            return None
 
    def d2f_dx2(self,z,y,x): #3
        if x<self.zVoxel-1 and x>0:
            return ((a[z][y][x+1]-2*float(a[z][y][x])+(a[z][y][x-1]))/(self.delta_x*self.delta_x))
        else:
            return None
     
     
    def d2f_dzdy(self,z,y,x): #4
        if y<self.zVoxel-1 and y>0 and z<self.xVoxel-1 and z>0:
            return ((a[z+1][y+1][x]-float(a[z+1][y-1][x])-float(a[z-1][y+1][x])+float(a[z-1][y-1][x]))/(4*self.delta_z*self.delta_y))
        else:
            return None

    
    def d2f_dzdx(self,z,y,x): #5   
        if x<self.yVoxel-1 and x>0 and z<self.xVoxel-1 and z>0:
            return (((a[z+1][y][x+1]-float(a[z+1][y][x-1])-float(a[z-1][y][x+1])+float(a[z-1][y][x+1]))/(4*self.delta_z*self.delta_y)))
        else:
            return None

  
    def d2f_dydx(self,z,y,x): #6    
        if y<self.zVoxel-1 and y>0 and x<self.zVoxel-1 and x>0:
            return (((a[z][y+1][x+1]-float(a[z][y+1][x-1])-float(a[z][y-1][x+1])+float(a[z][y-1][x+1]))/(4*self.delta_y*self.delta_x)))
        else:
            return None
    
    
    def d2f_dydz(self,z,y,x): #7
        if y<self.zVoxel-1 and y>0 and z<self.xVoxel-1 and z>0:
            return  (((a[z+1][y+1][x]-float(a[z-1][y+1][x])-float(a[z+1][y-1][x])+float(a[z+1][y-1][x]))/(4*self.delta_y*self.delta_z)))
        else:
            return None
    
    
    def d2f_dxdz(self,z,y,x): #8
        if x<self.yVxel-1 and x>0 and z<self.xVoxel-1 and z>0:
            return (((a[z+1][y][x+1]-float(a[z-1][y][x+1])-float(a[z+1][y][x-1])+float(a[z+1][y][x-1]))/(4*self.delta_x*self.delta_z)))
        else:
            return None
    
    def d2f_dxdy(self,z,y,x): #9
        if x<self.yVoxel-1 and x>0 and y<self.zVoxel-1 and y>0:
            return ((a[z][y+1][x+1]-float(a[z][y-1][x+1])-float(a[z][y-1][x-1])+float(a[z][y-1][x-1]))/(4*self.delta_x*self.delta_y))
        else:
            return None    
        
        
    def hessian(self,z,y,x):
       	return self.d2f_dx2(z,y,x),self.d2f_dxdy(z,y,x),self.d2f_dzdx(z,y,x),self.d2f_dxdy(z,y,x),self.d2f_dy2(z,y,x),self.d2f_dzdy(z,y,x),self.d2f_dxdz(z,y,x),self.d2f_dydz(z,y,x),self.d2f_dz2(z,y,x)
    
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
        
        
    def Calculate_hessian_for_each_point(self):
    	for z in range(0,len(a)):
        	for y in range(0,len(a[z])):
          		for x in range(0,len(a[z][y])):
            		 print>>self.densityDerivativesFile,self.hessian(z,y,x)
        
        

if __name__ == '__main__':
	
	DC=DensityCube()
	DC.Proccess_cube_file()
	print DC.title
	sys.exit()
   	DC.readDensityCubeFile()
    #sys.exit()
	DC.createCoordsDensityFile()
	DC.Calculate_hessian_for_each_point()
	DC.IntegrateDensity()
                