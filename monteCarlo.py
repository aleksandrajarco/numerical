'''
Created on Jan 18, 2014

@author: ola
'''


from math import cos
import argparse
from random import uniform
import os
import tempfile


def argp():
    ''' parse arguments given by user in command line'''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   
    parser.add_argument("-x", dest="xrange", required=True, help="boundary of integration range", type=float, action='append')
    #parser.add_argument("-xb", dest="xb", required=True, help="right boundary of integration range", type=float)
    parser.add_argument("-y", dest="yrange", required=True, help="minimum and maximum of function in [xa,xb] range", type=float, action='append')
    #parser.add_argument("-yb", dest="yb", required=True, help="maximum of function in [xa,xb] range", type=float)
    parser.add_argument("-n", dest="nsteps", required=False, help="number of steps", type=int, default=50000)
    parser.add_argument("-s", dest="steps", required=False, help="step size", type=int, default=1000)

    args = parser.parse_args()
    return args

def function(x):
        
    return cos(x)

class MonteCarlo():
    def __init__(self):
        #self.dataFile=open()
        self.dataFile=tempfile.NamedTemporaryFile(delete=False)
        self.PlotInput=tempfile.NamedTemporaryFile(delete=False)

    def createPlot(self): #argument as a string
        
        
        #OpenInput=PlotInput
        #print PlotInput.name
        #OpenInput=open(PlotInput,'w')
        command=""" set term postscript eps
        set title "Monte Carlo intergration"
        set output "MCIntegration.ps"
        plot "{0}" using 1:2 with lines""".format(self.dataFile.name)
        
        execute ="gnuplot {0}".format(self.PlotInput.name)
        
        self.PlotInput.write(command)
        self.PlotInput.close()
        self.dataFile.close()
        
        #self.PlotInput.close()
        os.system(execute)
        #return execute
    
    
    def monteCarloIntegration(self):
       
        #read arguments from command line
        args=argp()
        
        nsteps=args.nsteps
        steps=args.steps
        x1=args.xrange[0]
        x2=args.xrange[1]
        ymin=args.yrange[0]
        ymax=args.yrange[1]
        
        for n in range(steps,nsteps,steps):
            p=0.0
            for i in range(n):
                x=uniform(x1,x2)
                y=uniform(ymin,ymax)
                #print x,function(x),y
                if abs(function(x))>=abs(y):
                    p+=1
                #print p    
            r=p/n
            #print ymax-ymin,x2-x1,r
            #print r
            integral= r*(x2-x1)*(ymax-ymin)
        
            self.dataFile.write("{0} {1}\n".format(n,integral))
           
    
    def apply(self):
        
        self.monteCarloIntegration()
        self.createPlot()
        os.unlink(self.PlotInput.name)
        os.unlink(self.dataFile.name)
        #os.system(self.createPlot())

    
    
if __name__ == '__main__':
    MC=MonteCarlo()
    MC.apply()
    