#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 12:07:21 2021

@author: jean
"""

import numpy as np
#import LINinterpol as interp
import pchipmodule as interp
import matplotlib.pyplot as plt

class pchipmodule(object):
    
    def __init__(self,x_new,x,y,ilastval=-999,verbosity=0):
        self.x = x
        self.y = y
        self.x_new = x_new
        
        #self.y_new,ilastval = interp.pchipmodule.lininterpol(self.x_new,self.x,self.y,ilastval,verbosity)
        incfd = 1
        n     = self.x.size
        #d , ierr = interp.pchipmodule.dpchim(self.x,self.y)
        #print(self.y_new)
    
    def TestRoutine():
        print("... Cerate a sine wave")
        x = np.arange(0,2.0*np.pi,0.2)
        y = np.sin(x)

        # ny_value = lininterpol(nx_value,xold_vec,yold_vec,ilastval,[verbosity])
        ilastval = -999
        verbosity = 0
        
        x_new = [x[0]-0.1,x[0],0.5,1.0,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1]
        x_new = np.array( x_new, dtype=float )
        #y_new,ilastval = interp.lininterpol(x_new,x,y,ilastval,verbosity)

        incfd = 1
        n = x.size

        y_ = np.zeros([1,y.size], dtype=float)
        y_[0,:] = y[:]
        d, ierr = interp.pchipmodule.dpchim( x, y_ )

        next  = 0
        y_new = np.zeros([x_new.size], dtype=float)

        for i in enumerate(y_new):
            #findloc_min = minloc( ARRAY=auxvecxx-auxvecyy(indexing), DIM=1, MASK=auxvecxx>=auxvecyy(indexing) )
            #findloc_max = maxloc( ARRAY=auxvecxx-auxvecyy(indexing), DIM=1, MASK=auxvecxx<=auxvecyy(indexing) )
            #print()

            dif = x-x_new[i[0]]
            findloc_min = np.array( np.where( dif >= 0.0 ) , dtype=int ).flatten()
            findloc_max = np.array( np.where( dif >= 0.0 ) , dtype=int ).flatten()
            #print(findloc_min.size)
            #try:
            #    findloc_min = np.where( dif >= 0.0 ) #.argmin()
            #    findloc_max = np.where( dif <= 0.0 ) #.argmax()

            if findloc_min.size > 0 and findloc_max.size > 0:
                print(findloc_min[0],findloc_max[0])

            #except:
            #    y_new[i[0]] = 0.0

            #fe, next, ierr = dchfev( x1, x2, f1, f2, d1, d2, xe )

        print("... Interpolation of array: {0:}".format(x_new))
        
        x_new = [x[0]-0.5,x[0]-0.1,x[0],0.5,1.0,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1,x[-1]+0.5]
        #y_new,ilastval = interp.lininterpol(x_new,x,y,ilastval,verbosity)
        
        print("... Interpolation of array: {0:}".format(x_new))
        
        #print(ilastval)
        
        plt.plot(x,y)
        #plt.plot(x_new,y_new,color='orange')
        #plt.scatter(x_new,y_new,color='red')
        plt.xlabel("radians")
        plt.ylabel("sin(x)")
        
def main():
    i_object = pchipmodule
    i_object.TestRoutine()
    
    x = np.arange(0,2.0*np.pi,0.2)
    y = np.sin(x**2)
    
    x_new = [0.5,1.0,2.3,3.22222,4.66,5.2,6.0]
    
    #i_object = pchipmodule(x_new,x,y)
    #print( i_object.y_new )

if __name__ == "__main__":
    main()
