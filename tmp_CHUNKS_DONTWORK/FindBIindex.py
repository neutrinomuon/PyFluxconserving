#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 12:07:21 2021

@author: jean
"""

import numpy as np
import FindBIindex as interp
import matplotlib.pyplot as plt
import random

class FindBIindex(object):
    
    def __init__(self,x_new,x,y,IsKeepOn=0,verbosity=0):
        self.x = x
        self.y = y
        self.x_new = x_new
        
        #self.y_new,ilastval = interp.findbiindex(self.x_new,self.x,self.y,ilastval,verbosity)
        
        #print(self.y_new)
    
    def TestRoutine():
        print("... Cerate a sine wave")
        
        x   = np.arange(0,2.0*np.pi,0.2)
        y   = np.sin(x)
        ind = np.arange(0,np.size(x),1) 
        random.shuffle(ind)
        ind = np.sort(ind)

        x = x[ind]
        y = y[ind]

        # ny_value = lininterpol(nx_value,xold_vec,yold_vec,ilastval,[verbosity])
        verbosity = 1
        IsKeepOn  = 1
        
        x_new = [x[0]-1.0,x[0]-0.1,x[0],0.5,1.0,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1,x[-1]+0.3]
        
        w_new = interp.findbiindex(x, x_new, IsKeepOn, verbosity)
        y_new = interp.finterpolar(y,w_new)
        
        # Or you can call directly the subroutine
        y_int = interp.indexinterp(x_new,x,y, IsKeepOn, verbosity )
        
        print("... Interpolation of array: {0:}".format(x_new))
        
        #x_new = [x[0]-0.5,x[0]-0.1,x[0],0.5,1.0,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1,x[-1]+0.5]
        #y_new = interp.lininterpol(x_new,x,y,ilastval,verbosity)
        
        #print("... Interpolation of array: {0:}".format(x_new))
        
        #print(ilastval)
        
        plt.scatter(x,y)
        plt.plot(x,y)
        plt.plot(x_new,y_new,color='orange')
        plt.scatter(x_new,y_new,color='red')
        
        plt.plot(x_new,y_int,color='purple',linestyle='-.')
        plt.xlabel("radians")
        plt.ylabel("sin(x)")
        
def main():
    i_object = FindBIindex
    i_object.TestRoutine()
    
    x = np.arange(0,2.0*np.pi,0.2)
    y = np.sin(x**2)
    
    x_new = [0.5,1.0,2.3,3.22222,4.66,5.2,6.0]
    
    i_object = FindBIindex(x_new,x,y)
    #print( i_object.y_new )

if __name__ == "__main__":
    main()
