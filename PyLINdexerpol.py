"""
Revised interface on Sat Jan 30 12:07:21 2021

@author: Jean Gomes

RESUME :  Interpolation script calling LINdexerpol in Fortran and Scipy!

Version: v01 beta

PYTHON : Python compatibility using f2py revised. Better usage  with numpy.  

Written: Jean Michel Gomes © Copyright
"""

# Import libraries
# Be careful with importation of pylab - It needs to be first otherwise problems with other features
from pylab import *
import matplotlib.pyplot as plt

import numpy as np
import LINdexerpol as interp

import AkimaSpline as interp1

import scipy.interpolate as interpolate
from scipy.interpolate import CubicSpline
from scipy.interpolate import BSpline

class PyLINdexerpol(object):
    """Created on 
Last version on Wed Sep 23 14:33:51 2020

@author: Jean Gomes Copyright (c)

*** Number : 002
1) LINdexerpol
2) author_LINdexerpol

######################################################################
SUBROUTINE LINdexerpol( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,                                        &
                                                 nold_vec,IsKeepOn,verbosity )
lindexerpol in python

RESUME : For a given array of values 'xx_value' for the               
                  abscissa, then return the ordinate array values of           
                  'yy_value' based on a linear interpolation within a          
                  table of pair values [xold_vec, yold_vec]. This              
                  subroutine assumes that the values in the array              
                  xold_vec increases monotonically with yold_vec.                     
                                                                           
Input                    arguments = 5                                         
Output                 arguments = 2                                         
Optional              arguments = 1                                        
Total number of arguments = 98                                        
    
   INPUT  : 01) xx_value   -> Interpolate into x array                            
                  02) nxyvalue   -> # of array elements in vector x and o               
                  03)  xold_vec   -> x vector (abcissas)                                 
                  04) yold_vec    -> y vector (ordenadas)                                
                  05) nold_vec    -> # of array elements in vector t and y               
                  06) verbosity   -> Print & Check screen (OLD fortran way)              
                                                                           
OUTPUT : 01) yy_value   -> Interpolated y array                                
                  02) IsKeepOn -> Flag, if == 0 then there's a problem      
                                                                           
PYTHON : Python compatibility using f2py revised. Better usage        
                  with numpy.                                                  
                                                                           
Written: Jean Michel Gomes © Copyright ®                              
Date: Sun Oct 15 10:00:52 WEST 2006                                
Checked: Tue May  1 16:09:13 WEST 2012                                
                  Fri Dec 28 14:55:10 WET  2012                                
 ######################################################################

######################################################################  
SUBROUTINE author_LINdexerpol( a )
author_lindexerpol in python
######################################################################  

######################################################################             
                 HOW TO USE
                 
                 i_object = LINdexerpol
                 i_object.LINdexerpol(x_new,x,y,[verbosity])
                 
                 where
                 
                 x_new is the new array of points (abscissa) to be flux-conserved
                 x,y are the original arrays (abscissa and ordinate)
                 
                 [verbosity] => verbosity mode [0: No, 1:Yes] 
                 
                 Test ROUTINE
                 
                 i_object = LINdexerpol
                 i_object.TestRoutine()
                 
                 It will produce a plot with the sine interpolation
                 
                 time
                 
                 i_object = LINdexerpol
                 i_object.time()
                 
                 It will test the time it takes to interpolate as a function of elements in the array
######################################################################                  
"""  
     
    def __init__(self,x_new,x,y,verbosity=0):
        self.x = x
        self.y = y
        self.x_new = x_new
        
        self.verbosity = verbosity
        
        self.y_new, self.IsKeepOn = interp.lindexerpol(self.x_new,self.x,verbosity=self.verbosity)
        
        #print(self.y_new)
    
    def TestRoutine():
        print("... Create a sine wave function")
        print("... x and y values for x with equal step size")
        
        x = np.arange(0,2.3*np.pi,1.0)
        y = np.sin(x)
        
        #print(x.size)

        print("... step x: {0:}".format(x[1]-x[0]))

        # No verbosity 
        verbosity = 0
        
        print("... x and y values for x with unequal step size")
        
        #x_new = [x[0]-0.1,x[0],0.5,1.0,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1]
        x_new = [x[0]-0.2,x[0]-0.1,x[0],0.5,1.0,1.5,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1,x[-1]+0.5]
        x_new = np.array( x_new, dtype=float )
        x_new = np.unique(x_new)
        y_new = np.sin(x_new)
        
        dxnew = x_new[1:] - x_new[:-1]
        
        print("... step x_new: {0:}".format(dxnew))
                
        #x_new = [x[0]-0.5,x[0]-0.1,x[0],0.5,1.0,1.5,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1,x[-1]+0.5]
        
        x_interp          = np.arange(x[0]-1.89,x_new[-1]+1.5,0.01)
        #x_interp = np.array( [-0.7,-0.7,-0.7,-0.7,-0.7099999,0.2,0.3,0.3,0.9,1.2,1.2,1.5,2.3,2.3,3.4,3.5], dtype=float )
        x_new_interp = x_interp
    
        print("... Interpolation of array: {0:}".format([x_interp[0],x_interp[-1]]))

        # Calling directly the f2py subroutine
        #b,c,d,IsKeepOn = interp.spline3dcoe(x,y)
        #print(b,c,d)
        #b_new,c_new,d_new,IsKeepOn = interp.spline3dcoe(x_new,y_new)
        #print(b_new,c_new,d_new)
        
        y_interp,IsKeepOn = interp.lindexerpol(x_interp,x,y,verbosity=verbosity)
        y_new_interp,IsKeepOn = interp.lindexerpol(x_new_interp,x_new,y_new,verbosity=verbosity)
        
        # Calling Akima Spline interpolation
        delta_x = 0.1
        y_interp_Akima,IsKeepOn = interp1.akimaspline( x_interp, x, y, delta__x=delta_x, verbosity=verbosity )
        y_new_interp_Akima,IsKeepOn = interp1.akimaspline( x_new_interp, x_new, y_new, delta__x=delta_x, verbosity=verbosity )

        # Cubic Spline using scipy
        cs          = CubicSpline( x,y,bc_type='natural' )
        #spl        = interpolate.splrep(x,y)
        #print(spl[0])
        #print(spl[1])
        #print(spl[2])
        
        k            =  3
        t, c, k         = interpolate.splrep(x, y, s=0, k=k)
        
        # x_new must be striclty increasing
        x_new_ = np.unique(x_new)
        y_new_ = np.sin(x_new_)
        cs_new = CubicSpline( x_new_,y_new_ )
        spl_new = BSpline( x,y,k)
        #print('y_new_interp:', y_new_interp)

        Figure = plt.figure( figsize=(12,10),dpi=120,facecolor='w',edgecolor='w' )
        plt.subplots_adjust(bottom=.02, left=.06, right=.95, top=.98, wspace=0.0, hspace=0.0) 
        
        ax = plt.subplot(111)

        # Top plot ###########################################################
        ax1_top = subplot2grid( (20,20), (0,0), colspan=20, rowspan=10 )                                                 # Sets the position and size of the panel for Plot #01
        ax1_top.axis('on')
        ax1_top.axes.get_xaxis().set_visible(False)

        ylimmin = -3.9
        ylimmax =+4.9 # max( y_interp ) * 1.3
        #ax1_top.set_xlim( xmin,xmax )
        ax1_top.set_ylim( ylimmin,ylimmax )

        # Original values
        ax1_top.plot(x,y,color='darkorange',linewidth=10,label='Original')
        
        # Interpolated function Fortran - f2py
        ax1_top.plot(x_interp, y_interp, color='blue',label='LINdexerpol Interpolation',linewidth=10,alpha=0.3)
        ax1_top.plot(x_interp, y_interp_Akima, color='magenta',label='Akima Spline Interpolation',linestyle='--',linewidth=4,alpha=1.)
        ax1_top.scatter(x_interp, y_interp*0.)

        # Interpolated function using scipy
        #ax1_top.plot(x_interp, cs(x_interp), color='black',linestyle='-',linewidth=1.2,label='Scipy Cubic Spline')
        spline = BSpline(t, c, k, extrapolate=True )
        ax1_top.plot(x_interp, spline(x_interp), color='black',linestyle='-',linewidth=1.2,label='Scipy Cubic Spline')

        ax1_top.legend()
        ax1_top.set_xlabel("Radians")
        ax1_top.set_ylabel("f(x) = sin(x)")
        # Top plot ###########################################################

        # Bottom plot ########################################################
        ax2_top = subplot2grid( (20,20), (10,0), colspan=20, rowspan=9 )                                                # Sets the position and size of the panel for Plot #02
        ax2_top.axis('on')
        
        # Original Different Point Values
        ax2_top.plot(x_new,y_new,color='darkorange',label='Original Diff Points',linewidth=10)
        
        # Interpolated function Fortran - f2py
        ax2_top.plot(x_new_interp, y_new_interp, color='blue',label='LINdexerpol Interpolation',linewidth=10,alpha=0.3)
        ax2_top.plot(x_new_interp, y_new_interp_Akima, color='magenta',label='Akima Spline Interpolation',linestyle='--',linewidth=4,alpha=1.)
        ax2_top.scatter(x_new_interp, y_new_interp*0.)

        # # Interpolated function Python3 interface
        # ax2_top.plot(x_new_interp, y_new_interp_Akima, color='magenta',linestyle='--',linewidth=4,label='Akima Python3 Interface')

        # Interpolated function using scipy
        ax2_top.plot(x_new_interp, cs_new(x_new_interp), color='black',linestyle='-',linewidth=1.2,label='Scipy Cubic Spline')
        #ax2_top.plot(x_new_interp, spl_new(x_new_interp), color='black',linestyle='-',linewidth=1.2,label='Scipy Cubic Spline')

        # ax2_top.legend()
        # ax2_top.set_xlabel("Radians")
        # ax2_top.set_ylabel("f(x) = sin(x)")
        
        plt.show()
        # Bottom plot ########################################################

    def time():
        print("... Create a sine wave function")
        print("... x and y values for x with equal step size")
     
        x = np.arange(0,2.3*np.pi,0.05)
        y = np.sin(x)
     
        print("... step x: {0:}".format(x[1]-x[0]))

        # No verbosity 
        verbosity = 0
     
        delta_x = 0.1

        #i_ = 25

        #sizepoints = np.zeros(i_)
        #timearray  = np.zeros([i_,3])

        np.random.seed(10)
        
        #x_interp_ = np.random.normal(0.0,1.0,2**i_)
        #x_interp_ = np.unique(x_interp_)
        
        val1 = np.arange(1,10,1)
        val2 = np.arange(10,100,2)
        val3 = np.arange(100,1000,10)
        val4 = np.arange(1000,10000,100)
        val5 = np.arange(10000,100000,1000)
        val6 = np.arange(100000,1000000,10000)
        val7 = np.arange(1000000,10000000,1000000)

        val_array = np.concatenate([val1,val2,val3,val4,val5,val6,val7])
        val_array = np.unique(val_array)
        
        #print(val_array)
        #return
    
        x_interp_ = 0. + 2.3*np.pi * (  np.random.normal(0.0,1.0,val_array[-1])  )
        x_interp_ = np.unique(x_interp_)
       
        i_ = val_array.size
        sizepoints = np.zeros(i_)
       
        population = 10
        timearray  = np.zeros([i_,3,population])
        timearray_average  = np.zeros([i_,3])
        standard_deviation = np.zeros([i_,3])
       
        Is_Index = 0 
       
        value = 1
        for i in enumerate(val_array):
           
            for pop in range(0,population):
                value = val_array[i[0]]
                x_interp = x_interp_[0:value]
                
                print("... iteration: {}".format(i[0]))
                print("... Interpolation of array: {0:} {1:} with size == {2:}".format(x_interp[0],x_interp[-1],x_interp.size))
                
                # Calling directly the f2py subroutine
                auxil1time = time.time()
                y_interp,IsKeepOn = interp.lindexerpol(x_interp,x,y,verbosity=verbosity)
                #y_interp, IsKeepOn = interp.akimaspline( x_interp, x, y, delta__x=delta_x, verbosity=verbosity )
               
                # Calling Akima Spline interpolation
                auxil2time = time.time()
                #y_interp_Akima,IsKeepOn = interp1.akimaspline( x_interp, x, y, delta__x=delta_x, verbosity=verbosity )
                scipy_akima_object = interpolate.Akima1DInterpolator( x,y )
                y_interp_scipy = scipy_akima_object( x_interp )
                
                # Cubic Spline using scipy
                auxil3time = time.time()
                cs          = CubicSpline( x,y )
                y_interp_scipy = cs(x_interp)
                
                auxil4time = time.time()
                
                sizepoints[i[0]] = x_interp.size
                
                timearray[i[0],0,pop] = auxil2time-auxil1time
                timearray[i[0],1,pop] = auxil3time-auxil2time
                timearray[i[0],2,pop] = auxil4time-auxil3time
                
                timearray_average[i[0],0]  += timearray[i[0],0,pop]
                timearray_average[i[0],1]  += timearray[i[0],1,pop]
                timearray_average[i[0],2]  += timearray[i[0],2,pop]
                
                print("... Time    LINdexerpol: {0:} s ".format(auxil2time-auxil1time))
                print("... Time Akima   Spline: {0:} s ".format(auxil3time-auxil2time))
                print("... Time Scipy     Spline: {0:} s ".format(auxil4time-auxil3time))
                
            timearray_average[i[0],0]  /= population
            timearray_average[i[0],1]  /= population 
            timearray_average[i[0],2]  /= population 
                
            standard_deviation[i[0],0] = np.std( timearray[i[0],0,:] )
            standard_deviation[i[0],1] = np.std( timearray[i[0],1,:] )
            standard_deviation[i[0],2] = np.std( timearray[i[0],2,:] )
                
            #for i in range(population):
            #    plt.plot(sizepoints,timearray[:,0,i],linewidth=2,color='darkorange',alpha=0.5)
            #    lt.plot(sizepoints,timearray[:,1,i],linewidth=2,color='darkgreen',alpha=0.5)
            #    plt.plot(sizepoints,timearray[:,2,i],linewidth=2,color='#1f77b4',alpha=0.5)

           
        plt.plot(sizepoints,timearray_average[:,0],color='darkorange',label='Time LINdexerpol Fortran', zorder=1)
        #plt.errorbar(sizepoints,timearray_average[:,0], yerr=3.*standard_deviation[:,0],color='orange',ecolor = 'orange', elinewidth = 1, capsize=2, zorder=2)
        #plt.scatter(sizepoints,timearray_average[:,0],color='red', s=2, zorder=3)

        plt.plot(sizepoints,timearray_average[:,1],color='darkgreen',label='Time Akima Scipy', zorder=4)
        #plt.errorbar(sizepoints,timearray_average[:,1], yerr=3.*standard_deviation[:,1],color='green',ecolor = 'lightgreen', elinewidth = 1, capsize=2, zorder=5)
        #plt.scatter(sizepoints,timearray_average[:,1],color='green', s=2, zorder=6)
       
        plt.plot(sizepoints,timearray_average[:,2],color='#1f77b4',label='Time Scipy  Spline')
        #plt.errorbar(sizepoints,timearray_average[:,2], yerr=3.*standard_deviation[:,2],color='blue',ecolor = 'blue', elinewidth = 1, capsize=2, zorder=2)
        #plt.scatter(sizepoints,timearray_average[:,2],color='cyan', s=2, zorder=3)

        plt.xlabel('Number of points in array')
        plt.ylabel('Time [s]')        

        plt.xscale("log")
        plt.yscale("log")

        plt.legend() 
        plt.show()
        
def main():
    i_object = PyLINdexerpol
    i_object.TestRoutine()
    
    #i_object.time()
    
    #x = np.arange(0,2.0*np.pi,0.2)
    #y = np.sin(x**2)
    
    #x_new = [0.5,1.0,2.3,3.22222,4.66,5.2,6.0]
    
    #i_object = Interpolado(x_new,x,y)
    #print( i_object.y_new )

if __name__ == "__main__":
    main()
