#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 12:07:21 2021

RESUME : Interpolation script calling Akima Spline in Fortran!
    
Version: v01 beta

@author: Jean Gomes Copyright (c)

@email: jean@iastro.pt

Written: Jean Michel Gomes © Copyright
"""

# Import libraries
# Be careful with importation of pylab - It needs to be first otherwise problems with other features
from pylab import *
import matplotlib.pyplot as plt

import numpy as np

from PyFluxconserving import flib as interp
import scipy.interpolate as interp1

import scipy.interpolate as interpolate
from scipy.interpolate import CubicSpline
from scipy.interpolate import BSpline

# Class PyAkimaSpline
class PyAkimaSpline( object ):
    """Created on 
Last version on Wed Sep 23 14:33:51 2020

@author: Jean Gomes Copyright (c)

*** Number : 005
1) ABS__smooth
2) Akima_Setup
3) Akimainterp
4) AkimaSpline
5) author_AkimaSpline

RESUME : Class for evaluating Akima Spline.
                  Akima Spline interpolation. The Akima spline is a C1
                  differentiable function (that is, has a continuous first derivative)
                  but, in general, will have a discontinuous second derivative at the
                  knot points.

           [1]: Akima, H. (1970). A New Method of Interpolation and Smooth Curve
                  Fitting Based on Local Procedures. Journal of the ACM, 17(4),
                  589-602. Association for Computing
                  Machinery. doi:10.1145/321607.321609

######################################################################
SUBROUTINE ABS__smooth( x,y,delta_x )
abs__smooth in python 

RESUME : Absolute value function with small quadratic parabola-like
                  in the valley so that it is C1 continuous. C1 continuous
                 - Continuous first derivative.
                                                                            
Input                    arguments = 2                                         
Output                 arguments = 1                                         
Optional              arguments = 0                                         
Total number of arguments = 3                                         
  
  INPUT  : 01) x                -> Old x vector (abcissas)               
                  02) delta_x     -> Valley window point to be quadratic           
                                                                            
OUTPUT : 01) y               -> Output value                                  
                                                                            
PYTHON : Python compatibility using f2py revised. Better usage        
                  with numpy.                                                  
                                                                            
Written: Jean Michel Gomes © Copyright ®                              
Checked: Tue May  1 16:09:13 WEST 2012                                
                  Fri Dec 28 14:55:10 WET  2012                                
######################################################################

######################################################################
SUBROUTINE Akima_Setup( xold_vec,yold_vec,nold_vec,p0,p1,p2,p3,delta__x )
akima_setup in python
    
RESUME : Compute the Spline Coefficients for the Akima Spline         
                  function.                                                    
                                                                           
Input                    arguments = 4                                         
Output                 arguments = 4                                         
Optional              arguments = 0                                         
Total number of arguments = 8                                         
                                                                          
   INPUT  : 01) xold_vec -> Original Abscissa                            
                   02) yold_vec -> Original Ordinate                            
                   03) nold_vec -> # of Points                                  
                   04) delta__x  -> Valley window point to be quadratic          
                                                                           
OUTPUT : 01) p0            -> Spline Coefficient                           
                  02) p1            -> Spline Coefficient                           
                  03) p2            -> Spline Coefficient                           
                  04) p3            -> Spline Coefficient                           

PYTHON : Python compatibility using f2py revised. Better usage        
                  with numpy.                                                  
                                                                           
Written: Jean Michel Gomes © Copyright ®                              
Checked: Wed May  2 10:00:52 WEST 2012                                
                  Fri Dec 28 13:53:36 WET  2012                                
                 Sun Mar 10 10:05:03 WET  2013                                
######################################################################  

######################################################################  
SUBROUTINE Akimainterp( xx_value,yy_value,nx_value,xold_vec,nold_vec,p0,p1,                            &
                                                 p2,p3     ,dp0dxold,dp1dxold,dp2dxold,dp3dxold,                                  &
                                                dp0dyold,dp1dyold,dp2dyold,dp3dyold,dy____dx,                                  &
                                                dy_dxold,dy_dyold,Is_Deriv,Is_Index )

Akimainterp in python

RESUME : Evaluate the Akima Spline and its derivatives.  The          
                  Akima spline is a C1 differentiable unction (that is,        
                  has a continuous first derivative) but, in general,                        
                  will have a discontinuous second derivative at the knot      
                  points.                                                      
                                                                           
 Input                   arguments = 16
 Output                arguments = 4
 Optional              arguments = 1
 Total number of arguments = 21 

   INPUT  : 01) xx_value  -> New interpolated x array with points         
                  02) nxyvalue  -> # of elements in xx_value and yy_value       
                  03) xold_vec   -> Old x vector (abcissas)                      
                  04) nold_vec  -> # of elements in xold_vec and yold_vec       
                  05) p0             -> Spline Coefficient                           
                  06) p1             -> Spline Coefficient                          
                  07) p2             -> Spline Coefficient                          
                  08) p3             -> Spline Coefficient                           
                  09) dp0dxold -> Derivatives of coefficients                  
                  10) dp1dxold -> Derivatives of coefficients                  
                  11) dp2dxold -> Derivatives of coefficients                  
                  12) dp3dxold -> Derivatives of coefficients                  
                  13) dp0dyold -> Derivatives of coefficients                   
                  14) dp1dyold -> Derivatives of coefficients                  
                  15) dp2dyold -> Derivatives of coefficients                  
                  16) dp3dyold -> Derivatives of coefficients                  
                  17) Is_Deriv   -> Is derivative on:1 or off:0                  
                                                                           
OUTPUT : 01) yy_value -> New interpolated yy_value array              
                  02) dy____dx -> Derivative of yy_value w.r.t. xx_value       
                  03) dy_dxold -> Derivative of yy_value w.r.t. yold_vec       
                  04) dy_dyold -> Derivative of yy_value w.r.t. yold_vec       
                                                                           
PYTHON : Python compatibility using f2py revised. Better usage        
                  with numpy.                                                  
                                                                           
Written: Jean Michel Gomes © Copyright ®                              
Checked: Wed May  2 10:00:52 WEST 2012                                
                 Fri Dec 28 13:53:36 WET  2012                                
                 Sun Mar 10 10:05:03 WET  2013                                !
######################################################################  

######################################################################
SUBROUTINE AkimaSpline( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,                                       &
                                                  nold_vec,delta__x,IsKeepOn,verbosity )
akmaspline in python

RESUME : For a given array of values 'xx_value' for the               
              abscissa, then return the ordinate array values of           
              'yy_value' based on Akima Spline interpolation.               
                                                                          
   INPUT  : 01) xx_value  -> New interpolated x array with points        
                  02) nxyvalue  -> # of elements in xx_value and yy_value      
                  03) xold_vec  -> Old x vector (abcissas)                     
                  04) yold_vec  -> Old y vector (ordenadas)                    
                  05) nold_vec  -> # of elements in xold_vec and yold_vec      
                  06) delta__x  -> Value used in the absolute value function   
                                             with small quadratic parabola-like in the   
                                             valley so that it is C1 continuous. C1      
                                             continuous                                  
                  07) verbosity -> Print & Check screen                        
                                                                           
OUTPUT : 01) yy_value -> New interpolated y array with points         
                  02) IsKeepOn -> Flag, if == 0 then there's a problem         
                                                                           
PYTHON : Python compatibility using f2py revised. Better usage        
              with numpy.                                                  
                                                                           
Written: Jean Michel Gomes © Copyright ®                              
Checked: Wed May  2 10:00:52 WEST 2012                                
                 Fri Dec 28 13:53:36 WET  2012                                
                 Sun Mar 10 10:05:03 WET  2013                                
######################################################################  

######################################################################  
SUBROUTINE author_AkimaSpline( a )
author_akimaspline in python
######################################################################  

######################################################################             
                 HOW TO USE
                 
                 i_object = AkimaSpline
                 i_object.AkimaSpline(x_new,x,y,[delta_x],[verbosity])
                 
                 where
                 
                 x_new is the new array of points (abscissa) to be flux-conserved
                 x,y are the original arrays (abscissa and ordinate)
                 
                 [delta_x] => Value used in the absolute value function with small quadratic 
                                       parabola-like in the valley so that it is C1 continuous. 
                                       C1 continuous (continuous first derivative).
                 [verbosity] => verbosity mode [0: No, 1:Yes] 
                 
                 Test ROUTINE
                 
                 i_object = AkimaSpline
                 i_object.TestRoutine()
                 
                 It will produce a plot with the sine interpolation
                 
                 time
                 
                 i_object = AkimaSpline
                 i_object.time()
                 
                 It will test the time it takes to interpolate as a function of elements in the array
######################################################################                  
"""

    def __init__( self,x_new,x,y,delta_x=0.1,tipo=0,verbosity=0 ):
        self.x               = x
        self.y               = y
        self.x_new      = x_new

        self.delta_x = delta_x
        self.verbosity = verbosity
        
        if tipo == 0:
            self.y_new, self.IsKeepOn = interp.akimaspline( self.x_new, self.x, self.y, delta__x=self.delta_x, verbosity=self.verbosity )
            return
        if tipo == 1: # Akima using the Python3 interface
            self.y_new = akima_interp( self.x_new, self.x, self.y )
            return
        if tipo == 2: # Akima using the Scipy interface
            scipy_akima_object = interp1.Akima1DInterpolator(  self.x, self.y )
            self.y_new = scipy_akima_object( self.x_new )
            return
         
        # In case no tipo is chosen, then use tipo == 0
        self.y_new, self.IsKeepOn = interp.akimaspline( self.x_new, self.x, self.y, delta__x=self.delta_x,verbosity=self.verbosity )
        return
        
    def plot( self ):
        Figure = plt.figure( figsize=(12,10),dpi=120,facecolor='w',edgecolor='w' )
        plt.subplots_adjust(bottom=.02, left=.06, right=.95, top=.98, wspace=0.0, hspace=0.0) 
        
        ax = plt.subplot(111)

        # Main plot ##########################################################
        ax1_top = subplot2grid( (20,20), (0,0), colspan=20, rowspan=20 )                                                 # Sets the position and size of the panel for Plot #01
        ax1_top.axis('on')
        ax1_top.axes.get_xaxis().set_visible(True)
    
        ax1_top.plot( self.x,self.y,linewidth=4,color='darkorange' ,label='Original')
        ax1_top.plot( self.x_new,self.y_new,linewidth=3,color='#1f77b4' ,label='Interpolation')
        
        ax1_top.legend()
        
        ax1_top.set_xlabel('x')
        ax1_top.set_ylabel('y(x)')
        # Main plot ##########################################################
                                  
    def TestRoutine( delta_x=0.1 ):
        print("... Create a sine wave function")
        print("... x and y values for x with equal step size")
        
        x = np.arange(0,2.0*np.pi,1.0)
        y = np.sin(x)

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
        
        x_interp          = np.arange(x[0]-0.5,x_new[-1]+0.5,0.01)
        #x_interp = np.array( [-0.7,-0.7,-0.7,-0.7,-0.7099999,0.2,0.3,0.3,0.9,1.2,1.2,1.5,2.3,2.3,3.4,3.5], dtype=float )
        x_new_interp = x_interp
    
        print("... Interpolation of array: {0:}".format([x_interp[0],x_interp[-1]]))

        # Calling directly the f2py subroutine
        y_interp, IsKeepOn = interp.akimaspline( x_interp, x, y, delta__x=delta_x, verbosity=verbosity )
        y_new_interp, IsKeepOn = interp.akimaspline( x_new_interp, x_new, y_new, delta__x=delta_x, verbosity=verbosity )

        #print(y_new_interp)

        # Akima using the Python3 interface
        y_interp_Akima = akima_interp( x_interp, x, y, delta_x=delta_x )
        y_new_interp_Akima = akima_interp( x_new_interp, x_new, y_new, delta_x=delta_x )

        # From scipy Akima interpolator
        scipy_akima_object = interp1.Akima1DInterpolator( x,y )
        y_interp_scipy = scipy_akima_object( x_interp )
        
        # Transform into unique because scipy does not accept repetitive values
        x_new_unique = np.unique(x_new)
        y_new_unique = np.sin(x_new_unique)
        scipy_akima_object = interp1.Akima1DInterpolator( x_new_unique,y_new_unique )
        y_new_interp_scipy = scipy_akima_object( x_new_interp )

        #print(x_new)
        #print(x_new_unique)
        #print(x_new_interp)
        #for i in enumerate(x_new_interp):
         #   print(x_new_interp[i[0]],y_new_interp_Akima[i[0]])
        
        #print(x_new,y_new)
        #print(y_new_scipy)
        #print("... Interpolation of array: {0:}".format(x_new))
        
        Figure = plt.figure( figsize=(12,10),dpi=120,facecolor='w',edgecolor='w' )
        plt.subplots_adjust(bottom=.02, left=.06, right=.95, top=.98, wspace=0.0, hspace=0.0) 
        
        ax = plt.subplot(111)

        # Top plot ###########################################################
        ax1_top = subplot2grid( (20,20), (0,0), colspan=20, rowspan=10 )                                                 # Sets the position and size of the panel for Plot #01
        ax1_top.axis('on')
        ax1_top.axes.get_xaxis().set_visible(False)

        ylimmin = -1.2
        ylimmax = max( y_interp ) * 1.2
        #ax1_top.set_xlim( xmin,xmax )
        ax1_top.set_ylim( ylimmin,ylimmax )

        # Original values
        ax1_top.plot(x,y,color='darkorange',linewidth=10,label='Original')
        
        # Interpolated function Fortran - f2py
        ax1_top.plot(x_interp, y_interp, color='blue',label='Akima Interpolation',linewidth=10,alpha=0.3)

        # Interpolated function Python3 interface
        ax1_top.plot(x_interp, y_interp_Akima, color='magenta',linestyle='--',linewidth=4,label='Akima Python3 Interface')

        # Interpolated function using scipy
        ax1_top.plot(x_interp, y_interp_scipy,color='black',linewidth=1.2,label='Akima Scipy')

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
        ax2_top.plot(x_new_interp, y_new_interp, color='blue',label='Akima Interpolation',linewidth=10,alpha=0.3)
        
        # Interpolated function Python3 interface
        ax2_top.plot(x_new_interp, y_new_interp_Akima, color='magenta',linestyle='--',linewidth=4,label='Akima Python3 Interface')

        # Interpolated function using scipy
        ax2_top.plot(x_new_interp, y_new_interp_scipy,color='black',linewidth=1.2,label='Akima Scipy')
        
        ax2_top.legend()
        ax2_top.set_xlabel("Radians")
        ax2_top.set_ylabel("f(x) = sin(x)")
        
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
        
        value = 1
        for i in enumerate(val_array):
            
            for pop in range(0,population):
                value = val_array[i[0]]
                x_interp = x_interp_[0:value]
        
                print("... iteration: {}".format(i[0]))
                print("... Interpolation of array: {0:} {1:} with size == {2:}".format(x_interp[0],x_interp[-1],x_interp.size))

                # Calling directly the f2py subroutine
                auxil1time = time.time()
                #y_interp,IsKeepOn = interp.spline3darr(x_interp,x,y,verbosity=verbosity)
                y_interp, IsKeepOn = interp.akimaspline( x_interp, x, y, delta__x=delta_x, verbosity=verbosity )
                
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
                
                print("... Time Akima Fortran: {0:} s ".format(auxil2time-auxil1time))
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
        #    plt.plot(sizepoints,timearray[:,1,i],linewidth=2,color='darkgreen',alpha=0.5)
        #   plt.plot(sizepoints,timearray[:,2,i],linewidth=2,color='#1f77b4',alpha=0.5)

            
        plt.plot(sizepoints,timearray_average[:,0],color='darkorange',label='Time Akima Fortran', zorder=1)
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

# Python3 interface for calling the Akima functions
def akima_interp( x, xpt, ypt, delta_x=0.1 ):
    """Convenience method for those who don't want derivatives and don't want to evaluate the same spline multiple times"""
    #print("Python interface")
    
    i = np.argsort(xpt)

    xpt = np.array(xpt[i])
    ypt = np.array(ypt[i])
 
    # Take only unique values 
    # In case they repeat verify if they have the same value and give warning
    i_repeat = np.unique(xpt, return_index=True)[1]

    if i_repeat.size != xpt.size:
        #print("... Houston, we have a problem!")
        xpt = xpt[i_repeat]
        ypt = ypt[i_repeat]

    #print("i_repeat: ",i_repeat)

    p0, p1, p2, p3 = interp.akima_setup( xpt, ypt, delta__x=delta_x )
    #print( p0,p1,p2,p3)

    #y = x*0.0
    #return y       
 
    npt = len(xpt)
    zeros = np.zeros( [npt-1,npt] )

    i_repeat = np.unique(x, return_index=True)[1]
    
    #y = x*0.0
    #return y       
    
    # It should work as well, but the size of the matrices increase the memory exponentially
    #yy_value,dy____dx,dy_dxold,dy_dyold = akimainterp(xx_value,xold_vec,p0,p1,p2,p3   ,dp0dxold,dp1dxold,dp2dxold,dp3dxold,dp0dyold,dp1dyold,dp2dyold,dp3dyold,[is_deriv,is_index])
    #y, dydx, dydxpt, dydypt              = interp.akimainterp( x,            xpt,         p0, p1, p2, p3, zeros,        zeros,      zeros,       zeros,        zeros,       zeros,       zeros,       zeros )

    y, IsKeepOn = interp.akimaspline( x, xpt, ypt, delta__x=delta_x )

    return y

def akima_interp_with_derivs( xpt, ypt, x, delta_x=0.1 ):

    a = Akima( xpt, ypt, delta_x )

    return a.interp(x)

def main():
    i_object = PyAkimaSpline
    i_object.TestRoutine( delta_x=0.1 )
    
    #i_object = PyAkimaSpline    
    #i_object.time()
    
    # x = np.arange(0,2.0*np.pi,0.02)
    # y = np.sin(x**2)
    
    # x_new = [0.5,0.6,0.7,0.8,0.9,1.0,1.0,1.1,1.2,1.22,1.25,1.34,2.3,3.22222,4.66,5.2,5.6,5.7,5.8,6.0]
    
    # i_object = PyAkimaSpline( x_new,x,y,tipo=1 )
        
    # # #print( i_object.y_new )
    
    # o_object = PyAkimaSpline( x,y,delta_x=.400 )
    # a = o_object.interp( x_new )
        
    # #i_object.plot()
    
    # plt.plot( x,y )
    # #plt.plot( x_new,i_object.y_new,color='darkorange',linewidth=2 )
    # plt.plot( x_new,a[0],color='darkgreen',linewidth=2 )
    
if __name__ == "__main__":
    main()
