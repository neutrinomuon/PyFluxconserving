#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Revised interface on Sat Jan 30 12:07:21 2021

@author: Jean Gomes

RESUME : FluxConSpec subroutine is written in fortran 2003 standard. Flux -conserving 
                  interpolation script.

Version: v01 beta

PYTHON : Python compatibility using f2py revised. Better usage  with numpy.  

Written: Jean Michel Gomes © Copyright
"""

# Import libraries
# Be careful with importation of pylab - It needs to be first otherwise problems with other features
from pylab import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pylab as pl

import numpy as np
import random
import time

import sys
sys.path.insert(1, '/home/jean/Codes/Pynoptic/Interpolation/')
sys.path.insert(1, '/home/jean/Codes/Pynoptic/ColorMaps/')

import FluxConSpec as interp
import spectral_resampling_adam as sr
from spectres import spectral_resampling_numba as srn
import califa_cmap_alternative as califa_cmap

import scipy.interpolate as interpolate
from scipy.interpolate import CubicSpline
from scipy.interpolate import BSpline

class FluxConSpec(object):
    '''Created on 
Last version on Wed Sep 23 14:33:51 2020

    @author: Jean Gomes Copyright (c)

    *** Number : 002
    1) FluxConSpec
    2) author_FluxConSpec
    
######################################################################
SUBROUTINE FluxConSpec( Orlambda,Orfluxes,Nrlambda,O_lambda,O_fluxes,                               &
                        N_lambda,per_bins,slow_int,IsKeepOn,fill_val,                                                                 &
                        verbosity )
fluxconspec in  python

RESUME : Rebinning of input vector on a new grid and returns the
                  output array. This routine uses the derivative of a          
                  cumulative function to guarantee the conservation of        
                  density flux while rebinning the spectrum.                   
                                                                           
     OBS.  :  Central wavelengths must be specified to compute the         
                  wavelengths at the edges of the pixels. This routine         
                  assumes a monotonically increasing function for 'X',         
                  that can be unequally spaced.                                
                                                                           
                  The Cumulative function has several ways to be               
                  interpolated:                                                
                                                                           
                  Interpolation Schemes:                                       
                  00) slow_int -> LINinterpol                                  
                  01) slow_int -> SPLINE3DFor                                  
                  02) slow_int -> SPLINE1DArr                                  
                  03) slow_int -> pchipmodule                                  
                  04) slow_int -> AkimaSpline                                  
                  05) slow_int -> Interpolado                                  
                  06) slow_int -> LINdexerpol                                  
                  07) slow_int -> poly_interp
    
Input                    arguments = 7 
Output                 arguments = 2
Optional              arguments = 2
Total number of arguments = 11
   
  INPUT  : 01) Orlambda  -> New 'X' vector                      
                 02)  Nrlambda  -> Number of elements in new 'X' vector       
                 03) O_lambda   -> Old 'X' vector (Abscissa)                    
                 04) O_fluxes      -> Old 'Y' vector (Ordinate)
                 05) N_lambda   -> Number of elements in old 'X' vector         
                 06) per_bins     -> [0: Not conserve flux , 1: Conserve flux  ]  
                 07) slow_int      -> [0 -- 7]                                     
                 08) fill_val         -> Value to fill values outside boundary             
                 09) verbosity   -> Optional variable to print & check           
                                                                           
OUTPUT: 01) Orfluxes -> New 'Y' vector                               
                 02) IsKeepOn -> Flag, if == 0 then there's a problem         
                                                                           
PYTHON: Python compatibility using f2py revised. Better usage        
                  with numpy.                                                  
                                                                           
EXTRA ROUTINES : Those mentioned in slow_int interpolation scheme.    
                                                                           
       LOG : Modified to take into account Y negative values                  
                  and minor corrections in the recovered Y values                  
                  Two different corrections: One related to eps values in the      
                  cumulative function and another one related to the cumulative    
                  being set erroneously to zero.                                   
                                                                           
Checked: Sat Dec  8 12:03:44 WET  2012                                
                 Fri Dec 28 15:34:27 WET  2012                                
                 Wed Apr 29 18:13:43 WEST 2015                                
                 Thu Oct 20 21:55:02 WEST 2016
######################################################################             

######################################################################
SUBROUTINE author_FluxConSpec( a )
author_fluxconspec in python
######################################################################             

######################################################################             
                 HOW TO USE
                 
                 i_object = FluxConSpec
                 i_object.FluxConSpec(x_new,x,y,[per_bins],[slow_int],[fill_val],[verbosity])
                 
                 where
                 
                 x_new is the new array of points (abscissa) to be flux-conserved
                 x,y are the original arrays (abscissa and ordinate)
                 
                 [per_bins]  => optional - [0: Not conserve flux , 1: Conserve flux  ]
                 [slow_int]   => optional - (Type of interpolation - see above)
                 [fill_val]      => optional - Value to fill values outside boundary
                 [verbosity] => verbosity mode [0: No, 1:Yes] 
                 
                 Test ROUTINE
                 
                 i_object = FluxConSpec
                 i_object.TestRoutine()
                 
                 It will produce a plot with the different options for interpolation - Flux-conserving
                 
                 time
                 
                 i_object = FluxConSpec
                 i_object.time()
                 
                 It will test the time it takes to interpolate as a function of elements in the array
######################################################################
'''
    
    def __init__(self,x_new,x,y,per_bins=1,slow_int=0,fill_val=0.0,verbosity=0):
        self.x = x
        self.y = y
        self.x_new = x_new
        
        self.slow_int = slow_int
        self.fill_val = fill_val
        
        self.y_new,self.IsKeepOn =  interp.fluxconspec(self.x_new,self.x,self.y,per_bins,slow_int,fill_val,verbosity)
        
        #print(self.y_new)
        return
    
    def plot( self ):
        string_f = "slow_int = {0:}".format(self.slow_int)
        
        plt.plot(self.x,self.y,label='Original Values')
        plt.plot(self.x_new,self.y_new,label='Interpolated Values - ' + string_f)
        
        plt.xlabel('x')
        plt.ylabel('y')
        
        plt.legend()
        plt.title('Original versus New Interpolated Values')
        plt.show()
    
    def TestRoutine():
        print("... Create a sine wave function")
        x1 = np.arange(-0.5,np.pi/2.0,0.23)
        x2 = np.arange(np.pi/2.0,np.pi,0.010)
        x3 = np.arange(np.pi,2.0*np.pi+0.8+0.2,0.080)

        np.random.seed(100)
        x1 = -0.5  + (0.23 + 0.5) *  np.random.rand(10)
        x1 = np.sort(x1)

        #print("x1",x1)
        z = 5.2
        x = np.concatenate([x1,x2,x3])
        x = x / (1. + z)
        #x = np.concatenate([x1,x2])
        y = np.where( (x>=0.0) & (x <= 2.0*np.pi), np.sin(x)**4 + 2., 1.0)
        
        #print(x,y)
        #print()
        #ind = np.arange(0,np.size(x),1)
        #random.shuffle(ind)
        
        #x = x[ind]
        #y = y[ind]
        
        # ny_value = lininterpol(nx_value,xold_vec,yold_vec,ilastval,[verbosity])
        verbosity = 0
        IsKeepOn = 1
        slow_int = 0
        per_bins = 1
        
        fill_val = 1.256
        
        #x_new = [x[0]-0.1,x[0],0.5,1.0,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1]
        #y_new = interp.fluxconspec(x_new,x,y,per_bins,slow_int,IsKeepOn,verbosity)
        
        #print("... Interpolation of array: {0:}".format(x_new))
        
        x_new_ = np.arange(6.5,8.0,0.001)
        x_new_1 = np.arange(4.67,5.2,0.0123)

        x_new = [x[0]-0.7,x[0]-0.1,x[0],0.5,1.0,2.3,np.pi,3.22222,4.66,5.2,6.0,6.2]#,6.5]#,7.0,7.2,7.5,8.0] #,x[-1]+0.5,x[-1]+1.0,x[-1]+2.0]
        x_new = np.concatenate([x_new,x_new_1,x_new_])
        #x_new = np.arange(x[0]-1.0,x[-1]+1.0,0.3)
        x_new = np.array(x_new, dtype=float)
        x_new = np.unique(x_new)
        #x_new = np.arange(x[0]-0.5,x[-1]+0.5,0.5)
        
        begin_time = time.time()
        #print(fill_val)

        y_new_,IsKeepOn = interp.fluxconspec(x_new, x, y, per_bins, slow_int, fill_val=fill_val, verbosity=verbosity)

        auxil1time = time.time()
        #print(x.size,x)
        model_resampled = sr.spectres(x_new, x, y, fill=fill_val)
        #model_resampled = y_new
        auxil2time = time.time()
        
        print("... Time fluxconspec: {}s and spectral_resampling: {}s - Ratio fluxconspec/sample_resampling {}".format(auxil1time-begin_time,auxil2time-auxil1time,(auxil1time-begin_time)/(auxil2time-auxil1time)))
        print()
        
        t_new = np.linspace(x_new[0],x_new[-1], 31)
        h_new = np.interp(t_new,x_new,y_new_)
        
        int0 = np.trapz(y_new_,x=x_new)
        int1 = np.trapz(h_new,x=t_new,dx=(t_new[1]-t_new[0]))
        int2 = np.trapz(model_resampled,x=x_new)
        int3 = np.trapz(y,x=x)
 
        print("... Total areas in the resampling: Flux Conserving {} | Linear Interpolation {} | Spectral Resampling {} | Original Data {}".format(int0,int1,int2,int3))
        print("")
        
        # Types of interpolation used in the cumulative curve
        interpolation_names = ['Linear  Table Indexing    ' , 'Cubic Spline                   ', 'One Dimensional Spline ', 'Hermite Interpolation     ', 'Akima Interpolation        ', 'Simple Interpolation       ','Indexer Interpolation      ', 'Polynomial Interpolation ']
        N_slows = len(interpolation_names)
        
        # Compute Areas
        i_to_integrate = np.where( (x_new >= x[0]) & (x_new <= x[-1]) )
        
        i_to_integrate_original = np.where( (x >= x_new[0]) & (x <= x_new[-1]) )
        A_original = np.trapz(y[i_to_integrate_original],x=x[i_to_integrate_original])
        
        #for i in enumerate(d):
        # y_to_integrate = d[i[1]]
        # A = np.trapz(y_to_integrate[i_to_integrate],x=x_new[i_to_integrate])
        # print(A,int3,A/int3)

        d = {}
        for i in range(0,N_slows):
            slow_int = i
            #print(i,string_f)
            IsKeepOn = 1
            y_new,IsKeepOn = interp.fluxconspec(x_new,x,y,per_bins,slow_int,fill_val=fill_val,verbosity=verbosity)

            if IsKeepOn == 0:
                print("Deu Merda {}".format(slow_int))

            #if i == 3:
            #    print(x_new)
            #    print(y_new)

            y_to_integrate = y_new
            A = np.trapz(y_to_integrate[i_to_integrate],x=x_new[i_to_integrate])
            #print(A,int3,A/int3)
            
            string_ = "{0:} | A =".format(interpolation_names[i])
            #print(len(string_))
            
            string_f = "y_new_{0:} = {1:} {2:<12.6f}".format(i,string_,A)
            d[string_f] = y_new,A

        # Compute Area from Spectral Resampling
        A = np.trapz(model_resampled[i_to_integrate],x=x_new[i_to_integrate])
        #print(A,int3,A/int3)
        string_f = "Spectral Resampling Python - Adam | A = {0:<12.6f}".format(A)
        d[string_f] = model_resampled,A

        #print(d)
        #print("... Interpolation of array:")
        #print(x_new)
        print()
        #print(d)
        
        Figure = plt.figure( figsize=(12,10),dpi=120,facecolor='w',edgecolor='w' )
        plt.subplots_adjust(bottom=.02, left=.06, right=.95, top=.98, wspace=0.0, hspace=0.0) 
        
        ax = plt.subplot(111)

        # Top plot ###########################################################
        ax1_top = subplot2grid( (20,20), (0,0), colspan=20, rowspan=20 )                                                 # Sets the position and size of the panel for Plot #01
        ax1_top.axis('on')
        ax1_top.axes.get_xaxis().set_visible(True)
        
        xmin  = -1.8
        xmax = +10.7
        
        #xmin=7
        #xmax=7.25
        
        ylimmin = -1.2
        ylimmax = +5.2
        ax1_top.set_xlim( xmin,xmax )
        ax1_top.set_ylim( ylimmin,ylimmax )

        # Original values
        ax1_top.plot(x,y,color='darkorange',linewidth=3,label='Original A = {}'.format(A_original),zorder=1)
        ax1_top.scatter(x,y,color='green',s=100.0,label='Original',zorder=2)
        
        # Interpolated function Fortran - f2py
        #ax1_top.plot(x_interp, y_interp, color='blue',label='Spline 3D Interpolation',linewidth=10,alpha=0.3)

        cmap =  califa_cmap.get_califa_velocity_cmap()
        norm = mpl.colors.Normalize(vmin=0, vmax=len(d)-1)
        colors = cmap(norm(np.arange(0, len(d)+1, 1)))
        for i in enumerate(d):
            #print(i[1])
            
            if i[0] < len(d)-1:
                ax1_top.plot(x_new,d[i[1]][0]+0.3*i[0],color=colors[i[0]],linewidth=3,label='{}'.format(i[1]),zorder=i[0]+2)
            else:
                ax1_top.plot(x_new,d[i[1]][0]+0.0*i[0],color='purple',linewidth=3,linestyle='--',label='{}'.format(i[1]),zorder=i[0]+2)

        ax1_top.fill_between(x_new, y_new_, 0.0, facecolor='orange', alpha=0.6, edgecolor='white', hatch='/',zorder=0)

        # Interpolated function Python3 interface
        #ax1_top.plot(t_new, h_new, color='green',linestyle='--',linewidth=4,label='Numpy Python3 Interface')

        # Interpolated function using scipy
        #ax1_top.plot(x_interp, y_interp_scipy,color='black',linewidth=1.2,label='Akima Scipy')

        # Spectral Resampling
        #plt.plot(x_new, model_resampled, color='purple',linewidth=3,linestyle='--',label='Spectral Resampling Python - Adam')

        # Plot points corresponding to the x values
        #ax1_top.scatter(x,y*0,color='cyan',s=100.0)
        #ax1_top.scatter(x_new,y_new*0,color='magenta',s=30.0)

        for i in enumerate(x):
             x2 = [x[i[0]],x[i[0]]]
             y2 = [0,y[i[0]]]
             ax1_top.plot( x2,y2, marker = 'o',color='green' )

        for i in enumerate(x_new):
            x1 = [x_new[i[0]],x_new[i[0]]]
            y1 = [0,y_new_[i[0]]]
            
            ax1_top.plot( x1,y1, marker = 'o',color='black',alpha=0.3 )

        #get handles and labels
        handles, labels = plt.gca().get_legend_handles_labels()

        #print(handles)
        #print(labels)

        #specify order of items in legend
        order1 = np.array( [0,1], dtype=int )
        order2 = np.arange(2,N_slows+3,1)
        order = np.concatenate( [order1,order2], dtype=int )
        
        #add legend to plot
        ax1_top.legend( [handles[idx] for idx in order],[labels[idx] for idx in order],title = "Flux Conserving Interpolations",fancybox=True, shadow=True,loc='upper right' ).set_zorder(1200)

        #ax1_top.legend()
        ax1_top.set_xlabel("Radians")
        ax1_top.set_ylabel("f(x) = sin(x)")
        # Top plot ###########################################################
        
        return
    
    def time():
        print("... Create a sine wave function")
        print("... x and y values for x with equal step size")
     
        x = np.arange(0,4.3*np.pi,0.1)
        y = np.sin(x)
     
        print("... step x: {0:}".format(x[1]-x[0]))

        # No verbosity 
        verbosity = 0
     
        delta_x = 0.1

        #i_ = 25

        #sizepoints = np.zeros(i_)
        #timearray  = np.zeros([i_,3])

        np.random.seed(100)
        
        #x_interp_ = np.random.normal(0.0,1.0,2**i_)
        #x_interp_ = np.unique(x_interp_)
        
        val1 = np.arange(2,10,1)
        val2 = np.arange(10,100,2)
        val3 = np.arange(100,1000,10)
        val4 = np.arange(1000,10000,100)
        val5 = np.arange(10000,100000,1000)
        val6 = np.arange(100000,1000000,10000)
        val7 = np.arange(1000000,10000000,1000000)

        val_array = np.concatenate([val1,val2,val3,val4,val5])#,val6,val7])
        val_array = np.unique(val_array)
        val_array = np.sort(val_array)
    
        x_interp_ = 0. + 2.3*np.pi * (  np.random.normal(0.0,1.0,val_array[-1])  )
        #x_interp_ = np.arange(0,2.3*np.pi,0.001)
        x_interp_ = np.unique(x_interp_)
        x_interp_ = np.sort(x_interp_)
       
        per_bins = 1
        fill_val = -2.3
       
        population = 30
        
        i_ = val_array.size
        sizepoints = np.zeros(i_)
        
        interpolation_names = ['Linear  Table Indexing    ' , 'Cubic Spline                   ', 'One Dimensional Spline ', 'Hermite Interpolation     ', 'Akima Interpolation        ', 'Simple Interpolation       ', 'Indexer Interpolation      ','Polynomial Interpolation','Akima  Spline Scipy', 'Cubic Spline Scipy','Spectral Resampling','Spectral Resampling Numba','Numpy Interpolation'] 
        N_slows = len(interpolation_names)
        
        colors          = pl.cm.jet(np.linspace(0,1,N_slows))
        
        timearray  = np.zeros([i_,N_slows,population])
        timearray_average  = np.zeros([i_,N_slows])
        standard_deviation = np.zeros([i_,N_slows]) 
        
        for i in enumerate(val_array):
            value = val_array[i[0]]
            x_interp = x_interp_[:value]
            print(x_interp[:value].size)
            sizepoints[i[0]] = x_interp.size         
              
            for j in range(0,N_slows-3):
                for pop in range(0,population): 
                    slow_int = j
                    auxil1time = time.time()
                    y_interp,IsKeepOn = interp.fluxconspec(x_interp,x,y,per_bins,slow_int,fill_val=fill_val,verbosity=verbosity)
                    auxil2time = time.time()

                    timearray[i[0],j,pop] = auxil2time-auxil1time
                    timearray_average[i[0],j]  += timearray[i[0],j,pop]
                    #print("{}".format(auxil2time-auxil1time))
            
                timearray_average[i[0],j]  /= population
                standard_deviation[i[0],j] = np.std( timearray[i[0],j,:] )
                
            for pop in range(0,population): 
                x_interp = x_interp_[:value]

                # Calling Akima Spline interpolation
                auxil1time = time.time()
                #y_interp_Akima,IsKeepOn = interp1.akimaspline( x_interp, x, y, delta__x=delta_x, verbosity=verbosity )
                scipy_akima_object = interpolate.Akima1DInterpolator( x,y )
                y_interp_scipy = scipy_akima_object( x_interp )
                auxil2time = time.time()

                # Cubic Spline using scipy
                auxil3time = time.time()
                cs          = CubicSpline( x,y )
                y_interp_scipy = cs(x_interp)
                auxil4time = time.time()
                
                # Spectres
                auxil5time = time.time()
                model_resampled = sr.spectres(x_interp, x, y, fill=fill_val)
                auxil6time = time.time()
                
                # Spectres Numba
                auxil7time = time.time()
                model_resampled = srn.spectres_numba(x_interp, x, y, fill=fill_val)
                auxil8time = time.time()
                
                # Numpy Interpolation
                auxil9time = time.time()
                numpy_interp = np.interp(x_interp, x, y)
                auxil10time = time.time()
                
                j = N_slows-5
                timearray[i[0],j,pop] = auxil2time-auxil1time
                timearray_average[i[0],j]  += timearray[i[0],j,pop]
                
                j = N_slows-4
                timearray[i[0],j,pop] = auxil4time-auxil3time
                timearray_average[i[0],j]  += timearray[i[0],j,pop]

                j = N_slows-3
                timearray[i[0],j,pop] = auxil6time-auxil5time
                timearray_average[i[0],j]  += timearray[i[0],j,pop]

                j = N_slows-2
                timearray[i[0],j,pop] = auxil8time-auxil7time
                timearray_average[i[0],j]  += timearray[i[0],j,pop]
                
                j = N_slows-1
                timearray[i[0],j,pop] = auxil10time-auxil9time
                timearray_average[i[0],j]  += timearray[i[0],j,pop]

            timearray_average[i[0], N_slows-5]  /= population
            timearray_average[i[0], N_slows-4]  /= population
            timearray_average[i[0], N_slows-3]  /= population
            timearray_average[i[0], N_slows-2]  /= population
            timearray_average[i[0], N_slows-1]  /= population
            
            standard_deviation[i[0], N_slows-5] = np.std( timearray[i[0], N_slows-5,:] )
            standard_deviation[i[0], N_slows-4] = np.std( timearray[i[0], N_slows-4,:] )
            standard_deviation[i[0], N_slows-3] = np.std( timearray[i[0], N_slows-3,:] )
            standard_deviation[i[0], N_slows-2] = np.std( timearray[i[0], N_slows-2,:] )
            standard_deviation[i[0], N_slows-1] = np.std( timearray[i[0], N_slows-1,:] )

            
        for j in range(0,N_slows):
            #if j==N_slows-1:
            #    plt.plot(sizepoints[:],timearray_average[:,j],label='Time : {0:} {1:}'.format(j,interpolation_names[j]), zorder=1,color='black')
            #else:
            plt.plot(sizepoints[:],timearray_average[:,j],label='Time : {0:} {1:}'.format(j,interpolation_names[j]), zorder=1,color=colors[j])

                
        plt.xlabel('Number of points in array')
        plt.ylabel('Time [s]')        

        plt.xscale("log")
        plt.yscale("log")

        plt.legend( fontsize=7 ) 
        plt.show() 
        
        return 

def main():
    i_object = FluxConSpec
    
    # TestRoutine
    #i_object.TestRoutine()
    
    # time
    i_object.time()
    
    #x = np.arange(0,2.0*np.pi,0.2)
    #y = np.sin(x**2)
    
    #x_new = [0.5,1.0,2.3,3.22222,4.66,5.2,6.0]
    
    #i_object = FluxConSpec(x_new,x,y,slow_int=4)
    #i_object.plot()
    #print( i_object.y_new )

    # plt.plot(x,y,color='darkorange')
    # plt.plot(x_new,i_object.y_new)

if __name__ == "__main__":
    main()
