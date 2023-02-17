#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 18:52:15 2021

@author: jean
"""

import numpy as np
import scipy.interpolate as sinterp
import LINinterpol as interp1
import SPLINE1DArr as interp2
import Interpolado as interp3
import LINdexerpol as interp4

#import interp            as interp5

import matplotlib.pyplot as plt

import time

class Linear__int(object):
    
    def __init__(self,x_new,x,y,e=1.0e-8,verbosity=0):
        self.x = x
        self.y = y
        self.x_new = x_new
        
        self.y_new = interp2.spline1darr(self.x_new,self.x,self.y,e,verbosity)
        
        #print(self.y_new)
    
    def TestRoutine():
        print("... Cerate a sine wave")
        
        x1 = np.arange(0,np.pi,0.2)
        x2 = np.arange(np.pi,2.0*np.pi,0.33)
        #x = np.arange(0,2.0*np.pi,0.2)
        x = np.concatenate([x1,x2])
        y = np.sin(x)

        # ny_value = lininterpol(nx_value,xold_vec,yold_vec,ilastval,[verbosity])
        ilastval = -999
        verbosity = 0
        e         = 1.0e-12
        checking  = 0

        x_new = [x[0]-0.1,x[0],0.5,1.0,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1]

        y_new_interp1, ilatval = interp1.lininterpol(x_new, x, y, ilastval, verbosity)
        y_new_interp2          = interp2.spline1darr(x_new, x, y,        e, verbosity)
        y_new_interp3          = interp3.interpolado(x_new, x, y          , verbosity)
        y_new_interp4          = interp4.lindexerpol(x_new, x, y          , verbosity)
        #y_new_interp5          = interp5.interp_1lin(x_new, x, y, checking, verbosity)
        finterp = sinterp.interp1d(x, y, kind="slinear", fill_value="extrapolate")
        y_new_scipy = finterp(x_new)

        print("... Interpolation of array: {0:}".format(x_new))
        print("... Interpolation LINinterpol: {}".format(np.max(np.abs(y_new_scipy - y_new_interp1))))
        print("... Interpolation SPLINE1DArr: {}".format(np.max(np.abs(y_new_scipy - y_new_interp2))))
        print("... Interpolation Interpolado: {}".format(np.max(np.abs(y_new_scipy - y_new_interp3))))
        print("... Interpolation LINdexerpol: {}".format(np.max(np.abs(y_new_scipy - y_new_interp4))))
        #print("... Interpolation interp_1lin: {}".format(np.max(np.abs(y_new_scipy - y_new_interp5))))

        x_new = [x[0]-3.167,x[0]-2.5,x[0]-0.5,x[0]-0.1,x[0],0.5,1.0,2.3,3.22222,4.66,5.2,6.0,x[-1],x[-1]+0.1,x[-1]+0.5,x[-1]+1.75]

        checking = 0
        ilastval = -999
        y_new_interp1, ilastval = interp1.lininterpol(x_new, x, y, ilastval, verbosity)
        y_new_interp2           = interp2.spline1darr(x_new, x, y,        e, verbosity)
        y_new_interp3           = interp3.interpolado(x_new, x, y          , verbosity)
        y_new_interp4           = interp4.lindexerpol(x_new, x, y          , verbosity)
        #y_new_interp5           = interp5.interp_1lin(x_new, x, y, checking, verbosity)

        # Test with numpy
        #y_new_numpy = np.interp(x_new,x,y,period=2*np.pi)
        finterp = sinterp.interp1d(x,y,kind="slinear",fill_value="extrapolate")
        y_new_scipy = finterp(x_new)

        print("... Interpolation of array: {0:}".format(x_new))
        
        print("... Interpolation LINinterpol: {}".format(np.max(np.abs(y_new_scipy-y_new_interp1))))
        print("... Interpolation SPLINE1DArr: {}".format(np.max(np.abs(y_new_scipy-y_new_interp2))))
        print("... Interpolation Interpolado: {}".format(np.max(np.abs(y_new_scipy-y_new_interp3))))
        print("... Interpolation LINdexerpol: {}".format(np.max(np.abs(y_new_scipy-y_new_interp4))))
        #print("... Interpolation interp_1lin: {}".format(np.max(np.abs(y_new_scipy-y_new_interp5))))

        plt.plot(x,y,linewidth=10)

        plt.plot(x_new, y_new_interp1,color='magenta',linewidth=6)
        plt.plot(x_new, y_new_interp2,color='orange' ,linewidth=9,alpha=0.7)
        plt.plot(x_new, y_new_interp3,color='green'  ,linewidth=5,alpha=0.5)
        plt.plot(x_new, y_new_interp4,color='blue'  ,linewidth=2,alpha=0.5)
        #plt.plot(x_new, y_new_interp5, color='violet', linewidth=2, alpha=0.5)

        #plt.plot(x_new,y_new_numpy,'cyan')
        plt.plot(x_new, y_new_scipy, 'cyan',alpha=0.2)

        plt.scatter(x_new,y_new_interp1,c='red',s=200,alpha=1.0)

        plt.xlabel("radians")
        plt.ylabel("sin(x)")

        plt.show()

    def TestCPU( Nruns=10 ):
        x1 = np.arange(0, np.pi, 0.2)
        x2 = np.arange(np.pi, 2.0 * np.pi, 0.33)
        # x = np.arange(0,2.0*np.pi,0.2)
        x = np.concatenate([x1, x2])
        y = np.sin(x)

        delta = 0.000001
        x_new = np.arange(-1.0,x2[-1]+1.0,delta)

        verbosity = 0
        e = 1.0e-8

        t = np.zeros([7], dtype='float')

        i = np.random.randint(0,x_new.size,)


        start_lininterpol = time.time()
        for i in range(Nruns):
            ilastval = -999
            y_new_interp1, ilastval, IsKeepOn = interp1.lininterpol(x_new, x, y, ilastval, verbosity)
        end_lininterpol = time.time()
        print("... Time LINinterpol: {:12.5}s".format(end_lininterpol-start_lininterpol))
        t[0] = end_lininterpol-start_lininterpol

        start_spline1darr = time.time()
        for i in range(Nruns):
            y_new_interp2, IsKeepOn           = interp2.spline1darr(x_new, x, y, e, verbosity)
        end_spline1darr = time.time()
        print("... Time SPLINE1DArr: {:12.5}s".format(end_spline1darr-start_spline1darr))
        t[1] = end_spline1darr-start_spline1darr

        start_interpolado = time.time()
        for i in range(Nruns):
            y_new_interp3, IsKeepOn           = interp3.interpolado(x_new, x, y, verbosity)
        end_interpolado = time.time()
        print("... Time Interpolado: {:12.5}s".format(end_interpolado - start_interpolado))
        t[2] = end_interpolado - start_interpolado

        start_lindexerpol = time.time()
        for i in range(Nruns):
            y_new_interp4, IsKeepOn           = interp4.lindexerpol(x_new, x, y, verbosity)
        end_lindexerpol = time.time()
        print("... Time LINdexerpol: {:12.5}s".format(end_lindexerpol - start_lindexerpol))
        t[3] = end_lindexerpol - start_lindexerpol

        start_interp_1lin = time.time()
        checking = 0
        #for i in range(Nruns):
        #    y_new_interp4 = interp5.interp_1lin(x_new, x, y, checking, verbosity)
        end_interp_1lin = time.time()
        print("... Time interp_1lin: {:12.5}s".format(end_interp_1lin - start_interp_1lin))
        t[4] = end_interp_1lin - start_interp_1lin

        start_scipy = time.time()
        for i in range(Nruns):
            finterp = sinterp.interp1d(x, y, kind="slinear", fill_value="extrapolate")
            y_new_scipy = finterp(x_new)
        end_scipy = time.time()
        print("... Time scipyinterp: {:12.5}s".format(end_scipy - start_scipy))
        t[5] = end_scipy - start_scipy

        start_numpy = time.time()
        for i in range(Nruns):
            y_new_numpy = np.interp(x_new,x,y)
        end_numpy = time.time()
        print("... Time numpyinterp: {:12.5}s".format(end_numpy - start_numpy))
        t[6] = end_numpy - start_numpy

        print(t[1:]/t[0])

        print(x.size,y.size)

        plt.plot(x, y, linewidth=15)

        plt.plot(x_new, y_new_interp1, color='magenta', linewidth=6)
        plt.plot(x_new, y_new_interp2, color='orange', linewidth=9, alpha=0.7)
        plt.plot(x_new, y_new_interp3, color='green', linewidth=5, alpha=0.5)
        plt.plot(x_new, y_new_interp4, color='blue', linewidth=2, alpha=0.5)
        #plt.plot(x_new, y_new_interp4, color='violet', linewidth=2, alpha=0.5)

        # plt.plot(x_new,y_new_numpy,'cyan')
        #plt.plot(x_new, y_new_scipy, 'cyan', alpha=0.8)

        #plt.scatter(x_new, y_new_interp1, c='red', s=200, alpha=1.0)

        plt.xlabel("radians")
        plt.ylabel("sin(x)")

        plt.show()


def main():
    i_object = Linear__int
    #i_object.TestRoutine()
    
    i_object.TestCPU( Nruns=500 )

if __name__ == "__main__":
    main()
