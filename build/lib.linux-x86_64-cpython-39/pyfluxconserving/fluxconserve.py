import numpy as np
from scipy.interpolate import interp1d

def fluxconserve(wavelength, flux_density, new_wavelength, kind='linear', fill_val=0.0):
    '''Created on Last version on Wed Sep 23 14:33:51 2020

@author: Jean Gomes Copyright (c)

Python function called fluxconserve that interpolates a given flux density
array onto a new wavelength array and adjusts the flux density values to
conserve the total flux density point-by-point.

The fluxconserve function takes four input arguments:

wavelength: A 1D numpy array containing the wavelengths of the original flux density values.

flux_density: A 1D numpy array containing the flux density values corresponding to the wavelength array.

new_wavelength: A 1D numpy array containing the new wavelengths to interpolate the flux_density array onto.

kind: A string or integer specifying the kind of interpolation to use. The default value is 'linear'.

The function first uses the interp1d function from the scipy.interpolate module to interpolate the flux_density array onto the new new_wavelength array using the specified interpolation kind.

Next, the function calculates a scale factor for each point in the
new_wavelength array to adjust the interpolated flux_density values to
conserve the point-by-point flux density. This is done by dividing the
interpolated flux_density by the original flux_density values at the
corresponding wavelength points. Points where either the interpolated or
original flux_density is not finite are excluded from the scale factor
calculation.

Finally, the function returns the adjusted flux density values using the
scale factor, or a fill_val value for any NaN values.

Overall, this function is useful for conserving the total flux density
when interpolating flux density values onto a new wavelength grid, which
is a common task in many areas of astronomy.'''

    # interpolate the flux density to the new wavelength array
    interp_func = interp1d(wavelength, flux_density, kind=kind, bounds_error=False)
    new_flux_density = interp_func(new_wavelength)

    # interpolate the original flux density to the new wavelength array
    orig_flux_density = interp_func(new_wavelength)

    # adjust the interpolated flux density to conserve point-by-point flux density
    scale = np.ones_like(new_flux_density)
    idx = (np.isfinite(new_flux_density) & np.isfinite(orig_flux_density))
    scale[idx] = new_flux_density[idx] / orig_flux_density[idx]
    return np.where( np.isnan(orig_flux_density * scale), fill_val, orig_flux_density * scale )

# kind : str or int, optional Specifies the kind of interpolation as a string
# or as an integer specifying the order of the spline interpolator to use. The
# string has to be one of ‘linear’, ‘nearest’, ‘nearest-up’, ‘zero’,
# ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, or ‘next’. ‘zero’, ‘slinear’,
# ‘quadratic’ and ‘cubic’ refer to a spline interpolation of zeroth, first,
# second or third order; ‘previous’ and ‘next’ simply return the previous or
# next value of the point; ‘nearest-up’ and ‘nearest’ differ when
# interpolating half-integers (e.g. 0.5, 1.5) in that ‘nearest-up’ rounds up
# and ‘nearest’ rounds down. Default is ‘linear’.

# generate some test data
#wavelength = np.linspace(4000, 8000, 1000)
#flux_density = np.random.normal(size=1000)
#
## define the new wavelength array to interpolate to
#new_wavelength = np.linspace(1000, 12000, 500)
#
## interpolate the flux density to the new wavelength array while conserving point-by-point flux density
#new_flux_density = interpolate_flux(wavelength, flux_density, new_wavelength)
##new_flux_density = np.where( np.isnan(new_flux_density), 0.0, new_flux_density )
#
#
#print(new_wavelength.size,new_wavelength[0],new_wavelength[-1])
#print(new_wavelength.size,new_flux_density[0],new_flux_density[-1])
#
## plot the results
#import matplotlib.pyplot as plt
#plt.plot(wavelength, flux_density, label='Original')
#plt.plot(new_wavelength, new_flux_density, label='Interpolated')
#plt.legend()
#plt.show()
