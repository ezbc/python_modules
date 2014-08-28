#!/usr/bin/python

import numpy as np
from astropy.io import fits as pf
from mycoords import make_velocity_axis

def calc_moment(data = None, moment = 2, header = None, spectral_axis = 2):

    ''' Calculates a moment of an image along the spectral axis.
    '''

    if data is None:
    	raise ValueError('No data provided.')
    if type(data) is not np.ndarray:
    	raise ValueError('Data must be of type numpy.ndarray')
    if len(data.shape) != 3:
    	raise ValueError('Data must have 3 dimensions.')

    try:
        n = data.shape[spectral_axis]
    except IndexError:
        raise IndexError('Spectral axis not in data.')

    # Avoid Nans
    data = np.ma.array(data, mask = (data != data))

    # Moments about mean
    mean = data.mean(axis = spectral_axis)

    if moment == 0:
        mom_image = data.sum(axis = spectral_axis)
    if moment == 1:
        mom_image = (data - mean).sum(axis = spectral_axis) / n
    if moment == 2:
        mom_image = ((data.sum(axis = spectral_axis) - mean)**2) / n

    return mom_image

def polyfit2d(x, y, z, degree=1, p0=None):

    ''' 2D polynomial fit to z at the given x and y coordinates. z(x,y) =
    p[0] * x**deg + p[1] * y**deg + ... + p[deg].
    '''

    from scipy.optimize import curve_fit

    poly2d = make_poly2d_fit(degree=degree)

    def poly2d((x, y), a, b, c):
        coeffs = (a,b,c)
        z = coeffs[0] * np.ones(x.shape)

        z += coeffs[1] * x**degree + coeffs[2] * y**degree

        return z

    print poly2d((x,y), 0, 1, 1,)

    popt, pcov = curve_fit(poly2d, (x, y), z, p0=None)

    return popt, pcov

def poly2d((x, y), degree=1, *coeffs):

    if len(coeffs) != 2 * degree + 1:
        raise ValueError('Coeffslength must be 2*degree + 1')
    else:
    	pass

    z = coeffs[0] * np.ones(x.shape)

    for i in xrange(1, degree + 1, 2):
        z += coeffs[i] * x**i + coeffs[i + 1] * y**i

    return z.ravel()

def make_poly2d_fit(degree=1):

    ''' Make a function with 2 * degree + 1 coefficients as individual
    arguments. Useful for supplying to scipy.optimize.curve_fit.

    '''

    def poly2d((x, y), *coeffs):
        try:
            if x.shape == y.shape:
            	pass
            else:
            	raise ValueError('x and y must be the same shape.')
        except AttributeError:
            raise ValueError('x and y must be numpy.ndarray')

        z = coeffs[0] * np.ones(x.shape)

        for i in xrange(1, degree + 1, 2):
            z += coeffs[i] * x**i + coeffs[i + 1] * y**i

        return z

    return poly2d





