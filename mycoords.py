#!/usr/bin/python

''' Module for basic coordinate creation.
'''

def make_velocityAxis(h):
    """ Creates the velocity axis given a pyfits header. Assumes the third
    axis is the velocity axis in km/s using the radio definition.
    """

    from numpy import arange

    array = (arange(h['NAXIS3']) - h['CRPIX3'] + 1) * h['CDELT3'] + h['CRVAL3']

    return array / 1000.
















