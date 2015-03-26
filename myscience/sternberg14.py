#!/usr/bin/python

'''
Module for using the model of Sternberg et al. (2014)
'''

def calc_rh2(h_sd, alphaG=1.5, Z=1.0, phi_g=1.0, radiation_type='isotropic',
        return_fractions=True):

    '''
    Calculates ratio of molecular hydrogen to atomic hydrogen surface
    density using the model from Sternberg et al. (2014) given a total hydrogen
    surface density.

    Parameters
    ----------
    h_sd : array-like
        Hydrogen surface density in units of solar mass per parsec**2-
    return_fractions : bool
        Return f_H2 and f_HI?

    Returns
    -------
    rh2_fit : array-like
        Model ratio between molecular and atomic hydrogen masses.
    f_H2, f_HI : array-like, optional
        If return_fractions=True then the following is returned:
        f_H2 = mass fraction of molecular hydrogen
        f_HI = mass fraction of atomic hydrogen
    '''

    import numpy as np

    if radiation_type == 'isotropic':
        hi_sd = 9.5 / (Z * phi_g) * np.log(alphaG / 3.2 + 1) # Msun pc^-2
    elif radiation_type == 'beamed':
        hi_sd = 11.9 / (Z * phi_g) * np.log(alphaG / 2.0 + 1) # Msun pc^-2
    else:
        raise ValueError('radiation_type must be "beamed" or "isotropic"')

    f_H2 = 1.0 - hi_sd / h_sd
    f_HI = 1.0 - f_H2

    # Keep fractions within real fractional value range
    f_HI[f_HI > 1] = 1.0
    f_HI[f_HI < 0] = 0.0
    f_H2[f_H2 > 1] = 1.0
    f_H2[f_H2 < 0] = 0.0

    # ratio of molecular to atomic fraction
    R_H2 = f_H2 / f_HI

    if not return_fractions:
        return R_H2
    elif return_fractions:
        return R_H2, f_H2, f_HI

