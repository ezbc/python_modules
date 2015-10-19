#!/usr/bin/python

'''
Module for using the model of Sternberg et al. (2014)
'''

def calc_rh2(h_sd, alphaG=1.5, Z=1.0, phi_g=1.0, sigma_g21=1.9,
        radiation_type='isotropic', return_fractions=True):

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
        hi_sd = 6.7 / (Z * phi_g) * np.log(alphaG / 3.2 + 1) # Msun pc^-2
        if 0:
            hi_sd = 6.71 * (1.9 / sigma_g21) * \
                    np.log(alphaG / 3.2 + 1) # Msun pc^-2

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

def calc_n_H(I_UV=1, alphaG=1, phi_g=1.0, Z_g=1.0, Z_CO=1.0):

    ''' Equation 5 in Bialy et al. (2015)
    '''

    phi_cnm = calc_phi_cnm(alphaG=alphaG, Z_g=Z_g, Z_CO=Z_CO, phi_g=phi_g)

    n_H = 22.7 * I_UV * 4.1 / (1 + 3.1 * Z_g**0.365) * Z_g / Z_CO * phi_cnm / 3.

    return n_H

def calc_phi_cnm(alphaG=1.0, Z_g=1.0, Z_CO=1.0, phi_g=1.0,):

    phi_cnm = 2.58 * (1 + 3.1 * Z_g**0.365) / 4.1 * Z_CO / Z_g * \
              3.0 / alphaG * 2.62 / (1 + (2.63 * phi_g * Z_g)**0.5) * phi_g

    return phi_cnm

