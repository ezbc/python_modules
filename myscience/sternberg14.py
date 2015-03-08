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

def calc_T_cnm(phi_cnm, Z=1.0):

    ''' Calculates T_cnm from equation 19 in Krumholz et al. (2009).

    Parameters
    ----------
    phi_cnm : float, array-like
        Phi_cnm parameter, n_cnm = phi_cnm * n_cnm,min .
    Z : float, array-like
        Metallicity normalized by solar value.
    G0 : float, array-like
        FUV radiation field normalized by solar value.

    '''

    import numpy as np

    T_cnm2 = np.arange(1, 200, 0.01) / 100.0
    Z = np.asarray(Z)

    numerator = 20.0 * T_cnm2**-0.2 * np.exp(1.5 / T_cnm2) * \
                (1 + 3.1 * Z**0.365)

    denominator = 31.0 * (1.0 + 2.6 * (T_cnm2**0.5 * Z)**0.365)

    phi_cnm_interp = numerator / denominator

    T_cnm2_interp = np.interp(phi_cnm,
                       phi_cnm_interp,
                       T_cnm2
                       )

    T_cnm2_interp = T_cnm2[np.argmin(np.abs(phi_cnm_interp - phi_cnm))]

    T_cnm = T_cnm2_interp * 100.0

    return T_cnm
