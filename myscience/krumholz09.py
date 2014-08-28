#!/usr/bin/python

'''
Module for using the model of Krumholz et al. (2009)
'''

def calc_rh2(h_sd, phi_CNM = None,
        Z = 1.0, # metallicity
        a = 0.2, # ?
        f_diss = 0.1, # fraction of absorbing H2 which disociates
        phi_mol = 10.0, # molecular gas fraction
        mu_H = 2.3e-24, # molecular weight of H, g
        return_fractions=False
        ):
    '''
    Calculates ratio of molecular hydrogen to atomic hydrogen surface
    density using the model from Krumholz et al. 2008 given a total hydrogen
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

    # Constants
    c = 3.0e10 # speed of light, cm/s

    # solar values
    sigma_d_solar = 1e-21 # solar dust grain cross section, cm^2
    R_d_solar = 10**-16.5 # solar cloud radius, cm
    E_0_solar = 7.5e-4 # Radiation field, erg/s

    # cloud values
    sigma_d = sigma_d_solar * Z # dust grain cross section, cm^2
    R_d = R_d_solar * Z # cloud radius, cm

    # normalized radiation field strength, EQ 7
    chi = ((f_diss * sigma_d_solar * c * E_0_solar) \
            * (1.0 + (3.1 * Z**0.365))) \
            / (31.0 * phi_CNM * R_d_solar)

    # dust-adjusted radiation field, EQ 10
    psi = chi * (2.5 + chi) / (2.5 + (chi * np.e))

    # cloud optical depth, EQ 21
    tau_c = (3.0 * h_sd * sigma_d) / (4.0 * (3.1 * Z**0.365) * mu_H)

    tau_c = (3.0 * h_sd * 2.0 * 10.0**33 * sigma_d) / \
            (4.0 * (3.1 * 10**18)**2 * mu_H)

    # cloud optical depth, EQ 21
    tau_c = 0.067 * Z * h_sd

    f_H2_sub1 = (3.0 * psi) / (4.0 * tau_c)
    f_H2_sub2 = (4.0 * a * psi * phi_mol) / ((4.0 * tau_c) + (3.0 * (phi_mol \
            - 1.0) * psi))

    f_H2 = 1.0 - (f_H2_sub1 / (1.0 + f_H2_sub2))
    f_HI = 1.0 - f_H2

    # Keep fractions within real fractional value range
    f_HI[f_HI > 1] = 1.0
    f_HI[f_HI < 0] = 0.0
    f_H2[f_H2 > 1] = 1.0
    f_H2[f_H2 < 0] = 0.0

    # ratio of molecular to atomic fraction, EQ 17 Lee et al. 2012
    R_H2 = 4 * tau_c / (3 * psi) \
            * (1+ 0.8 * psi * phi_mol \
                / (4 * tau_c + 3 * (phi_mol - 1) * psi)) -1

    R_H2 = f_H2 / f_HI

    if not return_fractions:
        return R_H2
    elif return_fractions:
        return R_H2, f_H2, f_HI


