#!/usr/bin/python

'''
Module for using the model of Krumholz et al. (2009)
'''

def calc_rh2(h_sd, phi_cnm = None,
        Z = 1.0, # metallicity
        a = 0.2, # ?
        f_diss = 0.1, # fraction of absorbing H2 which disociates
        phi_mol = 10.0, # molecular gas fraction
        mu_H = 2.3e-24, # molecular weight of H, g
        G_0 = 1.0, # Radiation field
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
            / (31.0 * phi_cnm * R_d_solar)

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

def calc_phi_cnm(T_cnm, Z=1.0):

    ''' Calculates phi_cnm from equation 19 in Krumholz et al. (2009).

    Parameters
    ----------
    T_cnm : float, array-like
        Temperature of cold neutral medium in K.
    Z : float, array-like
        Metallicity normalized by solar value.

    '''

    import numpy as np

    T_cnm2 = np.asarray(T_cnm) / 100.0
    Z = np.asarray(Z)

    numerator = 20.0 * T_cnm2**-0.2 * np.exp(1.5 / T_cnm2) * \
                (1 + 3.1 * Z**0.365)

    denominator = 31.0 * (1.0 + 2.6 * (T_cnm2**0.5 * Z)**0.365)

    phi_cnm = numerator / denominator

    return phi_cnm

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

    T_cnm2 = np.arange(1, 300, 0.01) / 100.0
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

def calc_n_min(G_0=1.0, Z=1.0):

    ''' Calculates minimum volume density of CNM to maintain pressure balance
    with the WNM. See equation (5). Returns n_min in cm^-3.

    Parameters
    ----------
    G_0 : float, array-like
        Incident radiation field of FUV photons normalized to solar.
    Z : float, array-like
        Metallicity normalized to solar value.

    '''

    n_min = 31.0 * G_0 / (1 + 3.1 * Z**(0.365))

    return n_min

def calc_n_cnm(G_0=1.0, T_cnm=70.0, Z=1.0):

    ''' Calculates volume density of CNM. See equation (18). Returns n_cnm in
    cm^-3.

    Parameters
    ----------
    G_0 : float, array-like
        Incident radiation field of FUV photons normalized to solar.
    Z : float, array-like
        Metallicity normalized to solar value.

    '''

    import numpy as np

    T_cnm2 = T_cnm / 100.0

    numerator = 20.0 * G_0 * T_cnm2**-0.2 * np.exp(1.5 / T_cnm2)
    denominator = 1 + 2.6 * (T_cnm2**0.5 * Z)**0.365

    n_cnm = numerator / denominator

    return n_cnm

