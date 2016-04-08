#!/usr/bin/python

'''
Module for using the model of Krumholz et al. (2009)
'''

def calc_rh2(h_sd,
        phi_cnm = 3.0,
        Z = 1.0, # metallicity
        a = 0.2, # ?
        f_diss = 0.1, # fraction of absorbing H2 which disociates
        phi_mol = 10.0, # molecular gas fraction
        sigma_d=1.0,
        return_fractions=False,
        remove_helium=True,
        print_nbad=False,
        ):

    '''
    Calculates ratio of molecular hydrogen to atomic hydrogen surface
    density using the model from Krumholz et al. 2008 given a total hydrogen
    surface density.

    Parameters
    ----------
    h_sd : array-like
        Hydrogen surface density in units of solar mass per parsec**2-
    phi_cnm : float
        Ratio of CNM number density and minimum number density to maintain
        pressure balance with the WNM.
    Z : float
        Gas-phase metallicity relative to solar.
    a : float
        ?
    phi_mol : float
        Molecular gas fraction.
    sigma_d : float
        Dust cross section in units of 10^-21 cm^-2 (solar).
    return_fractions : bool
        Return f_H2 and f_HI?
    remove_helium : bool
        Remove contribution from Helium to total gas surface density? If
        True, assumes h_sd does not include Helium.

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

    # caluclate normalized radiation field
    chi = calc_chi(phi_cnm=phi_cnm, Z=Z, sigma_d=sigma_d)

    # dust-adjusted radiation field, EQ 10
    psi = chi * (2.5 + chi) / (2.5 + (chi * np.e))

    # cloud optical depth, EQ 21, include scaling of dust cross-section relative
    # to solar
    tau_c = 0.067 * Z * h_sd * sigma_d

    # remove contribution from Helium?
    if remove_helium:
        tau_c *= 1.4

    # calculate fractions, Equations 35 and 36
    f_H2_sub1 = (3.0 * psi) / (4.0 * tau_c)
    f_H2_sub2 = (4.0 * a * psi * phi_mol) / ((4.0 * tau_c) + (3.0 * (phi_mol \
            - 1.0) * psi))
    f_H2 = 1.0 - (f_H2_sub1 / (1.0 + f_H2_sub2))
    f_HI = 1.0 - f_H2

    # Keep fractions within real fractional value range
    # get number of poor data points:

    if 0:
        mask = (
            (f_HI > 1) | \
             (f_HI < 0) | \
             (f_H2 > 1) | \
             (f_H2 < 0))

        #if print_nbad:
        if np.sum(mask) > 0:
            print('')
            print('number of total points:', np.size(f_HI))
            print('number of f_HI < 0:', np.sum(f_HI < 0))
            print('number of f_HI > 1:', np.sum(f_HI > 1))
            print('number of f_H2 > 1:', np.sum(f_H2 > 1))
            print('number of f_H2 < 0:', np.sum(f_H2 < 0))
    if 1:
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

def calc_chi(phi_cnm=3.0, Z=1.0, sigma_d=1.0, phi_cnm_error=(0.0,0.0),
        Z_error=(0.0,0.0), sigma_d_error=(0.0,0.0)):

    ''' Calculates the normalized radiation field strength, EQ 7 of Krumholz et
    al. (2009).

    Parameters
    ----------
    phi_cnm : float
        Ratio of CNM number density and minimum number density to maintain
        pressure balance with the WNM.
    Z : float
        Gas-phase metallicity relative to solar.
    sigma_d : float
        Dust cross section in units of 10^-21 cm^-2 (solar).
    phi_cnm_error : float
        Uncertainty on phi_cnm.
    Z_error : float
        Uncertainty on Z.
    sigma_d_error : float
        Uncertainty on sigma_d.

    Returns
    -------
    chi : float
        Normalized radiation field strength.
    chi_error : float, optional
        Uncertainty on normalized radiation field strength. Returned if any
        input parameter errors are greater than 0.

    '''

    import numpy as np

    # Constants
    c = 3.0e10 # speed of light, cm/s
    R_d_solar = 10**-16.5 # solar cloud radius, cm
    E_0_solar = 7.5e-4 # solar radiation field, erg/s
    mu_H = 2.34e-24, # molecular weight of H + He, g

    # normalized values
    R_d = R_d_solar / 10**-16.5 * Z # formation rate of H2 on dust, cm^3 / s

    # normalized radiation field strength, EQ 7
    chi = 2.3 * (sigma_d / R_d) * (1 + 3.1 * Z**0.365) / phi_cnm

    if np.sum((phi_cnm_error, Z_error, sigma_d_error)) > 0:
        phi_cnm_error = np.array(phi_cnm_error)
        Z_error = np.array(Z_error)
        sigma_d_error = np.array(sigma_d_error)

        # https://www.wolframalpha.com/input/?i=partial+2.3+*+sigma+%2F+Z+*+(1+%2B+3.1+*+Z%5E0.365)+%2F+phi+with+respect+to+sigma
        sigma_d_comp = 2.3 * (1 + 3.1 * Z**0.365) / phi_cnm * sigma_d_error

        # https://www.wolframalpha.com/input/?i=partial+2.3+*+sigma+%2F+Z+*+(1+%2B+3.1+*+Z%5E0.365)+%2F+phi+with+respect+to+Z
        Z_comp = -(sigma_d *(4.52755 *Z**0.365+2.3))/(Z**2 * phi_cnm) * Z_error

        # https://www.wolframalpha.com/input/?i=partial+2.3+*+sigma+%2F+Z+*+(1+%2B+3.1+*+Z%5E0.365)+%2F+phi+with+respect+to+phi
        phi_cnm_comp =  (sigma_d *(-7.13* Z**0.365-2.3))/(Z *phi_cnm**2) * \
                        phi_cnm_error

        chi_error = (sigma_d_comp**2 + Z_comp**2 + phi_cnm_comp**2)**0.5

        return chi, chi_error
    else:
        return chi


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

def calc_T_cnm(phi_cnm, Z=1.0, phi_cnm_error=(0.0,0.0), calc_error=False):

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

    if calc_error != (0.0, 0.0):

        phi_cnm_low = phi_cnm - phi_cnm_error[0]
        phi_cnm_hi = phi_cnm + phi_cnm_error[1]

        T_cnm2_low = T_cnm2[np.argmin(np.abs(phi_cnm_interp - phi_cnm_low))]

        T_cnm2_hi = T_cnm2[np.argmin(np.abs(phi_cnm_interp - phi_cnm_hi))]

        T_cnm2_error_low = abs((T_cnm2_interp - T_cnm2_low) / \
                               (phi_cnm - phi_cnm_low))  * phi_cnm_error[0]

        T_cnm2_error_hi = abs((T_cnm2_hi - T_cnm2_interp) / \
                              (phi_cnm_hi - phi_cnm)) * phi_cnm_error[1]

        T_cnm_error = 100 * np.array((T_cnm2_error_low, T_cnm2_error_hi))

        return T_cnm, T_cnm_error

    T_cnm = T_cnm2_interp * 100.0

    return T_cnm

def calc_n_min(G_0=1.0, G_0_error=(0.0,0.0), Z=1.0, calc_error=False):

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

    if calc_error:
        n_min_error = 31.0 * G_0_error / (1 + 3.1 * Z**(0.365))

        return n_min, n_min_error

    return n_min

def calc_n_cnm(G_0=1.0, T_cnm=70.0, Z=1.0, G_0_error=(0.0,0.0),
        T_cnm_error=(0.0,0.0), calc_error=False):

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

    if calc_error:
        T_cnm_comp = - 9.49 * np.exp(1.5 / T_cnm2) * G_0 * Z / \
                     (T_cnm2**0.7 * (T_cnm2**0.5 * Z)**0.635 * \
                        (1 + 2.6 * (T_cnm2**0.5 * Z)**0.365 )**2) - \
                     30.0 * np.exp(1.5 / T_cnm2) * G_0 / \
                     (T_cnm2**2.2 * (1 + 2.6 * (T_cnm2**0.5 * Z)**0.365 )) - \
                     4.0 * np.exp(1.5 / T_cnm2) * G_0 / \
                     (T_cnm2**1.2 * (1 + 2.6 * (T_cnm2**0.5 * Z)**0.365 ))

        G_0_comp = 20.0 * T_cnm2**-0.2 * np.exp(1.5 / T_cnm2) / \
                    (1 + 2.6 * (T_cnm2**0.5 * Z)**0.365)

        n_cnm_error = np.sqrt((T_cnm_comp**2 * (T_cnm_error / 100)**2 + \
                               G_0_comp**2 * G_0_error**2))

        return n_cnm, n_cnm_error

    return n_cnm

