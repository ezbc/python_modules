#!/usr/bin/python

import numpy as np

'''
Module for using the model of Sternberg et al. (2014)
'''

def calc_rh2(h_sd, alphaG=1.5, Z=1.0, phi_g=1.0, sigma_g21=1.9,
        radiation_type='isotropic', return_fractions=True, remove_helium=True,):

    '''
    Calculates ratio of molecular hydrogen to atomic hydrogen surface
    density using the model from Sternberg et al. (2014) given a total hydrogen
    surface density.

    Parameters
    ----------
    h_sd : array-like
        Hydrogen surface density in units of solar mass per parsec**2-
    alphaG : float
        Fundamental dimensionless parameter for the HI-to-H2 transition.
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
        if remove_helium:
            # use the constant 6.71 instead of 9.5 to disclude Helium contribution
            hi_sd = 6.73684 / (Z * phi_g) * np.log(alphaG / 3.2 + 1) # Msun pc^-2
        else:
            hi_sd = 9.5 / (Z * phi_g) * \
                    np.log(alphaG / 3.2 + 1) # Msun pc^-2
    elif radiation_type == 'beamed':
        if remove_helium:
            hi_sd = 11.9 / 1.4101 / \
                    (Z * phi_g) * np.log(alphaG / 2.0 + 1) # Msun pc^-2
        else:
            hi_sd = 11.9 / (Z * phi_g) * np.log(alphaG / 2.0 + 1) # Msun pc^-2
    else:
        raise ValueError('radiation_type must be "beamed" or "isotropic"')

    f_H2 = 1.0 - hi_sd / h_sd
    f_HI = 1.0 - f_H2

    # Keep fractions within real fractional value range
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

def calc_n_H(I_UV=1, alphaG=1, phi_g=1.0, Z_g=1.0, Z_CO=1.0, I_UV_error=0.,
        alphaG_error=0., phi_g_error=0., Z_g_error=0., calc_errors=True):

    ''' Equation 5 in Bialy et al. (2015)
    '''

    #phi_cnm = calc_phi_cnm(alphaG=alphaG, Z_g=Z_g, Z_CO=Z_CO, phi_g=phi_g)

    #n_H = 22.7*I_UV * 4.1 / (1 + 3.1 * Z_g**0.365) * Z_g / Z_CO * phi_cnm / 3.

    n_H = 1.54 * I_UV / alphaG * phi_g / (1 + (2.64 * phi_g * Z_g)**0.5) * 100.

    if 0:
        calc_errors = any(((any(alphaG_error != 0),
                            any(phi_g_error != 0),
                            any(I_UV_error != 0))))

    if calc_errors:

        # Calculated from Mathematica
        n_H_error = \
        np.sqrt((I_UV_error**2 * alphaG**2 * phi_g**2 * \
        (3402.78 + 8983.33 * Z_g * phi_g + 11057.7 * np.sqrt(Z_g * phi_g)) + \
        I_UV**2 * (alphaG_error**2 * phi_g**2 * (3402.78 + 8983.33 * Z_g * phi_g +\
        11057.7 * np.sqrt(Z_g * phi_g)) + alphaG**2 * \
        (3402.78 + 2245.83 * Z_g * phi_g + 5528.86 * np.sqrt(Z_g * phi_g)) * \
        phi_g_error**2)) / (alphaG**4 * (0.615457 + np.sqrt(Z_g * phi_g))**4))


        alphaG_comp = - 154. * I_UV * alphaG_error * phi_g / \
                      (alphaG**2 * (1 + 1.62481 * (Z_g * phi_g)**0.5))

        I_UV_comp =  154. * I_UV_error * phi_g / \
                      (alphaG * (1 + 1.62481 * (Z_g * phi_g)**0.5))


        phi_g_comp = (- 125.11 * I_UV * Z_g * phi_g / \
                        (alphaG * (Z_g * phi_g)**0.5 * \
                            (1 + 1.62481 * (Z_g * phi_g)**0.5)**2) + \
                     154. * I_UV / \
                        (alphaG * (1 + 1.62481 * (Z_g * phi_g)**0.5))) * \
                     phi_g_error

        Z_g_comp = 154. * I_UV * phi_g / alphaG * \
                   (-0.812404 * phi_g / \
                    ((phi_g * Z_g)**0.5 * \
                     (1.62481 * (phi_g * Z_g)**0.5 + 1)**2)) * Z_g_error

        n_H_error_2 = (I_UV_comp**2 + alphaG_comp**2 + phi_g_comp**2 + \
                       Z_g_comp**2)**0.5

        #np.testing.assert_almost_equal(n_H_error, n_H_error_2, decimal=3)
        #np.testing.assert_almost_equal(n_H_error_1, n_H_error_2, decimal=3)

        return n_H, n_H_error

    return n_H

def calc_phi_cnm(alphaG=1.0, Z_g=1.0, Z_CO=1.0, phi_g=1.0, alphaG_error=0.,
        phi_g_error=0.):

    phi_cnm = 2.58 * (1 + 3.1 * Z_g**0.365) / 4.1 * Z_CO / Z_g * \
              3.0 / alphaG * 2.62 / (1 + (2.63 * phi_g * Z_g)**0.5) * phi_g

    return phi_cnm

def calc_w(phi_g=1.0, Z_g=1.0, phi_g_error=(0.0,0.0), Z_g_error=(0.0,0.0)):

    ''' Calculates the normalized H2-dust-limited dissociation bandwidth. alphaG
    can be related to the Krumholz et al. (2009) normalized radiation field
    strength, chi, by alphaG = w * chi. See Equation 93 of Sternberg et al.
    (2014).

    Parameters
    ----------
    phi_g : float
        Order-unity grain composition factor.
    Z_g : float
        Gas-phase metallicity.
    phi_g_error : float
        Uncertainty on phi_g.
    Z_g_error : float
        Uncertainty on Z_g.

    Returns
    -------
    w : float
        Normalized H2-dust-limited dissociation bandwidth
    w_error : float, optional
        Uncertainty on w. Returned if any input parameter errors are greater
        than 0.

    '''

    # calculate bandwidth
    w = 1.0 / (1 + (2.64 * phi_g * Z_g)**0.5)

    # calculate errors if any are greater than 0
    if np.sum((phi_g_error, Z_g_error)) > 0:
        phi_g_error = np.array(phi_g_error)
        Z_g_error = np.array(Z_g_error)

        # https://www.wolframalpha.com/input/?i=partial+derivative+1+%2F+(1+%2B+(2.64+*+phi+*+Z)**0.5)+with+respect+to+phi
        phi_g_comp = -(0.812404*Z_g)/((Z_g*phi_g)**0.5*\
                     (1.62481 *(Z_g * phi_g)**0.5+1)**2)
        phi_g_comp *= phi_g_error

        # https://www.wolframalpha.com/input/?i=partial+derivative+1+%2F+(1+%2B+(2.64+*+phi+*+Z)**0.5)+with+respect+to+Z
        Z_g_comp = -(0.812404 * phi_g)/((Z_g *phi_g)**0.5*\
                   (1.62481*(Z_g*phi_g)**0.5+1)**2)
        Z_g_comp *= Z_g_error

        w_error = (phi_g_comp**2 + Z_g_comp**2)**0.5

        return w, w_error

    return w

def calc_chi(alphaG=1.0, w=1.0/(1+2.64**0.5), alphaG_error=(0.0,0.0),
        w_error=(0.0,0.0)):

    ''' Calculates the normalized radiation field equivalent to the Krumholz
    et al. (2009) parameter chi from Sternberg et al. (2014) parameters. See
    equation 92 of Sternberg et al. (2014).

    Parameters
    ----------
    alphaG : float
        Fundamental dimensionless parameter for the HI-to-H2 transition.
    w : float
        Normalized H2-dust-limited dissociation bandwidth
    alphaG_error : float
        Uncertainty on alphaG.
    w_error : float
        Uncertainty on w.

    Returns
    -------
    chi : float
        Normalized radiation field from Krumholz et al. (2009) = alphaG / w. See
        equation 93 of Sternberg et al. (2014).
    chi_error : float
        Uncertainty on chi. Returned if any input parameter errors are greater
        than 0.

    '''

    # calculate chi
    chi = alphaG / w

    # calculate errors if any are greater than 0
    if np.sum((alphaG_error, w_error)) > 0:
        alphaG_error = np.array(alphaG_error)
        w_error = np.array(w_error)

        alphaG_comp = 1.0 / w
        alphaG_comp *= alphaG_error

        w_comp = - alphaG / w**2
        w_comp *= w_error

        chi_error = (alphaG_comp**2 + w_comp**2)**0.5

        return chi, chi_error

    return chi


