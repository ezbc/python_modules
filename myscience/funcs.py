#!/usr/bin/python

import numpy as np

def calc_radiation_field(T_dust, a=100, T_dust_error=0., beta=2.0,
        beta_error=0.0, field_type='draine'):

    ''' Calculates ambient interstellar radiation field from dust temperature.
    Dust is assumed to be in thermal equilibrium and absorbs all incident
    radiation. See Equation 5 in Lee et al. (2012).

    Parameters
    ----------
    T_dust : float
        Average dust temperature in K.
    a : float
        Average dust grain radius in 100 nanometers. Typical radius = 100 nm
    T_dust_error : float, optional
        Average dust temperature uncertainty in K.
    beta : float, optional
        Slope of modified black body spectrum.
    beta_error : float, optional
        Uncertainty of slope of modified black body spectrum.
    field_type : str, optional
        Radiation field type. Options are 'draine', 'mathis', and 'habing'.

    Returns
    -------
    I_UV : float
        Ambient FUV radiation field relative to Draine Field.

    Notes
    -----
    Min's Description:

    In Lee+12, we used a modified black-body (emissivity spectral index beta of
    2) + IRAS 60/100 to derive the dust temperature.

    In this case, the ISRF heating dust grains (from *all directions*) can be
    actually easily calculated by

    U(M83) = (T_dust / T_norm K)^6

    --- This relation is based on many observations showing that large
    interstellar silicates have an equilibrium temperature of T_norm K in the
    solar neighborhood.

    --- Here U(M83) is the ISRF by Mathis+83 integrated from 0.09 micron to 8
    micron.

    To be specific, U(M83) = 4pi * \int J(nu) d(nu) (0.09 micron ~ 8 micron)
    and U(M83) = 1 corresponds to 2.2e-5 W m-2.

    --- For Perseus with T_dust ~ 17 K, U(M83) ~ 0.8.

    --- But what we need to know is I_UV, the scaling factor for the Draine
    field.

    For the FUV range of 6 ~ 13.6 eV, the Draine field is a factor of ~1.5
    stronger than the Mathis field.  Based on all these estimates, I would
    expect I_UV ~ (1/1.5) * 0.8 ~ 0.5 for Perseus.

    --- Another thing I want to mention is that the "empirical" relation I used
    above is independent from the dust (LW photon) absorption cross section.

    For more information visit:
    http://ezbc.me/research/2016/02/25/radiation-field/

    '''

    # Calculate Mathis field for dust in equilibrium
    T_norm = 17.5 # K
    U = (T_dust / T_norm)**(4.0 + beta)

    if np.any(T_dust_error > 0.0) or np.any(beta_error > 0.0):

        # derive for T_dust
        # http://www.wolframalpha.com/input/?i=d%2FdT+(T%2Fc)%5E(4%2BB)
        T_dust_comp = (beta + 4.0) * (T_dust / T_norm)**(beta + 3.0) * \
                      T_dust_error

        # derive for beta
        # http://www.wolframalpha.com/input/?i=d%2FdB+(T%2Fc)%5E(4%2BB)
        beta_comp = (T_dust / T_norm)**(beta + 4.0) * np.log(T_dust / T_norm) \
                    * beta_error

        U_error = (T_dust_comp**2 + beta_comp**2)**0.5

        if type(U_error) is np.float64 or type(U_error) is float:
            U_error = np.array((U_error, U_error))

        # Scale Draine Field by the Draine / Mathis field in the FUV
        if field_type == 'draine':
            U *= 1.48
            U_error *= 1.48
        elif field_type == 'habing':
            U *= 1.14
            U_error *= 1.14

        return U, U_error

    # Scale Draine Field by the Draine / Mathis field in the FUV
    if field_type == 'draine':
        U *= 1.48
    elif field_type == 'habing':
        U *= 1.14

    return U

def calc_temperature(n_H=1.0, pressure=3800.0, pressure_error=(100,100),
        n_H_error=0, calc_error=True):

    ''' Calculates temperature of atomic hydrogen assuming thermal equilibrium.
    P/k = n_H * T --> T = (P / k) / n_H

    Parameters
    ----------
    n_H : float
        Atomic hydrogen number density.
    pressure : float
        P / k of atomic hydrogen in K / cm^3

    Returns
    -------
    T : float
        Temperature of atomic hydrogen in K.

    '''

    T = pressure / n_H

    if calc_error:
        n_H_error = np.array(n_H_error, dtype=float)
        pressure_error = np.array(pressure_error, dtype=float)
        n_H_comp = pressure / n_H**2 * n_H_error
        pressure_comp =  1.0 / n_H * pressure_error
        T_error = (n_H_comp**2 + pressure_comp**2)**0.5
        return T, T_error

    return T

def calc_density(T_H=1.0, pressure=3800.0, pressure_error=(100,100),
        T_H_error=0, calc_error=True):

    ''' Calculates pressure of atomic hydrogen assuming thermal equilibrium.
    P/k = T_H * n --> n = (P / k) / T_H

    Parameters
    ----------
    T_H : float
        Atomic hydrogen number density.
    pressure : float
        P / k of atomic hydrogen in K / cm^3

    Returns
    -------
    n : float
        Temperature of atomic hydrogen in K.
    '''

    n = pressure / T_H

    if calc_error:
        T_H_error = np.array(T_H_error, dtype=float)
        pressure_error = np.array(pressure_error, dtype=float)
        T_H_comp = pressure / T_H**2 * T_H_error
        pressure_comp =  1.0 / T_H * pressure_error
        n_error = (T_H_comp**2 + pressure_comp**2)**0.5
        return n, n_error

    return n

def Tkin_to_FWHM(Tkin):

    '''
    Tkin <= m_H * FWHM^2 / (8 k_B * ln(2))
    FWHM = (Tkin * (8 k_B * ln(2)) / m_H)^2
    '''

    import constants

    FWHM = (Tkin * 8 * constants.cgs.k * np.log(2) / constants.cgs.mh)**0.5

    return FWHM

def std_to_FWHM(std):

    return 2.355 * std

def FWHM_to_std(FWHM):

    return FWHM / 2.355

def gaussian(amp, fwhm, mean):
    return lambda x: amp * np.exp(-(x-mean)**2/4./(fwhm/2.355)**2)

