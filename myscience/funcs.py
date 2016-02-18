#!/usr/bin/python

import numpy as np

def calc_radiation_field(T_dust, a=100, T_dust_error=0.):

    ''' Calculates ambient interstellar radiation field from dust temperature.
    Dust is assumed to be in thermal equilibrium and absorbs all incident
    radiation. See Equation 5 in Lee et al. (2012).

    Parameters
    ----------
    T_dust : float
        Average dust temperature in K.
    a : float
        Average dust grain radius in 100 nanometers. Typical radius = 100 nm

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

    U(M83) = (T_dust / 17.5 K)^6

    --- This relation is based on many observations showing that large
    interstellar silicates have an equilibrium temperature of 17.5 K in the
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

    '''

    # Calculate Draine field for dust in equilibrium
    U = (T_dust / 17.5)**6.0

    # Scale Draine Field by the Draine / Mathis field in the FUV
    I_UV = U * 1.5

    if T_dust_error > 0:
        U_error = np.sqrt((1/6.0 * (T_dust)**5 / 17.5**6 * T_dust_error)**0.5)
        I_UV_error = U_error * 1.5

        if type(I_UV_error) is np.float64 or type(I_UV_error) is float:
            I_UV_error = np.array((I_UV_error, I_UV_error))

        return I_UV, I_UV_error

    return I_UV

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

def calc_density(T_H=1.0, pressure=3700.0, pressure_error=(1200,1200),
        T_H_error=0, calc_error=True):

    ''' Calculates density of atomic hydrogen assuming thermal equilibrium.
    P/k = n_H * T --> n_H = (P / k) / T_H. Assumes pressure from Jenkins & Tripp
    (2011).

    Parameters
    ----------
    T_H : float
        Atomic hydrogen kinetic temperature.
    pressure : float
        P / k of atomic hydrogen in K / cm^3
    T_H_error : float, array-like
        Atomic hydrogen kinetic temperature error.

    Returns
    -------
    n_H : float
        Density of atomic hydrogen in cm^-3.
    n_H_error : float, optional
        Density error of atomic hydrogen in cm^-3.

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

