#!/usr/bin/python

import numpy as np

def calc_radiation_field(T_dust, a=100):

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
    U : float
        Ambient radiation field in Draine Field.

    '''

    G = 4.6e-11 * a * T_dust**6 # erg cm^-2 s^-1
    D_0 = 5.8e-11 # I_UV s^-1
    G_0 = G / D_0 / 2.07e7

    if T_dust > 0:
        G_0 = (15 / (T_dust * (a / 100)**-(1/15)))**6
    else:
        return np.nan

    return G_0

def calc_temperature(n_H=1.0, pressure=3800.0):

    ''' P/k = n_H * T --> T = (P / k) / n_H

    '''

    T = pressure / n_H

    return T


