#!/usr/bin/python

import numpy as np

def calc_sf_threshold(alphaG, sigma_g=1.0):

    ''' Calculates star formation threshold gas surface density from Bialy et
    al. (2016).

    Parameters
    ----------
    alphaG : float
        alphaG parameter from Sternberg et al. (2014). Dimensionless.
    sigma_g : float, optional
        Dust absorption cross section. Units relative to galactic value of 1.9 x
        10^-21 cm^-2.

    Returns
    -------
    h_sd : float
        Gas surface density threshold. Units of solar masses pc^-2.

    '''

    h_sd = 15.8 / sigma_g * np.log((alphaG / 2.0)**1.43 + 1)

    return h_sd


