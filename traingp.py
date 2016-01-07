#!/usr/bin/python

import constants
from myscience import funcs
import numpy as np

def create_train_data(velocity_axis=None,
                      velocity_range=None,
                      nchannels=None,
                      width_range=[5, 20],
                      temp_range=[40, 8000],
                      amp_range=None,
                      mean_range=None,
                      ncomponents=4,
                      nspectra=100,
                      rms=0.01):

    '''

    Parameters
    ----------
    velocity_axis : array-like
        Array of velocity value at each channel. Units of km/s.
    velocity_range : array-like
        Used if velocity_axis is None. Lower and upper limit to velocity range.
        Paired with nchannels. Units of km/s.
    nchannels : int, optional
        Used if velocity_axis is None. Number of channels in spectra. Paired
        with velocity_range.
    width_range : array-like
        Lower and upper limit to standard deviation of Gaussians expected in
        spectra. Units of km/s.
    temp_range : array-like
        Lower and upper limit to kinetic temperatures of gas, e.g. CNM and WNM
        temperatures. Temperatures converted to be width_range. Units of K.
    amp_range : array-like
        Lower and upper limit to amplitudes of Gaussian components in K.
    mean_range : array-like
        Lower and upper limit to centers of Gaussian components in km/s.
        Default is (0.25 and 0.75 of (maximum velocity - minimum velocity)) +
        minimum velocity.
    RMS : float
        RMS of spectra. Units of K / (km / s)
    ncomponents : int
        Number of Gaussian components in spectra.
    nspectra : int
        Number of spectra to generate.

    Returns
    -------
    agd_train : dict
        AGD training data set to be used as input for gausspy.gp.train function.

    '''

    # Get FHWM range
    # --------------
    if temp_range is not None:
        FWHM_range = funcs.Tkin_to_FWHM(np.array(temp_range))
    elif width_range is not None:
        FWHM_range = funcs.std_to_FWHM(np.array(width_range))
    else:
        temp_range = [40, 8000]
        FWHM_range = funcs.Tkin_to_FWHM(np.array(temp_range))

    # Check velocity axis
    # -------------------
    if velocity_axis is None and velocity_range is not None:
        velocity_axis = np.linspace(velocity_range[0],
                                    velocity_range[1],
                                    nchannels)
        print 'Created velocity axis'
    else:
        nchannels = len(velocity_axis)

    # Gaussian mean range
    if mean_range is not None:
        #velocity_diff = np.max(velocity_axis) - np.min(velocity_axis)
        #mean_range = np.array((0.25 * velocity_diff, 0.75 * velocity_diff))
        mean_range += np.min(velocity_axis)

    agd_data = format_train_data(RMS=rms,
                                 NCOMPS=ncomponents,
                                 NCHANNELS=nchannels,
                                 NSPECTRA=nspectra,
                                 AMP_lims=amp_range,
                                 FWHM_lims=FWHM_range,
                                 MEAN_lims=mean_range)

    return agd_data

def format_train_data(RMS=0.05, NCOMPS=4, NCHANNELS=512, NSPECTRA=10,
        AMP_lims=None,
        FWHM_lims=None,
        MEAN_lims=None,):

    # Component properties
    AMP_lims = [RMS * 5, RMS * 25]
    FWHM_lims = [10, 35] # channels
    MEAN_lims = [0.25 * NCHANNELS, 0.75 * NCHANNELS]

    # Initialize
    agd_data = {}
    chan = np.arange(NCHANNELS)
    errors = chan * 0. + RMS # Constant noise for all spectra

    # Begin populating data
    for i in range(NSPECTRA):
        spectrum_i = np.random.randn(NCHANNELS) * RMS

        # Sample random components:
        amps = np.random.rand(NCOMPS) * (AMP_lims[1] - AMP_lims[0]) + \
                    AMP_lims[0]
        fwhms = np.random.rand(NCOMPS) * (FWHM_lims[1] - FWHM_lims[0]) + \
                    FWHM_lims[0]
        means = np.random.rand(NCOMPS) * (MEAN_lims[1] - MEAN_lims[0]) + \
                    MEAN_lims[0]

        # Create spectrum
        for a, w, m in zip(amps, fwhms, means):
            spectrum_i += funcs.gaussian(a, w, m)(chan)

        # Enter results into AGD dataset
        agd_data['data_list'] = agd_data.get('data_list', []) + [spectrum_i]
        agd_data['x_values'] = agd_data.get('x_values', []) + [chan]
        agd_data['errors'] = agd_data.get('errors', []) + [errors]

        # If training data, keep answers
        agd_data['amplitudes'] = agd_data.get('amplitudes', []) + [amps]
        agd_data['fwhms'] = agd_data.get('fwhms', []) + [fwhms]
        agd_data['means'] = agd_data.get('means', []) + [means]

    return agd_data







