#!/usr/bin/python

def check_residuals(xarr,profile,threshold):
    """ Checks tile at to find any significant peaks in the
    residuals that should be fitted. Threshold is in units of standard
    deviations.
    """

    from numpy import where

    profileMax = profile.max()

    if profileMax is not None:
        if profileMax >= threshold:
            return (profileMax, xarr[where(profile == profileMax)[0][0]])
        else:
            return None
    else:
        return None

def make_2gaussGuesses(xarr,profile,threshold):
    """ Checks tile at to find any significant peaks in the
    residuals that should be fitted. Threshold is in units of standard
    deviations.
    """

    from numpy import where

    profileMax = profile.max()

    guesses = []

    if profileMax >= threshold:
        

        return (profileMax, xarr[where(profile == profileMax)[0][0]])

def guess_width(xarr,profile,peakPos):
    ''' Guesses the width of a residual peak.
    peak : int
        In xarr units.
    '''
    from numpy import where

    # Find amplitude of maximum residual
    peakVal = profile[peakPos]
    # Define width as maxpos minus position where residual is 0.5 maxpos
    try:
        minPoses = where(profile[peakPos:-1] < 0.5 * peakVal)[0]
        widthPos = abs(min(minPoses) - peakPos)
        width = abs(xarr[widthPos] - xarr[peakPos])
    except ValueError:
        minPoses = where(profile[0:peakPos] < 0.5 * peakVal)[0]
        widthPos = abs(max(minPoses) - peakPos)
        width = abs( -xarr[widthPos] + xarr[peakPos])
    return width[0]

def make_velocityAxis(h):
    """ Creates the velocity axis given a header.
    """
    from numpy import arange

    array = (arange(h['NAXIS3']) - h['CRPIX3'] + 1) * h['CDELT3'] + h['CRVAL3']

    return array / 1000.

def modelData(x,amp,offset,width):

    profile = np.zeros(len(x))

    for i in xrange(len(x)):
        profile[i] = amp*np.exp(-(x[i]-offset)**2/(2*width**2))

    return profile









