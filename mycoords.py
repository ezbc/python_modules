#!/usr/bin/python

''' Module for basic coordinate creation.
'''

def make_velocity_axis(h):
    """ Creates the velocity axis given a pyfits header. Assumes the third
    axis is the velocity axis in km/s using the radio definition.
    """

    from numpy import arange

    array = (arange(h['NAXIS3']) - h['CRPIX3'] + 1) * h['CDELT3'] + h['CRVAL3']

    return array / 1000.

def calc_image_origin(x_limits=None, y_limits=None, delta_x=0.01, delta_y=0.01,
        coord_type='equatorial', ref_wcs=(0.0, 0.0)):

    ''' Calculates pixel coordinates of an image at x-coord, y-coord = (0.0,
    0.0) degrees in either equatorial or galactic coordinates.

    Parameters
    ----------
    x_limits : tuple
        Upper and lower limits of x-coord in degrees.
    y_limits : tuple
        Upper and lower limits of y-coord in degrees.'
    delta_x : float
        Arclength of x-coord in degrees.
    delta_y : float
        Arclength of y-coord in degrees.
    coord_type : str, optional
        'equatorial' or 'galactic'
    ref_wcs : tuple, optional
        Reference x and y coordinates in degrees. Default is (0.0, 0.0), which
        is needed to be read by unsophisticated programs like Kvis.

    Returns
    -------
    ref_pixel : tuple
        Reference right ascension and declination in pixels for reference
        coordinates.
    npix : tuple
        Number of pixels along right ascension and declinations axes.

    '''

    from astropy import wcs

    # Initialize WCS object
    # CRPIX should be at 0, 0 so that unsophisticated software like Kvis can
    # read the coordinates properly.
    w = wcs.WCS(naxis=2, relax=True)
    w.wcs.cdelt = [delta_x, delta_y]
    w.wcs.crval = [x_limits[0], y_limits[0]]
    w.wcs.crpix = [0, 0]

    # Define coordinate type
    #if coord_type.lower() == 'equatorial':
    #    w.wcs.ctype = ['RA---CAR', 'DEC--CAR']
    #elif coord_type.lower() == 'galactic':
    #    w.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']

    ref_pixel = w.wcs_world2pix((ref_wcs,), 1)[0]

    # Get pixel dimensions of image
    x_pix = w.wcs_world2pix(((x_limits[0], y_limits[0]),
                             (x_limits[1], y_limits[0])),
                            1)
    y_pix = w.wcs_world2pix(((x_limits[0], y_limits[0]),
                             (x_limits[0], y_limits[1])),
                            1)

    x_length_pix = x_pix[1, 0] - x_pix[0, 0]
    y_length_pix = y_pix[1, 1] - y_pix[1, 0]

    npix = (x_length_pix, y_length_pix)

    return ref_pixel, npix

    '''
        dec_range = (21.3, 30.3)
        ra_range = (60.0, 73.0)


Ordering converted to RING
WCSAXES =                    2 / Number of coordinate axes                      CRPIX1  =               7300.0 / Pixel coordinate of reference point            CRPIX2  =              -2130.0 / Pixel coordinate of reference point            CDELT1  =                -0.01 / [deg] Coordinate increment at reference point  CDELT2  =                 0.01 / [deg] Coordinate increment at reference point  CUNIT1  = 'deg'                / Units of coordinate increment and value        CUNIT2  = 'deg'                / Units of coordinate increment and value        CTYPE1  = 'RA---CAR'           / Right ascension, plate caree projection        CTYPE2  = 'DEC--CAR'           / Declination, plate caree projection            CRVAL1  =                  0.0 / [deg] Coordinate value at reference point      CRVAL2  =                  0.0 / [deg] Coordinate value at reference point      LONPOLE =                  0.0 / [deg] Native longitude of celestial pole       LATPOLE =                 90.0 / [deg] Native latitude of celestial pole        END
    '''





