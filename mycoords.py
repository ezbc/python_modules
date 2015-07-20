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


def convert_limit_coordinates(prop_dict,
        coords=('region_limit', 'co_noise_limits', 'plot_limit'), header=None):

    # Initialize pixel keys
    for coord in coords:
        prop_dict[coord].update({'pixel': []})

        if coord in ('region_limit',
                     'plot_limit',
                     'region_limit_bin',
                     'plot_limit_bin'):
            limit_wcs = prop_dict[coord]['wcs']

            for limits in limit_wcs:
                # convert centers to pixel coords
                limit_pixels = get_pix_coords(ra=limits[0],
                                             dec=limits[1],
                                             header=header)[:2].tolist()

                prop_dict[coord]['pixel'].append(limit_pixels[0])
                prop_dict[coord]['pixel'].append(limit_pixels[1])
        elif coord == 'co_noise_limits':
            region_limits = prop_dict[coord]['wcs']

            # Cycle through each region, convert WCS limits to pixels
            for region in region_limits:
                region_pixels = []
                for limits in region:
                    # convert centers to pixel coords
                    limit_pixels = get_pix_coords(ra=limits[0],
                                                  dec=limits[1],
                                                  header=header)[:2].tolist()
                    region_pixels.append(limit_pixels)

                # Append individual regions back to CO noise
                prop_dict[coord]['pixel'].append(region_pixels)

    return prop_dict

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec), or Ra in degrees
    and dec in degrees.
    '''

    import pywcsgrid2 as wcs
    import pywcs

    # convert to degrees if ra and dec are array-like
    try:
        if len(ra) == 3 and len(dec) == 3:
            ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
        else:
            raise ValueError('RA and Dec must be in (hrs,min,sec) and' + \
                    ' (deg,arcmin,arcsec) or in degrees.')
    except TypeError:
        ra_deg, dec_deg = ra, dec

    wcs_header = pywcs.WCS(header)
    pix_coords = wcs_header.wcs_sky2pix([[ra_deg, dec_deg, 0]], 0)[0]

    return pix_coords

def hrs2degs(ra=None, dec=None):
    ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
    '''

    ra_deg = 15*(ra[0] + ra[1]/60. + ra[2]/3600.)
    dec_deg = dec[0] + dec[1]/60. + dec[2]/3600.

    return (ra_deg, dec_deg)

def load_ds9_region(props, filename=None, header=None, key='regions'):

    import pyregion as pyr

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    props[key] = {}

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        if region.comment is not None:
            tag = region.comment
            region_name = tag[tag.find('text={')+6:tag.find('}')].lower()

            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            poly_verts_pix = []
            for i in xrange(0, len(poly_verts)):
                poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                                dec=poly_verts[i][1],
                                                header=header)[:-1][::-1].tolist())

            props[key][region_name] = {}
            props[key][region_name]['poly_verts'] = {}
            props[key][region_name]['poly_verts']['wcs'] = poly_verts
            props[key][region_name]['poly_verts']['pixel'] = poly_verts_pix

    return props






