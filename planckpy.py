#!/usr/bin/python

''' Extracts a region from the Planck data in HEALPix format, and converts to
galactic or equatorial coordinates.

    Original author: Dr. Robert Benjamin - bobbenjamin@me.com

    Edited and ported by: Elijah Bernstein-Cooper - ezbc@astro.wisc.edu

Code hosted at:
git@bitbucket.org:ezbc/planckpy.git

Module requires the following libraries:
    numpy - http://www.scipy.org/scipylib/download.html
    astropy - http://www.astropy.org/
    healpy - https://pypi.python.org/pypi/healpy

'''

import numpy as np
from astropy.io import fits as pf
import healpy

def build_header(header=None, axis_grids=None, reverse_xaxis=True, field=0,
        coord_type='equatorial', xy_res=None, wcs_header=None):

    ''' Builds a header for the extracted Planck image.

    Parameters
    ----------
    header : pyfits.Header instance
        Header from Planck data as a pyfits.Header object instance.
    axis_grids : tuple, array-like
        Tuple of axis grids.
    reverse_xaxis : bool
        x-axis maximum at the origin?
    field : int
        Column to extract from Planck data.
    coord_type : str
        Only option is 'equatorial'. Galactic coordinates will be implemented
        soon.

    Returns
    -------
    header_new : pyfits.Header instance
        Header for extracted region from Planck data as a pyfits.Header object
        instance.

    '''

    # Number of axes
    naxis = len(axis_grids)

    # Basic FITS info
    header_new = pf.Header()
    items = [('SIMPLE', True),
             ('BITPIX', -32),
             ('NAXIS', naxis),
             ]
    for item in items:
        header_new.append(item)

    # length of each axis, x and y correspond to dec and ra lengths
    for i in range(naxis):
        header_new.append(('NAXIS%s' % (i+1), axis_grids[0].shape[-(i+1)]))

    header_new.append(('EXTEND', True))

    column = field + 1

    # CTYPE, CRPIX, CRDELT, CRDELT
    for card in wcs_header.cards:
        header_new.append((card[0], card[1], card[2]))

    # Extended info
    items = [('CELLSCAL', 'CONSTANT'),
             ('BSCALE', 1.0),
             ('BZERO', 0.0),
             ('BLANK', -1),
             ('BUNIT', header['TUNIT%s' % column]),
             ('TELESCOP', 'Planck'),
             ]
    for item in items:
        header_new.append(item)

    # characteristics of data type
    items = [('TYPE', header['TTYPE%s' % column],
                header.comments['TTYPE%s' % column]),
             ('FORM', header['TFORM%s' % column],
                header.comments['TFORM%s' % column]),
            ]
    for item in items:
        header_new.append(item)

    items = [('BMAJ', 5/3600.), # Planck has 5' beam
             ('BMIN', 5/3600.),
             ('EPOCH', 2000.),
             ('EQUINOX', 2000.),
             ]
    for item in items:
        header_new.append(item)

    # Add rest frequency for CO images
    try:
        if header['AST-COMP'] in ['CO-TYPE1', 'CO-TYPE2', 'CO-TYPE3']:
            if field in (0,1,2,3): # then CO J 1-->0
                rest_freq = 115.2712018e9 # Hz
            elif field in (4,5,6,7): # then CO J 2-->1
                rest_freq = 230.5380000e9 # Hz
            elif field in(8,9,10,11): # then CO J 3-->2
                rest_freq = 345.796e9 # Hz
            header_new.append(('RESTFREQ', 115.2712018e9))
    except KeyError:
        pass

    return header_new

def build_wcs(grids=None, coord_reses=None, ranges=None, reverse_xaxis=True,
        coord_type='equatorial'):

    ''' Builds grid of WCS coordinates given x and y pixel coordinate grids.

    Parameters
    ----------
    grids : tuple
    coord_reses : tuple
    ranges : tuple
    reverse_xaxis : bool
    coord_type : str

    Returns
    -------
    x_coords, y_coords : array-like
        RA and dec grids of coordinates.
    wcs_header : astropy.io.fits.header.Header instance
        Header

    '''

    from astropy import wcs

    # Break up tuples to be more easily read
    x_grid, y_grid = grids[0], grids[1]
    x_coord_res, y_coord_res = coord_reses[0], coord_reses[1]
    x_range, y_range = ranges[0], ranges[1]

    if reverse_xaxis:
        x_coord_res *= -1.0

    # Initialize WCS object
    # CRPIX should be at 0, 0 so that unsophisticated software like Kvis can
    # read the coordinates properly.
    w = wcs.WCS(naxis=2)
    w.wcs.cdelt = np.array([x_coord_res, y_coord_res])
    w.wcs.crval = [0, 0]
    if reverse_xaxis:
        w.wcs.crpix = [-x_range[1] / w.wcs.cdelt[0],
                       -y_range[0] / w.wcs.cdelt[1]]
    else:
        w.wcs.crpix = [-x_range[0] / w.wcs.cdelt[0],
                       -y_range[0] / w.wcs.cdelt[1]]

    # Define coordinate type
    if coord_type.lower() == 'equatorial':
        w.wcs.ctype = ['RA---CAR', 'DEC--CAR']
    elif coord_type.lower() == 'galactic':
        w.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']

    # Get pixel coordinates as WCS coordinates
    x_coords, y_coords = w.all_pix2world(x_grid,
                                         y_grid,
                                         1,)

    # Write wcs object as an easily readable header object
    wcs_header = w.to_header()

    return x_coords, y_coords, wcs_header

def get_planck_filename(data_location='./', data_type=None,
        dr_version=1):

    ''' Finds the file name for data_type requested. The files are:
            COM_CompMap_dust-commrul_2048_R1.00.fits
            COM_CompMap_Lfreqfor-commrul_0256_R1.00.fits
            COM_CompMap_Lfreqfor-commrul_2048_R1.00.fits
            HFI_CompMap_CO-Type1_2048_R1.10.fits
            HFI_CompMap_CO-Type2_2048_R1.10.fits
            HFI_CompMap_CO-Type3_2048_R1.10.fits
            HFI_CompMap_DustOpacity_2048_R1.10.fits
            HFI_CompMap_ThermalDustModel_2048_R1.20.fits
            HFI_SkyMap_100_2048_R1.10_nominal.fits
            HFI_SkyMap_143_2048_R1.10_nominal.fits
            HFI_SkyMap_217_2048_R1.10_nominal.fits
            HFI_SkyMap_353_2048_R1.10_nominal.fits
            HFI_SkyMap_545_2048_R1.10_nominal.fits
            HFI_SkyMap_857_2048_R1.10_nominal.fits
            LFI_SkyMap_030_1024_R1.10_nominal.fits
            LFI_SkyMap_044_1024_R1.10_nominal.fits
            LFI_SkyMap_070_1024_R1.10_nominal.fits

    Parameters
    ----------
    data_location : str
        Location of Planck survey archive files.
    data_type : str
        Type of data from the Planck survey. See the following for the full
        descriptions of each data type.
        http://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/

        See the header of the get_data() function to see how to specify the
        data type.
    data_type_version : int
        Data release number.

    Returns
    -------
    data_file : str
        File path to relevant Planck data product.
    nside : int
        Number of sides for HEALPix format.

    '''

    import os

    CO_types = ['CO-Type1', 'CO-Type2', 'CO-Type3']
    freq_types = [ '030', '044', '070', '100', '143', '217', '353', '545',
            '857']
    dust_types = ['Thermal', 'Dust Opacity']
    dust_names = ['dust-commrul', 'ThermalDustModel']

    # Determine which file is chosen
    if data_type in CO_types:
        data_file = 'HFI_CompMap_%s_2048_R1.%s0.fits' % (data_type, dr_version)
        nside = 2048
    elif data_type in dust_types:
        if data_type == 'Dust Opacity':
            data_file = 'HFI_CompMap_%s_2048_R1.%s0.fits' % ('ThermalDustModel',
                    dr_version)
            nside = 2048
        if data_type == 'Thermal':
            data_file = 'COM_CompMap_%s_2048_R1.%s0.fits' % ('dust-commrul',
                    dr_version)
            nside = 2048
    elif data_type in freq_types:
        if data_type in freq_types[:3]:
            receiver_type = 'LFI'
            nside = 1024
        else:
            receiver_type = 'HFI'
            nside = 2048
        data_file = '%s_SkyMap_%s_%s_R1.%s0_nominal.fits' % \
                (receiver_type, data_type, nside, dr_version)
    else:
        raise LookupError('Invalid data type chosen.')

    data_file = data_location + data_file

    # Does file exist on disk?
    if not os.path.isfile(data_file):
        raise IOError('No such file: %s \n Check data release version and \
            the directory including the Planck data.' % data_file)

    return data_file

def gal2sphere(x_coords, y_coords):

    ''' Converts galactic coordinates to angular coordinates of a point on the
    sphere.

    Parameters
    ----------
    x_coords, y_coords : array-like
        N-dimensional galactic x and y coordinates

    Returns
    -------
    phi_grid : array-like
        Rotation angle in spherical coordinates
    theta_grid : array-like
        Azimuthal angle in spherical coordinates

    '''

    # Convert from phi / theta to l/b
    phi_grid = x_coords / 180. * np.pi
    theta_grid = (90. - y_coords) / 180. * np.pi

    return phi_grid, theta_grid

def switch_coords(x_coords, y_coords, coord_type='equatorial'):

    ''' Switches coordinates between equatorial and galactic.

    Parameters
    ----------
    x_coords, y_coords : array-like
        N-dimensional x and y coordinates
    coord_type : str
        Coordinate system of x_coords and y_coords. Options are 'equatorial'
        and 'galactic'. Default is to switch between 'equatorial' to 'galactic'.

    Returns
    -------
    x_coords_sw, y_coords_sw : array-like
        N-dimensional x and y coordinates in the switched coordinate system.

    '''

    from astropy.coordinates import ICRS as eq
    from astropy.coordinates import Galactic as gal
    from astropy import units

    # Convert coordinates to arrays
    x_coords, y_coords = np.copy(x_coords), np.copy(y_coords)

    if coord_type.lower() == 'galactic':
        coords = gal(l=x_coords,
                  b=y_coords,
                  unit=(units.degree, units.degree)
                  )
        x_coords_sw = coords.icrs.ra.deg
        y_coords_sw = coords.icrs.dec.deg
    elif coord_type.lower() == 'equatorial':
        coords = eq(ra=x_coords,
                  dec=y_coords,
                  unit=(units.degree, units.degree)
                  )
        x_coords_sw = coords.galactic.l.deg
        y_coords_sw = coords.galactic.b.deg

    return x_coords_sw, y_coords_sw

def get_data(data_location='./', data_type=None, x_range=(0,360), y_range=(-90,
    90), coord_type='equatorial', field=0, resolution=0.1, cut_last_pixel=False,
    verbose=True, return_header=True, reverse_xaxis=True, dr_version=1):

    ''' Extracts region from Planck data set. Region will be in equatorial
    coordinates. Planck data must be on disk. Visit
    http://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/ to
    download Planck archive data.

    Parameters
    ----------
    data_location : str
        Filepath to location of Planck data. Default is current directory.
    data_type : str
        Data type to choose from. Options are:
            Narrow-band: ['CO-Type1', 'CO-Type2', 'CO-Type3']
                'CO-Type1' fields:
                    0:  12CO J 1-->0 Intensity
                    1:  12CO J 1-->0 Intensity Error
                    2:  12CO J 1-->0 Null Test
                    3:  12CO J 1-->0 Mask
                    4:  12CO J 2-->1 Intensity
                    5:  12CO J 2-->1 Intensity Error
                    6:  12CO J 2-->1 Null Test
                    7:  12CO J 2-->1 Mask
                    8:  12CO J 3-->2 Intensity
                    9:  12CO J 3-->2 Intensity Error
                    10: 12CO J 3-->2 Null Test
                    11: 12CO J 3-->2 Mask
                'CO-Type2' fields:
                    0:  12CO J 1-->0 Intensity
                    1:  12CO J 1-->0 Intensity Error
                    2:  12CO J 1-->0 Null Test
                    3:  12CO J 1-->0 Mask
                    4:  12CO J 2-->1 Intensity
                    5:  12CO J 2-->1 Intensity Error
                    6:  12CO J 2-->1 Null Test
                    7:  12CO J 2-->1 Mask
                'CO-Type3' fields:
                    0:  12CO Intensity
                    1:  12CO Intensity Error
                    2:  12CO Null Test
                    3:  12CO Mask

            Broad-band (GHz): ['030', '044', '070', '100', '143', '217', '353',
                               '545', '857']
                Broad-band fields:
                    0: I stokes
                    1: Hits
                    2: II_cov

            Processed data products: ['Dust Opacity', 'Thermal']
                'Dust Opacity' fields:
                    0: Opacity 353GHz
                    1: Error on opacity
                    2: E(B-V)
                    3: Error on E(B-V)
                    4: T for high freq correction
                    5: Error on T
                    6: Beta for high freq correction
                    7: Error on Beta
                'Thermal' dust model fields:
                    0: Intensity
                    1: Intensity standard deviation
                    2: Intensity ??
                    3: Intensity ??

    x_coord_range : array-like
        Lower and upper longitude. Default is whole sky.
    y_range : array-like
        Lower and upper latitude. Default is whole sky.
    coord_type : str
        Only option is 'equatorial'. Galactic coordinates will be implemented
        soon.
    field : int
        Field in data type.
    resolution : float
        Pixel resolution in arcseconds.
    cut_last_pixel : bool
        Cuts off one pixel
    return_header : bool
        Return the header?
    verbose : bool
        Verbose?
    reverse_xaxis : bool
        The highest x-axis value begins at the origin.
    dr_version : int
        Data release version of data.

    Returns
    -------
    map : array-like
        Map of extracted region from Planck data.
    header : dict, optional
        FITS format header.

    Examples
    --------
    >>> import planckpy as pl
    >>> from astropy.io import fits
    >>> (data, header) = pl.get_data(data_type='857', x_range=(21.3, 30.3),
            y_range=(60.0, 73.0))
    >>> data.shape
    (131, 91)
    >>> header['TYPE']
    'I_STOKES'
    >>> fits.writeto('planck_region_857GHz.fits', data, header=header)

    '''

    if data_type is None:
        print('WARNING (get_data): No data type chosen. Returning None type.')
        return None

    if coord_type is 'galactic':
    	raise ValueError('Galactic coordinates not yet implemented.')

    if np.abs(y_range[0]) > 90.0 or np.abs(y_range[1]) > 90.0:
    	raise ValueError('y coordinates must be < +90 deg and > -90 deg')
    if x_range[0] < 0.0 or \
       x_range[1] < 0.0 or \
       x_range[0] > 360.0 or \
       x_range[1] > 360.0:
    	raise ValueError('x coordinates must be > 0 deg and < 360 deg')

    if x_range[0] >= x_range[1]:
        raise ValueError('x_range[0] must be < x_range[1]')
    if y_range[0] >= y_range[1]:
        raise ValueError('y_range[0] must be < y_range[1]')

    # Get the filename
    filename = get_planck_filename(data_type = data_type, data_location =
            data_location, dr_version = dr_version)

    if verbose:
        print('Reading file:\n%s' % (filename))

    # Read the map using healpy library
    map_data = healpy.read_map(filename, field = field, h = True)
    map_raw, header_raw = map_data[0], map_data[1]

    # Change format of healpy header to pyfits header
    #   ...an actually useful format
    header_pf = pf.Header()
    for item in header_raw:
        header_pf.append(item)

    # Get nside from HEALPix format, acts as resolution of HEALPix image
    nside = header_pf['NSIDE']


    # Set up longitude / latitude grid for extracting the region
    x_coord_res, y_coord_res = resolution, resolution
    x_pixel_count = np.ceil((x_range[1] - x_range[0]) / x_coord_res) + 1
    y_pixel_count = np.ceil((y_range[1] - y_range[0]) / y_coord_res) + 1

    # The counts should be positive! Problems will arise if southern
    # coordinates queried
    x_pixel_count = np.abs(x_pixel_count)
    y_pixel_count = np.abs(y_pixel_count)
    if cut_last_pixel:
        x_pixel_count -= 1
        y_pixel_count -= 1

    # Write pixel positions
    x_pix = np.arange(x_pixel_count)
    y_pix = np.arange(y_pixel_count)

    # Create grids of coordinate positions
    x_grid, y_grid = np.meshgrid(x_pix, y_pix)

    # Need x and y range in equatorial coordinates
    if coord_type.lower() == 'galactic':
        x_range, y_range = switch_coords(x_range, y_range, coord_type)

    # Convert pixel coordinates to WCS coordinates
    x_coords, y_coords, wcs_header = build_wcs(grids=(x_grid, y_grid),
                                        coord_reses=(x_coord_res, y_coord_res),
                                        ranges=(x_range, y_range),
                                        reverse_xaxis=reverse_xaxis,
                                        coord_type='equatorial')

    # Get galactic coordinates for angles
    if coord_type == 'equatorial':
        l_grid, b_grid = switch_coords(x_coords, y_coords, coord_type)
    elif coord_type == 'galactic':
    	l_grid, b_grid = x_coords, y_coords

    # Convert from phi / theta to l/b
    phi_grid, theta_grid = gal2sphere(l_grid, b_grid)

    # Convert from angle to pixel
    pixel_indices = healpy.ang2pix(nside = nside, theta = theta_grid,
            phi = phi_grid, )

    # Map the column data to a 2d array
    map_region = map_raw[pixel_indices]

    # Omit any degenerate axes
    map_region = np.squeeze(map_region)

    # Change bad data to be NaNs
    bad_data = header_pf['BAD_DATA']
    map_region[map_region == bad_data] = np.NaN

    # Rebuild header if galactic coordinates chosen
    if coord_type.lower() == 'galactic':
        x_range, y_range = switch_coords(x_range, y_range, 'equatorial')

        # Convert pixel coordinates to WCS coordinates
        x_coords, y_coords, wcs_header = build_wcs(grids=(x_grid,
                                                          y_grid),
                                                    coord_reses=(x_coord_res,
                                                                 y_coord_res),
                                                    ranges=(x_range,
                                                            y_range),
                                                    reverse_xaxis=reverse_xaxis,
                                                    coord_type=coord_type)

    # Build a header
    if return_header:
        header_region = build_header(header = header_pf,
                wcs_header=wcs_header,
                axis_grids = (x_coords, y_coords),
                reverse_xaxis = reverse_xaxis,
                field = field,
                coord_type = coord_type,
                xy_res = (x_coord_res, y_coord_res))
        return map_region, header_region
    else:
        return map_region


