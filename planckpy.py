#!/usr/bin/python

''' Extracts a region from the Planck data in HEALPix format, and converts to
galactic coordinates.

    Original author: Dr. Robert Benjamin - bobbenjamin@me.com

    Edited and ported by: Elijah Bernstein-Cooper - ezbc@astro.wisc.edu

Code hosted at:
git@bitbucket.org:ezbc/planckpy.git

Module requires the following libraries:
    numpy - http://www.scipy.org/scipylib/download.html
    pyfits - http://www.stsci.edu/institute/software_hardware/pyfits/Download
    healpy - https://pypi.python.org/pypi/healpy


'''

import numpy as np
from astropy.io import fits as pf
import healpy
#from kapteyn import wcs
from astropy.coordinates import ICRS as eq
from astropy.coordinates import Galactic as gal
from astropy import units
from astropy import wcs

def build_header(header = None, axis_grids = None, reverse_xaxis = True, field =
        0, coord_type = 'galactic', xy_res=None):

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
        Options are 'galactic' and 'equatorial'

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
             ('EXTEND', True),
             ]
    for item in items:
        header_new.append(item)

    # length of each axis
    for i in range(naxis):
        header_new.append(('NAXIS%s' % (i+1), len(axis_grids[i])))

    column = field + 1

    # Extended info
    items = [('BSCALE', 1.0),
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

    # Axes characteristics
    if coord_type == 'galactic':
        ctypes = ['GLON--CAR', 'GLAT--CAR', 'VELO-LSR']
    if coord_type == 'equatorial':
        ctypes = ['RA---CAR', 'DEC--CAR', 'VELO-LSR']

    for i, axis in enumerate(axis_grids):
        crpix = 0
    	ctype = ctypes[i]
    	crval = axis[0, 0]

        print 'axis'
        print axis[0:3, 0:3]

        #
        if i == 0: # ra
            # choose cdelt along line of constant dec
            dec_angle = np.cos(np.deg2rad(axis_grids[1][0, 0]))
            cdelt = (axis[1, 0] - axis[0, 0]) * dec_angle

            cdelt = axis_grids[1][0, 1] - axis_grids[1][0, 0]

            cdelt = xy_res[0]

            # if reverse_axis, put the highest RA value at the origin
            if reverse_xaxis:
                cdelt = -cdelt
                crval = axis[-1, 0]
                #crval = 82.5
            elif not reverse_xaxis:
                crval = axis[0, 0]
        elif i == 1: # dec
            cdelt = axis[0, 1] - axis[0, 0]
            cdelt = xy_res[1]

        print ctype, cdelt, crval

    	items = [('CRPIX%s' % (i+1), crpix),
    	         ('CDELT%s' % (i+1), cdelt),
    	         ('CRVAL%s' % (i+1), crval),
    	         ('CTYPE%s' % (i+1), ctype),
    	         ]
        for item in items:
            header_new.append(item)

    items = [#('CELLSCAL','CONSTANT'),
             ('BMAJ', 5/3600.),
             ('BMIN', 5/3600.),
             #('COORDSYS', 'GALACTIC'),
             #('EPOCH', 2000.),
             #('BAD_DATA', -1e28)
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
        nothing = None

    return header_new

def get_planck_filename(data_location = './', data_type = None,
        dr_version = 1):

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

def eq2gal(ra,dec):

    '''

    Original code from:
    http://community.dur.ac.uk/physics.astrolab/python_links.html

    '''

    rmat = numpy.array([
     [-0.054875539726,-0.873437108010,-0.483834985808],
     [+0.494109453312,-0.444829589425,+0.746982251810],
     [-0.867666135858,-0.198076386122,+0.455983795705]],dtype='d')
    cosb = math.cos(dec)
    v1 = numpy.array([math.cos(ra)*cosb,math.sin(ra)*cosb,math.sin(dec)])
    v2 = numpy.multiply(rmat,v1)
    x = v2[0]
    y = v2[1]
    z = v2[2]
    r = math.sqrt(x*x+y*y)
    if r == 0.0:
        l = 0.0
    else:
        l = math.atan2(y,x)
    if z == 0.0:
        b = 0.0
    else:
        b = math.atan2(z,r)
    ll = math.fmod(l,2.*math.pi)
    if ll < 0.0: ll = ll + 2.*math.pi

    bb = math.fmod(b,2*math.pi)
    if abs(bb) >= math.pi: print "Ugh!"

    return ll, bb

def get_data(data_location='./', data_type = None, x_range = (0,360),
        y_range = (-90, 90), coord_type = 'galactic', field = 0,
        resolution = 0.1, cut_last_pixel = False, verbose = True, return_header
        = True, reverse_xaxis = True, dr_version = 1):

    ''' Extracts region from Planck data set. Region will be in galactic
    coordinates. Planck data must be on disk.

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
        Options are 'galactic' and 'equatorial'
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
    >>> import pyfits as pf
    >>> (data, header) = pl.get_data(data_type = '857', x_range =
            (155,165), latitude_range = (-30, -15))
    >>> data.shape
    (151, 101)
    >>> header['TYPE']
    'I_STOKES'
    >>> pf.writeto('planck_region_857GHz.fits', data, header = header)

    '''

    if data_type is None:
        print('WARNING (get_data): No data type chosen. Returning None type.')
    	return None

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

    '''
    # Set up longitude / latitude grid for extracting the region
    dec_angle = np.cos(np.deg2rad(y_range[0]))
    x_coord_res = resolution * dec_angle
    print x_coord_res
    y_coord_res = resolution
    x_pixel_count = np.floor((x_range[1] - x_range[0]) / x_coord_res + 1)
    y_pixel_count = np.floor((y_range[1] - y_range[0]) / y_coord_res + 1)
    if cut_last_pixel:
    	x_pixel_count -= 1
    	y_pixel_count -= 1

    # Write axes of l/b positions
    x_center = (x_range[1] - x_range[0]) / 2.
    y_center = (y_range[1] - y_range[1]) / 2.
    x_axis = - x_center + x_coord_res * np.arange(x_pixel_count)
    #x_axis = x_range[0] + x_coord_res * np.arange(x_pixel_count)
    y_axis = y_range[0] + y_coord_res * np.arange(y_pixel_count)

    # Create grids of coordinate positions
    # Using the plate-caree projection
    #   constant arclengths of the image will be sampled
    #   right ascension grid will span a larger range in hours with
    #   increasing declination
    x_grid = np.zeros(shape = (x_pixel_count, y_pixel_count))
    y_grid = np.zeros(shape = (x_pixel_count, y_pixel_count))
    for y in range(len(y_axis)):
    	dec_angle = np.cos(np.deg2rad(y_axis[y]))
        x_grid[:, y] = x_range[0] + x_center + x_axis / dec_angle
        #x_grid[:, y] = x_axis / dec_angle
    for x in range(len(x_axis)):
        y_grid[x, :] = y_axis

    print('DEC')
    print(y_grid[0, -1])
    print('RA')
    print(x_grid[:, -1])
    print('DEC')
    print(y_grid[0, 0])
    print('RA')
    print(x_grid[:, 0])
    '''


    # Set up longitude / latitude grid for extracting the region
    x_coord_res, y_coord_res = resolution, resolution
    x_pixel_count = np.ceil((x_range[1] - x_range[0]) / x_coord_res) + 1
    y_pixel_count = np.ceil((y_range[1] - y_range[0]) / y_coord_res) + 1
    if cut_last_pixel:
    	x_pixel_count -= 1
    	y_pixel_count -= 1

    # Write pixel positions
    x_pix = np.arange(x_pixel_count)
    y_pix = np.arange(y_pixel_count)

    # Create grids of coordinate positions
    x_grid, y_grid = np.meshgrid(x_pix, y_pix)

    w = wcs.WCS(naxis=2)
    w.wcs.cdelt = np.array([-x_coord_res, y_coord_res])
    w.wcs.crpix = [0, 0]
    w.wcs.crval = [x_range[1], y_range[0]]
    w.wcs.ctype = ['RA---CAR', 'DEC--CAR']

    x_coords, y_coords = w.all_pix2world(x_grid,
                                         y_grid,
                                         1,)
                                         #ra_dec_order=True)

    print x_coords[0:3, 0:3]
    print x_coords[-4:-1, -4:-1]
    print y_coords[0:3, 0:3]

    # Convert to galactic coordinates
    if coord_type.lower() == 'equatorial':
        axes = eq(ra=x_coords,
                  dec=y_coords,
                  unit=(units.degree, units.degree)
                  )
        longitude_grid = axes.galactic.l.deg
        latitude_grid = axes.galactic.b.deg
    elif coord_type.lower() == 'galactic':
        longitude_grid = x_coords
        latitude_grid = y_coords

    # Convert from phi / theta to l/b
    phi_grid = longitude_grid / 180. * np.pi
    theta_grid = (90. - latitude_grid) / 180. * np.pi

    # Convert from angle to pixel
    pixel_indices = healpy.ang2pix(nside = nside, theta = theta_grid,
            phi = phi_grid, )

    # Map the column data to a 2d array
    map_region = map_raw[pixel_indices]

    # Omit any degenerate axes
    map_region = np.squeeze(map_region)

    # Move set bad_data to be >-1e30 to be read by kvis
    map_region[map_region < -1e30] = np.NaN

    # Reverse the array
    if reverse_xaxis:
        #map_region = map_region.T[::, ::-1]
        #map_region = map_region[::, ::-1]
        map_region = map_region
    elif not reverse_xaxis:
        #map_region = map_region.T
        map_region = map_region

    # Build a header
    if return_header:
    	header_region = build_header(header = header_pf,
    	        axis_grids = (x_coords, y_coords),
    	        reverse_xaxis = reverse_xaxis,
    	        field = field,
    	        coord_type = coord_type,
    	        xy_res = (x_coord_res, y_coord_res))
        return map_region, header_region
    else:
        return map_region


