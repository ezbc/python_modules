#!/usr/bin/python

''' Extracts a region from the Planck data in HEALPix format, and converts to
cartesian or similar projection.

Original Author: Dr. Bob Benjamin.


'''

import numpy as np
import pyfits as pf
import healpy

def build_header(header = None, axes = None, reverse_xaxis = True, field =
        0):

    ''' Builds a header for the extracted Planck image.

    Parameters
    ----------
    header : pyfits.Header instance
        Header from Planck data as a pyfits.Header object instance.
    axes : tuple, array-like
        Tuple of axis grids.
    reverse_xaxis : bool
        x-axis maximum at the origin?
    field : int
        Column to extract from Planck data.

    Returns
    -------
    header_new : pyfits.Header instance
        Header for extracted region from Planck data as a pyfits.Header object
        instance.

    '''

    # Number of axes
    naxis = len(axes)

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
        header_new.append(('NAXIS%s' % (i+1), len(axes[i])))

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
    ctypes = ['GLON', 'GLAT', 'VELO-LSR']
    for i, axis in enumerate(axes):
        crpix = 0
    	cdelt = axis[1] - axis[0]
    	ctype = ctypes[i]
    	crval = axis[0]
        if i == 0 and reverse_xaxis:
        	cdelt = -cdelt
        	crval = axis[-1]
        elif i == 0 and not reverse_xaxis:
    	    crval = axis[0]

    	items = [('CRPIX%s' % (i+1), crpix),
    	         ('CDELT%s' % (i+1), cdelt),
    	         ('CRVAL%s' % (i+1), crval),
    	         ('CTYPE%s' % (i+1), ctype),
    	         ]
        for item in items:
            header_new.append(item)

    items = [('CELLSCAL','CONSTANT'),
             ('BMAJ', 5/3600.),
             ('BMIN', 5/3600.),
             ('COORDSYS', 'GALACTIC'),
             ('EPOCH', 2000.),
             ('BAD_DATA', header['BAD_DATA'])]
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

def get_data(data_location='./', data_type = None, longitude_range =
        (0,360), latitude_range = (-90, 90), field = 0, resolution = 0.1,
        cut_last_pixel = False, verbose = True, return_header = True,
        reverse_xaxis = True, dr_version = 1):

    ''' Extracts region from Planck data set. Region will be in galactic
    coordinates.

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

    longitude : array-like
        Lower and upper longitude. Default is whole sky.
    latitude : array-like
        Lower and upper latitude. Default is whole sky.
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
    >>> (data, header) = pl.get_data(data_type = '857', longitude_range =
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

    # Set up longitude / latitude grid for extracting the region
    longitude_res, latitude_res = resolution, resolution
    pixel_count_long = (longitude_range[1] - longitude_range[0]) / \
        longitude_res + 1
    pixel_count_lat = (latitude_range[1] - latitude_range[0]) / \
        latitude_res + 1
    if cut_last_pixel:
    	pixel_count_long -= 1
    	pixel_count_lat -= 1

    # Write axes of l/b positions
    longitude_axis = longitude_range[0] + longitude_res * \
            np.arange(pixel_count_long)
    latitude_axis = latitude_range[0] + latitude_res * \
            np.arange(pixel_count_lat)

    # Create map of l/b positions
    longitude_grid = np.zeros(shape = (pixel_count_long, pixel_count_lat))
    latitude_grid = np.zeros(shape = (pixel_count_long, pixel_count_lat))
    for b in range(len(latitude_axis)):
        longitude_grid[:, b] = longitude_axis
    for l in range(len(longitude_axis)):
        latitude_grid[l, :] = latitude_axis

    # Convert from phi / theta to l/b
    phi_grid = longitude_grid / 180. * np.pi
    theta_grid = (90. - latitude_grid) / 180. * np.pi

    # Convert from angle to pixel
    pixel_indices = healpy.ang2pix(nside = nside, theta = theta_grid,
            phi = phi_grid, )

    # Map the column data to a 2d array
    map_region = map_raw[pixel_indices]

    # Omit the degenerate axes
    map_region = np.squeeze(map_region)

    # Reverse the array
    if reverse_xaxis:
        map_region = map_region.T[::, ::-1]
    elif not reverse_xaxis:
        map_region = map_region.T

    # Build a header
    if return_header:
    	header_region = build_header(header = header_pf, axes =
    	        (longitude_axis, latitude_axis), reverse_xaxis = reverse_xaxis,
    	        field = field)
        return map_region, header_region
    else:
        return map_region

def main():
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    if False:
        (data, header) = get_data(data_location=data_location,
                data_type='CO-Type3', longitude_range=(152,180),
                latitude_range=(-25,-3), field=0, resolution=0.01,
                cut_last_pixel=False, verbose = True)

        pf.writeto('/d/bip3/ezbc/taurus/data/planck/taurus_planck_region.fits',
                data, header = header, clobber = True, output_verify = 'fix')

    if 0:
        # Extract the data
        (data, header) = get_data(data_location=data_location,
                data_type='CO-Type3', longitude_range = (145,165),
                latitude_range=(-30,-5), field=0, resolution = 0.01,
                cut_last_pixel = False, verbose = True)

        print(header)

        # Write the data to FITS format
        pf.writeto('/d/bip3/ezbc/perseus/data/planck/perseus_planck_region.fits',
                data, header = header, clobber = True, output_verify = 'fix')

    if 1:
        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = get_data(data_location=data_location,
                data_type='Dust Opacity', longitude_range = (145,150),
                latitude_range=(-10,-5), field=0, resolution = 0.1,
                cut_last_pixel = False, verbose = True, dr_version = 2)

        print(header)

        # Write the data to FITS format
        pf.writeto(output_dir + 'tau353.fits', data, header = header, clobber
                = True, output_verify = 'fix')

if __name__ == '__main__':
	main()


