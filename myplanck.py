#!/usr/bin/python

''' Extracts a region from the Planck data in HEALpix format, and converts to
cartesian or similar projection.

Original Author: Dr. Bob Benjamin.


'''

import numpy as np
import pyfits as pf
import healpy

def build_header(header = None, axes = None, reverse_xaxis = True,):

    # Change format of healpy header to pyfits header
    #   ...an actually useful format
    header_edit = pf.Header()
    for item in header:
        header_edit.append(item)

    # Number of axes
    naxis = len(axes)

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

    items = [('BSCALE', 1.0),
             ('BZERO', 0.0),
             ('BLANK', -1),
             ('BUNIT', header_edit['TUNIT1']),
             ('TELESCOP', 'Planck'),
             ]
    for item in items:
        header_new.append(item)

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

        print('crval = %f' % crval)

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
             ('EPOCH', 2000.)]
    for item in items:
        header_new.append(item)

    # add rest frequency
    if header_edit['AST-COMP'] in ['CO-TYPE1', 'CO-TYPE2', 'CO-TYPE3']:
    	header_new.append(('RESTFREQ', 115.2712018e9))

    return header_new

def get_planck_filename(data_location = './', data_type = None,
        data_type_version = 10):

    ''' Finds the file name for data_type requested. The files are:
            COM_CompMap_dust-commrul_2048_R1.00.fits
            COM_CompMap_Lfreqfor-commrul_0256_R1.00.fits
            COM_CompMap_Lfreqfor-commrul_2048_R1.00.fits
            HFI_CompMap_CO-Type1_2048_R1.10.fits
            HFI_CompMap_CO-Type2_2048_R1.10.fits
            HFI_CompMap_CO-Type3_2048_R1.10.fits
            HFI_CompMap_DustOpacity_2048_R1.10.fits
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

        The options are:
            CO:
                'CO-Type1'
                'CO-Type2'
                'CO-Type3'
            Frequency bands (in GHz):
                '030'
                '044'
                '070'
                '100'
                '143'
                '217'
                '353'
                '545',
                '857'
            Processed data products:
                'Thermal'
                'Dust Opacity'
                'TAU353',
                'TAU353ERR',
                'EBV',
                'EBV_ERR',
                'T_HF',
                'T_HF_ERR',
                'BETAHF',
                'BETAHFERR'
    data_type_version : int
        Data release number in the format of a two digit integer. Data release
        one would be 10.

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
    dust_names = ['dust-commrul', 'DustOpacity']

    ttypeArray = ['TAU353',
                  'TAU353ERR',
                  'EBV',
                  'EBV_ERR',
                  'T_HF',
                  'T_HF_ERR',
                  'BETAHF',
                  'BETAHFERR']
    ttypeDescArray = ['Opacity 353GHz',
                      'Error on opacity',
                      'E(B-V)',
                      'Error on E(B-V)',
                      'T for high freq correction',
                      'Error on T',
                      'Beta for high freq correction',
                      'Error on Beta']

    # Determine which file is chosen
    if data_type in CO_types:
        data_file = 'HFI_CompMap_%s_2048_R1.10.fits' % data_type
        nside = 2048
    elif data_type in dust_types:
        data_name = dust_names[dust_types == data_type]
        if data_type == 'Dust Opacity':
            data_file = 'HFI_CompMap_%s_2048_R1.10.fits' % data_name
            nside = 2048
        elif data_type == 'Thermal':
            data_file = 'COM_CompMap_%s_2048_R1.10.fits' % data_name
            nside = 2048
    elif data_type in freq_types:
        if data_type in freq_types[:3]:
        	receiver_type = 'LFI'
        	nside = 1024
        else:
        	receiver_type = 'HFI'
        	nside = 2048
        	data_file = '%s_SkyMap_%s_%s_R1.10_nominal.fits' % \
                    (receiver_type, data_type, nside,)

    data_file = data_location + data_file

    if not os.path.isfile(data_file):
        raise IOError('No such file: %s' % data_file)

    return data_file, nside

def get_data(data_location='./', data_type = 'CO-Type3', longitude_range =
        (152,180), latitude_range = (-25,-3), field = 0, resolution = 0.01,
        cut_last_pixel = False, verbose = True, return_header = True,
        reverse_xaxis = True):

    ''' Extracts region from Planck data set.

    Parameters
    ----------
    data_location : str
        Filepath to location of Planck data. Default is current directory.
    data_type : str
        Data type to choose from. Options are:
            Narrow-band: 'CO-Type1', 'CO-Type2', 'CO-Type3'
            Broad-band (GHz): '030', '044', '070', '100', '143', '217', '353',
                        '545', '857'
            Processed data products: 'Thermal', 'Dust Opacity', 'TAU353',
                'TAU353ERR', 'EBV', 'EBV_ERR', 'T_HF', 'T_HF_ERR', 'BETAHF',
                'BETAHFERR'
    longitude : array-like
        Lower and upper longitude.
    latitude : array-like
        Lower and upper latitude.
    field : int
        Field in data.
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

    Returns
    -------
    map : array-like
        Map of extracted region from Planck data.
    header : dict, optional
        FITS format header.

    '''

    # Get the filename
    filename, nside = get_planck_filename(data_type = data_type, data_location =
            data_location)

    if verbose:
        print('Reading file:\n%s' % (filename))

    # Read the map using healpy library
    map_data = healpy.read_map(filename, field = field, h = True)
    map_raw, header_raw = map_data[0], map_data[1]

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

    if return_header:
    	header_region = build_header(header = header_raw, axes =
    	        (longitude_axis, latitude_axis), reverse_xaxis = reverse_xaxis)
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

    # Extract the data
    (data, header) = get_data(data_location=data_location, data_type='CO-Type3',
            longitude_range = (145,165), latitude_range=(-30,-5), field=0,
            resolution = 0.01, cut_last_pixel = False, verbose = True)

    # Write the data to FITS format
    pf.writeto('/d/bip3/ezbc/perseus/data/planck/perseus_planck_region.fits',
            data, header = header, clobber = True, output_verify = 'fix')


if __name__ == '__main__':
	main()


