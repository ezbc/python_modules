#!/usr/bin/python

''' Module for image analysis in astro applications.

'''

import numpy as np
import mymath
import mycoords
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import nanmedian
import warnings
warnings.filterwarnings('ignore')

def subtract_background(image, degree=1, fit_regions=None):

    ''' Subtracts a background from an image. Can supply either pixel
    coordinates of the image or DS9 regions.
    '''

    image = image.squeeze()
    if len(image.shape) != 2:
        raise ValueError('Image must be a 2D numpy.ndarray')



    # Create grid for the fitting
    x_len, y_len = image.shape[1], image.shape[0]

    x_axis = np.linspace(0, x_len, x_len)
    y_axis = np.linspace(0, y_len, y_len)

    x, y = np.meshgrid(x_axis, y_axis)

    # Mask the image with the fit regions
    if fit_regions is None:
    	fit_regions = (0, 0, image.shape[0] - 1, image.shape[1] - 1)

    mask = np.ones(image.shape)

    for region in fit_regions:
    	mask[(x > region[1]) & \
    	     (x < region[3]) & \
    	     (y > region[0]) & \
    	     (y < region[2])] = 0

    image = np.ma.array(image, mask=mask)

    #plt.imshow(image)
    #plt.show()

    popt, pcov = mymath.polyfit2d(x, y, image, degree=degree)

    image_fit = mymath.poly2d((x, y), degree=degree, *popt)

    image -= image_fit

    return image

def read_ds9_region(filename):

    ''' Converts DS9 region file into format for plotting region.

    Need the following format:
        angle : degrees
        xy : pixels
        width : pixels
        height : pixels

    Region file provides following format:
        # Region file format: DS9 version 4.1
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        fk5
        box(4:17:04.740,+29:20:31.32,5854.33",11972.7",130) # text={test}

    pyregion module reads DS9 regions:
    http://leejjoon.github.io/pyregion/users/overview.html


    '''

    # Import external modules
    import pyregion as pyr

    # Read region file
    region = pyr.open(filename)

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    return region[0].coord_list

def calc_region_mask(filename, data, header, region_name='',):

    ''' Masks all pixels which are not within the region.

    '''

    import mygeometry as myg

    regions = load_ds9_region(filename, header, region_name=region_name)

    region_vertices = regions[region_name]['poly_verts']['pixel']

    # block off region
    region_mask = np.logical_not(myg.get_polygon_mask(data,
                                                      region_vertices))

    return region_mask

def load_ds9_region(filename, header=None, region_name=''):

    import pyregion as pyr

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    ds9_regions = {}
    region_indices = []
    for i, region in enumerate(regions):
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        name = tag[tag.find('text={')+6:tag.find('}')].lower()

        if name == region_name:
            region_indices.append(i)

    for i in region_indices:
        region = regions[i]

        # Format vertices to be 2 x N array
        poly_verts = []
        for i in xrange(0, len(region.coord_list)/2):
            poly_verts.append((region.coord_list[2*i],
                               region.coord_list[2*i+1]))

        ds9_regions[region_name] = {}
        ds9_regions[region_name]['poly_verts'] = {}
        ds9_regions[region_name]['poly_verts']['wcs'] = poly_verts

        if header is not None:
            poly_verts_pix = []
            for i in xrange(0, len(poly_verts)):
                poly_verts_pix.append(\
                        get_pix_coords(ra=poly_verts[i][0],
                                       dec=poly_verts[i][1],
                                       header=header)[:-1][::-1].tolist())

            ds9_regions[region_name]['poly_verts']['pixel'] = poly_verts_pix

    return ds9_regions

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec).
    '''

    import pywcsgrid2 as wcs
    import pywcs
    from astropy.wcs import WCS

    # convert to degrees
    if type(ra) is tuple or type(ra) is list:
        ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
    else:
    	ra_deg, dec_deg = ra, dec

    wcs_header = WCS(header)
    pix_coords = wcs_header.wcs_world2pix([[ra_deg, dec_deg],], 0)[0]

    return np.hstack((pix_coords, -1))

def hrs2degs(ra=None, dec=None):
    ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
    '''

    ra_deg = 15*(ra[0] + ra[1]/60. + ra[2]/3600.)
    dec_deg = dec[0] + dec[1]/60. + dec[2]/3600.

    return (ra_deg, dec_deg)

def bin_image(ndarray, binsize=(1,1), header=None, origin='lower left',
        statistic=np.nansum, return_weights=False, weights=None,
        quick_bin=False):

    ''' Bins an image.

    Parameters
    ----------
    image : array-like
        2D image or 3D to be binned. First two axes will be binned.
    binsize : tuple
        Integer number of pixels to bin in x,y directions
    header : astropy.io.Header, optional
        Header
    origin : string
        Location to begin binning.
    statistic : function
        Function with which to operate on the pixels being binned.
    weights : array-like, optional
        Weights to be used in the statistic. The statistic must have an
        argument 'weights' if this argument is used.
    return_weights : bool, optional
        Return weights for image corresponding to the number of valid pixels
        used to create a binned pixel?

    Returns
    -------
    image_binned : array-like
        Binned image.
    header_binned : astropy.io.Header, optional
        If header is provided, edited header is returned.

    http://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array

    '''

    import numpy as np
    from astropy.io import fits
    from scipy.stats import binned_statistic_2d
    from scipy.stats import binned_statistic_dd

    #if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
    #    raise ValueError("Operation not supported.")

    binsize = np.array(binsize, dtype=int)

    # write new shape given binsize
    new_shape = np.array(ndarray.shape) / binsize
    shape_excess = np.array(ndarray.shape) % binsize

    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))

    # Trim image to allow for binning
    for axis in xrange(ndarray.ndim):
        axis_length = ndarray.shape[axis]
        excess = shape_excess[axis]

        # trim elements
        trim_indices = np.arange(excess)

        if excess % 2 == 1:
            end_indices = range(axis_length - excess / 2 - 1, axis_length)
        else:
            end_indices = range(axis_length - excess / 2, axis_length)

        trim_indices[excess / 2: excess] = end_indices

        ndarray = np.delete(ndarray, trim_indices, axis)
        #statistic = lambda x: [print(element) for element in x]

    if quick_bin:
        if len(binsize) == 2:
            # reshape trimmed array
            compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                          ndarray.shape)]
            flattened = [l for p in compression_pairs for l in p]

            ndarray = ndarray.reshape(flattened)
            for i in range(len(new_shape)):
                ndarray = statistic(ndarray, (-1*(i+1)))
            ndarray_bin = ndarray
        elif len(binsize) == 3:
            ndarray_bin = np.empty((ndarray.shape[0],
                                    new_shape[1],
                                    new_shape[2]))

            for i in xrange(ndarray.shape[0]):
                # reshape trimmed array
                compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                              ndarray.shape)]
                flattened = [l for p in compression_pairs for l in p]

                ndarray = ndarray.reshape(flattened)
                for i in range(len(new_shape)):
                    ndarray = statistic(ndarray, (-1*(i+1)))

                ndarray_bin[i, :, :] = ndarray
    else:
        if len(binsize) == 2:
            xx, yy = np.meshgrid(np.arange(0, ndarray.shape[0], 1),
                                 np.arange(0, ndarray.shape[1], 1))

            ndarray_bin = binned_statistic_2d(xx.T.ravel(),
                                              yy.T.ravel(),
                                              values=ndarray.ravel(),
                                              bins=new_shape,
                                              statistic=statistic)[0]

            if return_weights:
                count = lambda x: np.size(x[~np.isnan(x)])

                weights = binned_statistic_2d(xx.T.ravel(),
                                              yy.T.ravel(),
                                              values=ndarray.ravel(),
                                              bins=new_shape,
                                              statistic=count)[0]
        if len(binsize) == 3:
            xx, yy = np.meshgrid(np.arange(0, ndarray.shape[1], 1),
                                 np.arange(0, ndarray.shape[2], 1))

            ndarray_bin = np.empty((ndarray.shape[0],
                                    new_shape[1],
                                    new_shape[2]))

            for i in xrange(ndarray.shape[0]):
                ndarray_bin[i, :, :] = binned_statistic_2d(xx.T.ravel(),
                                            yy.T.ravel(),
                                            values=ndarray[i, :, :].ravel(),
                                            bins=new_shape[1:],
                                            statistic=statistic)[0]

            if return_weights:
                weights = np.empty((ndarray.shape[0],
                                    new_shape[1],
                                    new_shape[2]))
                count = lambda x: np.size(x[~np.isnan(x)])
                for i in xrange(ndarray.shape[0]):

                    weights[i, :, :] = binned_statistic_2d(xx.T.ravel(),
                                                yy.T.ravel(),
                                                values=ndarray[i, :, :].ravel(),
                                                bins=new_shape[1:],
                                                statistic=count)[0]
    if 0:
        if weights is not None:
            def new_statistic(slice):
                return statistic(ndarray[slice], weights=weights[slice])
        else:
            def new_statistic(slice):
                return statistic(ndarray[slice])

        ndarray_bin = np.empty(new_shape)

        if len(binsize) <= 2:
            if binsize.size == 1:
                binsize = np.array(binsize, binsize)

            for j in xrange(new_shape[0]):
                for k in xrange(new_shape[1]):
                    indices = (range(j*binsize[0], (j+1)*binsize[0]),
                               range(k*binsize[1], (k+1)*binsize[1]))
                    print indices
                    ndarray_bin[j,k] = new_statistic(indices)
        if len(binsize) == 3:
            for i in xrange(new_shape[0]):
                for j in xrange(new_shape[1]):
                    for k in xrange(new_shape[2]):
                        indices = (
                                   range(i*binsize[0], (i+1)*binsize[0]),
                                   range(j*binsize[1], (j+1)*binsize[1]),
                                   range(k*binsize[2], (k+1)*binsize[2]))
                        ndarray_bin[i, j,k] = new_statistic(indices)

    # Edit header
    if header is not None:
        header_bin = header.copy()

        header_bin['NAXIS1'] = new_shape[0]
        header_bin['NAXIS2'] = new_shape[1]
        #header_bin['CRVAL1'] = header['CRVAL1'] / binsize[-2] + \
        header_bin['CRVAL1'] = header['CRVAL1'] + \
                               (excess / 2) * header['CDELT1']
        #header_bin['CRVAL2'] = header['CRVAL2'] / binsize[-1] + \
        header_bin['CRVAL2'] = header['CRVAL2'] + \
                               (excess / 2) * header['CDELT2']
        header_bin['CDELT1'] = header['CDELT1'] * binsize[-2]
        header_bin['CDELT2'] = header['CDELT2'] * binsize[-1]
        header_bin['CRPIX1'] = header['CRPIX1'] / binsize[-2]
        header_bin['CRPIX2'] = header['CRPIX2'] / binsize[-1]

        if return_weights:
            result = ndarray_bin, header_bin, weights
        else:
            result = ndarray_bin, header_bin
    else:
        if return_weights:
            result = ndarray_bin, weights
        else:
            result = ndarray_bin

    return result

def _bin_image(image, binsize=(1,1), header=None, origin='lower left',
        func=np.nansum, return_weights=False):

    ''' Bins an image.

    Parameters
    ----------
    image : array-like
        2D image or 3D to be binned. First two axes will be binned.
    binsize : tuple
        Integer number of pixels to bin in x,y directions
    header : astropy.io.Header, optional
        Header
    origin : string
        Location to begin binning.
    func : function
        Function with which to operate on the pixels being binned. Must accept
        image as first argument and an axis to bin as the second argument.
    return_weights : bool, optional
        Return weights for image corresponding to the number of valid pixels
        used to create a binned pixel?

    Returns
    -------
    image_binned : array-like
        Binned image.
    header_binned : astropy.io.Header, optional
        If header is provided, edited header is returned.

    '''

    import numpy as np
    from astropy.io import fits



    # Cycle through slices of image, if 2D, only one slice
    if image.ndim == 2:
        # Write axes grids
        x_grid_bin = np.arange(0, image.shape[0], binsize[0], dtype=int)
        y_grid_bin = np.arange(0, image.shape[1], binsize[1], dtype=int)

        # Bin image
        image_binned = np.zeros((len(x_grid_bin), len(y_grid_bin)))

        for i in xrange(0, len(x_grid_bin) - 1):
            for j in xrange(0, len(y_grid_bin) - 1):
                image_binned[i, j] = \
                    func(image[x_grid_bin[i]:x_grid_bin[i+1],
                               y_grid_bin[j]:y_grid_bin[j+1]])

        # Get edge pixels
        for i in xrange(0, len(x_grid_bin) - 1):
            image_binned[i, -1] = \
                func(image[x_grid_bin[i]:x_grid_bin[i+1],
                           y_grid_bin[-1]:-1])
        for j in xrange(0, len(y_grid_bin) - 1):
            image_binned[i, j] = \
                    func(image[x_grid_bin[-1]:-1,
                               y_grid_bin[j]:y_grid_bin[j+1]])
        image_binned[-1, -1] = \
                func(image[x_grid_bin[-1]:-1,
                           y_grid_bin[-1]:-1])

        if return_weights:
            # Bin image
            image_weights = np.zeros((len(x_grid_bin), len(y_grid_bin)))

            for i in xrange(0, len(x_grid_bin) - 1):
                for j in xrange(0, len(y_grid_bin) - 1):
                    image_weights[i, j] = \
                        np.sum(~np.isnan(image[x_grid_bin[i]:x_grid_bin[i+1],
                                               y_grid_bin[j]:y_grid_bin[j+1]]))

            # Get edge pixels
            for i in xrange(0, len(x_grid_bin) - 1):
                image_weights[i, -1] = \
                    np.sum(~np.isnan(image[x_grid_bin[i]:x_grid_bin[i+1],
                                           y_grid_bin[-1]:-1]))
            for j in xrange(0, len(y_grid_bin) - 1):
                image_weights[i, j] = \
                        np.sum(~np.isnan(image[x_grid_bin[-1]:-1,
                                               y_grid_bin[j]:y_grid_bin[j+1]]))
            image_weights[-1, -1] = \
                    np.sum(~np.isnan(image[x_grid_bin[-1]:-1,
                                           y_grid_bin[-1]:-1]))

            image_weights /= np.nansum(image_weights)

    elif image.ndim == 3:
        # Write axes grids
        x_grid_bin = np.arange(0, image.shape[1], binsize[0], dtype=int)
        y_grid_bin = np.arange(0, image.shape[2], binsize[1], dtype=int)

        # Bin image
        image_binned = np.zeros((image.shape[0],
                                 len(x_grid_bin),
                                 len(y_grid_bin),)
                                )

        for k in xrange(0, image.shape[0]):
            for i in xrange(0, len(x_grid_bin) - 1):
                for j in xrange(0, len(y_grid_bin) - 1):
                    image_binned[k, i, j] = \
                            func(image[k, x_grid_bin[i]:x_grid_bin[i+1],
                                          y_grid_bin[j]:y_grid_bin[j+1]])

            # Get edge pixels
            for i in xrange(0, len(x_grid_bin) - 1):
                image_binned[k, i, -1] = \
                    func(image[k, x_grid_bin[i]:x_grid_bin[i+1],
                                  y_grid_bin[-1]:-1])
            for j in xrange(0, len(y_grid_bin) - 1):
                image_binned[k, -1, j] = \
                        func(image[k, x_grid_bin[-1]:-1,
                                      y_grid_bin[j]:y_grid_bin[j+1]])
                image_binned[k, -1, -1] = \
                    func(image[k, x_grid_bin[-1]:-1,
                                  y_grid_bin[-1]:-1])


        if return_weights:
            # Bin image
            image_weights = np.zeros((image.shape[0],
                                     len(x_grid_bin),
                                     len(y_grid_bin),)
                                    )

            for k in xrange(0, image.shape[0]):
                for i in xrange(0, len(x_grid_bin) - 1):
                    for j in xrange(0, len(y_grid_bin) - 1):
                        image_weights[k, i, j] = \
                            np.sum(~np.isnan(image[k, x_grid_bin[i]:x_grid_bin[i+1],
                                              y_grid_bin[j]:y_grid_bin[j+1]]))

                # Get edge pixels
                for i in xrange(0, len(x_grid_bin) - 1):
                    image_weights[k, i, -1] = \
                        np.sum(~np.isnan(image[k, x_grid_bin[i]:x_grid_bin[i+1],
                                      y_grid_bin[-1]:-1]))
                for j in xrange(0, len(y_grid_bin) - 1):
                    image_weights[k, -1, j] = \
                        np.sum(~np.isnan(image[k, x_grid_bin[-1]:-1,
                                          y_grid_bin[j]:y_grid_bin[j+1]]))
                    image_weights[k, -1, -1] = \
                        np.sum(~np.isnan(image[k, x_grid_bin[-1]:-1,
                                      y_grid_bin[-1]:-1]))

            image_weights /= np.nansum(image_weights)

    # Edit header
    if header is not None:
        header_bin = header.copy()

        header_bin['NAXIS1'] = len(y_grid_bin)
        header_bin['NAXIS2'] = len(x_grid_bin)
        header_bin['CDELT1'] = header['CDELT1'] * binsize[0]
        header_bin['CDELT2'] = header['CDELT2'] * binsize[1]
        header_bin['CRPIX1'] = header['CRPIX1'] / binsize[0]
        header_bin['CRPIX2'] = header['CRPIX2'] / binsize[1]

        if return_weights:
            result = image_binned, header_bin, image_weights
        else:
            result = image_binned, header_bin
    else:
        if return_weights:
            result = image_binned, image_weights
        else:
            result = image_binned

    return result

def calculate_nhi(cube=None, velocity_axis=None, velocity_range=None,
        return_nhi_error=False, noise_cube=None,
        velocity_noise_range=[90,100], Tsys=30., header=None,
        fits_filename=None, fits_error_filename=None, verbose=False):

    ''' Calculates an N(HI) image given a velocity range within which to
    include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used. If the cube is two dimensional,
        assumes the spatial dimension has been raveled into the 1st axis.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to
    fits_filename : str
        If specified, and a header is provided, the N(HI) image will be written
        as a fits file.
    header : pyfits.Header
        Header from cube.
    velocity_range : tuple, ndarray
        If velocity range is a tuple including images of the same
        shape as the spatial dimensions of the cube, then each velocity range
        will be applied to individual pixels
    velocity_noise_range : tuple
        Lower and upper velocities which do not have emission. Used to
        calculate the N(HI) noise map.

    Returns
    -------
    nhi_image : array-like
        N(HI) in units of 10^20 cm^-2
    nhi_image_error : array-like, optional
        N(HI) error in units of 10^20 cm^-2
    '''

    import numpy as np

    if velocity_range is None:
        velocity_range = np.min(velocity_axis), np.max(velocity_axis)

    if type(velocity_range) is not np.ndarray:
        velocity_range = np.asarray(velocity_range)

    #if velocity_range.ndim != 1 | velocity_range.ndim != 2:
    #    raise ValueError('velocity_range must be pair or a list of pairs' + \
    #            ' of HI velocities.')

    if cube is None:
    	raise TypeError('cube is required')
    if velocity_axis is None:
    	raise TypeError('velocity_axis is required')

    # Calculate NHI from cube if set
    # If the cube is 3 dimensional, one velocity, and two spatial axes
    if cube.ndim == 3:
        image = np.empty((cube.shape[1],
                          cube.shape[2]))
        image[:,:] = np.NaN

        # If multiple velocity ranges, then include velocities between each
        # one
        if velocity_range.ndim == 3:
            if velocity_range.shape[1] == cube.shape[1]:
                #raise ValueError('3D multi vel range not yet implemented')
                for i in xrange(0, cube.shape[1]):
                    for j in xrange(0, cube.shape[2]):
                        indices = (velocity_range[0, i, j] <= velocity_axis) & \
                                  (velocity_range[1, i, j] >= velocity_axis)
                        image[i, j] = np.nansum(cube[indices, i, j])
                # Calculate image error
                if return_nhi_error:
                    image_error = np.empty((cube.shape[1],
                                            cube.shape[2]))
                    for i in xrange(0, cube.shape[1]):
                        for j in xrange(0, cube.shape[2]):
                            indices = \
                                (velocity_range[0, i, j] <= velocity_axis) & \
                                (velocity_range[1, i, j] >= velocity_axis)
                        image_error[i, j] = \
                                (noise_cube[indices, i, j]**2).sum(axis=0)**0.5
            else:
                raise ValueError('Incorrect velocity range shape')

        # If only one velocity range provided, integrate from lower to upper
        elif velocity_range.ndim == 1:
            indices = np.where((velocity_axis >= velocity_range[0]) & \
                               (velocity_axis <= velocity_range[1]))[0]

            image[:,:] = np.nansum(cube[indices,:,:], axis=0)

            # Calculate image error
            if return_nhi_error:
                image_error = np.empty((cube.shape[1],
                                  cube.shape[2]))
                image_error[:,:] = np.NaN
                image_error[:,:] = (noise_cube[indices,:,:]**2).sum(axis=0)**0.5
    # If the cube is unraveled, one dimension velocity, and one dimension
    # spatial
    if cube.ndim == 2:
        image = np.empty((cube.shape[1]))
        # If only one velocity range provided, integrate whole cube with one
        # range
        if velocity_range.ndim == 1:
            image[:] = np.NaN
            indices = np.where((velocity_axis >= velocity_range[0]) & \
                               (velocity_axis <= velocity_range[1]))[0]
            image[:] = np.nansum(cube[indices,:], axis=0)
        # If an image of velocity ranges provided, integrate each pixel with a
        # different range
        elif velocity_range.shape[1] == cube.shape[1]:
            for i in xrange(0, cube.shape[1]):
                indices = (velocity_range[0, i] <= velocity_axis) & \
                          (velocity_range[1, i] >= velocity_axis)
                image[i] = np.nansum(cube[indices, i])

        # Calculate image error
        if return_nhi_error:
            image_error = np.empty((cube.shape[1]))
            image_error[:] = np.NaN
            image_error[:] = (noise_cube[indices,:]**2).sum(axis=0)**0.5

    # NHI in units of 1e20 cm^-2
    delta_v = velocity_axis[-1] - velocity_axis[-2]
    nhi_image = image * 1.823e-2 * delta_v

    if fits_filename is not None and header is not None:
        if verbose:
            print('Writing N(HI) image to FITS file %s' % fits_filename)
        header['BUNIT'] = '1e20 cm^-2'
        header.remove('CDELT3')
        header.remove('CRVAL3')
        header.remove('CRPIX3')
        header.remove('CTYPE3')
        header.remove('NAXIS3')
        header['NAXIS'] = 2

        fits.writeto(fits_filename, image*1.823e-2, header = header, clobber =
                True, output_verify = 'fix')

    if fits_error_filename is not None and header is not None:
        if verbose:
            print('Writing N(HI) error image to FITS file %s' % fits_filename)

        fits.writeto(fits_error_filename, image_error * 1.823e-2, header =
                header, clobber = True, output_verify = 'fix')

    if return_nhi_error:
        nhi_image_error = image_error * 1.823e-2 * delta_v
        return nhi_image, nhi_image_error
    else:
        return nhi_image

def calculate_noise_cube(cube=None, velocity_axis=None,
            velocity_noise_range=[-110,-90,90,110], header=None, Tsys=30.,
            filename=None):

    """ Calcuates noise envelopes for each pixel in a cube
    """

    import numpy as np
    from astropy.io import fits

    noise_cube = np.zeros(cube.shape)
    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            profile = cube[:,i,j]
            noise = calculate_noise(profile, velocity_axis,
                    velocity_noise_range)
            #noise = 0.1 # measured in line free region
            noise_cube[:,i,j] = calculate_noise_scale(Tsys,
                    profile, noise=noise)

    if filename is not None:
        fits.writeto(filename, noise_cube, header=header)

    return noise_cube

def calculate_noise(profile, velocity_axis, velocity_range):
    """ Calculates rms noise of Tile profile given velocity ranges.
    """
    import numpy as np

    std = 0

    # calculate noises for each individual region
    for i in range(len(velocity_range) / 2):
        velMin = velocity_range[2*i + 0]
        velMax = velocity_range[2*i + 1]

        noise_region = np.where((velocity_axis >= velMin) & \
                        (velocity_axis <= velMax))

        std += np.std(profile[noise_region])

    std /= len(velocity_range) / 2
    return std

def calculate_noise_scale(Tsys, profile, noise=None):
    """ Creates an array for scaling the noise by (Tsys + Tb) / Tsys
    """
    import numpy as np
    n = np.zeros(len(profile))
    n = (Tsys + profile) / Tsys * noise

    return n

def calculate_sd(image, sd_factor=1/1.25):

    ''' Calculates a surface density image given a velocity range within which
    to include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to '''

    import numpy as np

    # NHI in units of 1e20 cm^-2
    sd_image = image * sd_factor

    return sd_image

def calculate_nh2(nhi_image = None, av_image = None, dgr = 1.1e-1):

    ''' Calculates the total gas column density given N(HI), A_v and a
    dust-to-gas ratio.

    Parameters
    ----------
    '''

    import numpy as np

    nh2_image = 0.5 * (av_image / dgr - nhi_image)

    return nh2_image

def calculate_nh2_error(nhi_image_error=None, av_image_error=None, dgr=1.1e-1):

    ''' Calculates the total gas column density given N(HI), A_v and a
    dust-to-gas ratio.

    Parameters
    ----------
    '''

    import numpy as np

    nh2_image_error = 0.5 * ((av_image_error / dgr)**2 - \
            nhi_image_error**2)**0.5

    return nh2_image_error

def fit_background(av_data, background_mask=None, background_dim=1):

    from scipy.interpolate import interp2d
    from scipy.interpolate import SmoothBivariateSpline as spline

    if background_mask is None:
        background_mask = np.zeros(av_data.shape)

    if background_dim == 1:
        background = np.nanmean(av_data[~background_mask])

    if background_dim == 2:
        #av_data = np.ma.array(av_data, mask=background_mask)

        loc = np.where(~background_mask)

        x = loc[0]
        y = loc[1]

        z = av_data[~background_mask]

        #print av_data.size, z.shape
        assert z.shape == x.shape

        bbox = [0, av_data.shape[0], 0, av_data.shape[1]]

        result_interp = spline(x, y, z, bbox=bbox, kx=1, ky=1)

        x_grid, y_grid = np.where(av_data)

        background_flat = result_interp.ev(x_grid, y_grid)

        background = np.reshape(background_flat, av_data.shape)

    return background

def calc_spectrum(cube, statistic=nanmedian, mask=None):

    ''' Calculates the spectrum of a cube. Assumes velocity axis is 0th
    dimension of cube.

    '''

    if mask is None:
        mask = np.zeros(cube.shape[1:], dtype=bool)

    spectrum = statistic(cube[:, ~mask], axis=1)

    return spectrum

def get_radial_profile(image, center=None, stddev=False, binsize=1,
        mask=None, weights=None):

    ''' Calculates radial profiles of an image at the center.

    '''

    # import external modules
    import numpy as np
    from scipy.optimize import curve_fit
    from agpy import azimuthalAverage as radial_average

    if stddev and weights is not None:
        weights=None

    result = radial_average(image, binsize=binsize, center=center,
            stddev=stddev, mask=mask, interpnan=False, returnradii=True,
            weights=weights)

    return result

def crop_cube(cube, vel_axis, vel_range):

    ''' Crops a cube to include only velocities within vel_range.

    '''

    crop_indices = (vel_axis > vel_range[0]) & (vel_axis < vel_range[1])

    cube_crop = cube[crop_indices, :, :]
    vel_axis_crop = vel_axis[crop_indices]

    return cube_crop, vel_axis_crop

def calculate_moment(data, moment=0, weights=1):

    ''' Calculates moment of data. Weights needed for moment of 1 or 2.
    '''

    data, weights = np.squeeze(data), np.squeeze(weights)

    if moment == 0:
        result = np.nanmean(data)
    elif moment == 1:
        result = np.zeros(data.shape[1:])
        for i in xrange(result.shape[0]):
            for j in xrange(result.shape[1]):
                result[i,j] = \
                    np.nansum(data[:,i,j] * weights)/np.nansum(data[:,i,j])
    elif moment == 2:
        M1 = np.nansum(data*weights)/np.nansum(data),
        result = np.sqrt(np.nansum(data*(weights-M1)**2)/np.nansum(data))

    return result


