#!/usr/bin/python

''' Module for image analysis in astro applications.

'''

import numpy as np
import mymath
import mycoords
import matplotlib.pyplot as plt
import pyfits as pf

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

def load_ds9_region(filename, header=None):

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    for core in cores:
    	region = read_ds9_region(filename_base + core + '.reg')
        box_center_pixel = get_pix_coords(ra = region[0],
                                          dec = region[1],
                                          header = header)
        box_center_pixel = (int(box_center_pixel[1]), int(box_center_pixel[0]))
        box_height = region[2] / header['CDELT1']
        box_width = region[3] / header['CDELT2']
        cores[core].update({'box_center_pix': box_center_pixel})
        cores[core].update({'box_width': box_width})
        cores[core].update({'box_height': box_height})
        cores[core].update({'box_angle': region[4]})

    return cores

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec).
    '''

    import pywcsgrid2 as wcs
    import pywcs

    # convert to degrees
    if type(ra) is tuple and type(dec) is tuple:
        ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
    else:
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

def bin_image(image, binsize=(1,1), header=None, origin='lower left',
        func=np.nansum):

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
        Function with which to operate on the pixels being binned.

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
        x_grid_bin = np.arange(0, image.shape[0], binsize[0])
        y_grid_bin = np.arange(0, image.shape[1], binsize[1])

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
            image_binned[-1, j] = \
                    func(image[x_grid_bin[-1]:-1,
                               y_grid_bin[j]:y_grid_bin[j+1]])
        image_binned[-1, -1] = \
                func(image[x_grid_bin[-1]:-1,
                           y_grid_bin[-1]:-1])

    elif image.ndim == 3:
        # Write axes grids
        x_grid_bin = np.arange(0, image.shape[1], binsize[0])
        y_grid_bin = np.arange(0, image.shape[2], binsize[1])

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
    # Edit header
    if header is not None:
        header_bin = header.copy()

        header_bin['NAXIS1'] = len(y_grid_bin)
        header_bin['NAXIS2'] = len(x_grid_bin)
        header_bin['CDELT1'] = header['CDELT1'] * binsize[0]
        header_bin['CDELT2'] = header['CDELT2'] * binsize[1]
        header_bin['CRPIX1'] = header['CRPIX1'] / binsize[0]
        header_bin['CRPIX2'] = header['CRPIX2'] / binsize[1]

        result = image_binned, header_bin
    else:
        result = image_binned

    return result

def calculate_nhi(cube=None, velocity_axis=None, velocity_range=[],
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
        if velocity_range.ndim == 2:
            if velocity_range[0].ndim == cube.shape[1:]:
                raise ValueError('3D multi vel range not yet implemented')
                '''
                indices = np.zeros(cube.shape)
                for i in xrange(0, indices.shape[0]):
                    for j in xrange(0, indices.shape[1]):
                        indices[:, i,j] = (velocity_range[0, i, j] < \
                                        velocity_axis) & \
                                       (velocity_range[1, i, j] > \
                                        velocity_axis)
                indices = (velocity_range[0] < velocity_axis) & \
                          (velocity_range[1] > velocity_axis)
                indices = np.where(indices == 1)[0]
                '''
                image[:,:] = cube[indices].sum(axis=0)
            else:
                indices = np.zeros(cube.shape[0])
                for i in xrange(0, velocity_range.shape[0]):
                    indices += (velocity_axis > velocity_range[i, 0]) & \
                                (velocity_axis < velocity_range[i, 1])
                indices = np.where(indices == 1)[0]
                image[:,:] = cube[indices,:,:].sum(axis=0)
        # If only one velocity range provided, integrate from lower to upper
        elif velocity_range.ndim == 1:
            indices = np.where((velocity_axis > velocity_range[0]) & \
                               (velocity_axis < velocity_range[1]))[0]

            image[:,:] = cube[indices,:,:].sum(axis=0)

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
            indices = np.where((velocity_axis > velocity_range[0]) & \
                    (velocity_axis < velocity_range[1]))[0]
            image[:] = cube[indices,:].sum(axis=0)
        # If an image of velocity ranges provided, integrate each pixel with a
        # different range
        elif velocity_range[0].shape[0] == cube.shape[1]:
            image = np.zeros(cube.shape[1])
            for i in xrange(0, cube.shape[1]):
                indices = (velocity_range[0, i] < velocity_axis) & \
                          (velocity_range[1, i] > velocity_axis)
                image[i] = np.sum(cube[indices, i])

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

        pf.writeto(fits_filename, image*1.823e-2, header = header, clobber =
                True, output_verify = 'fix')

    if fits_error_filename is not None and header is not None:
        if verbose:
            print('Writing N(HI) error image to FITS file %s' % fits_filename)

        pf.writeto(fits_error_filename, image_error * 1.823e-2, header =
                header, clobber = True, output_verify = 'fix')

    if return_nhi_error:
        nhi_image_error = image_error * 1.823e-2
        return nhi_image, nhi_image_error
    else:
        return nhi_image

def calculate_noise_cube(cube=None, velocity_axis=None,
            velocity_noise_range=[-110,-90,90,110], header=None, Tsys=30.,
            filename=None):

    """ Calcuates noise envelopes for each pixel in a cube
    """

    import numpy as np
    import pyfits as pf

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
        pf.writeto(filename, noise_cube, header=header)

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

    nh_image = 0.5 * (av_image / dgr - nhi_image)

    return nh_image

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
