#!/usr/bin/python -W ignore::DeprecationWarning

# Import external modules
import numpy as np
import warnings
warnings.filterwarnings('ignore')

'''
Module Classes

'''

class Cloud():

    def __init__(self, av_filename,
            hi_filename,
            av_error_filename=None,
            hi_error_filename=None,
            cloud_prop_filename=None,
            dgr_grid=None,
            intercept_grid=None,
            width_grid=None,
            threshold_delta_dgr=0.0001,
            residual_width_scale=3.0,
            av_error=None,
            av_background=None,
            verbose=False,
            clobber_likelihoods=True,
            likelihood_filename='',
            hi_noise_vel_range=[0,1],
            binsize=1,
            vel_range_diff_thres=2.0,
            init_vel_range=[0,1],
            pixel_mask_increase_fraction=0.01,
            diagnostic_filename=None,
            use_bin_weights=False,
            use_only_binned_data=False,
            binned_data_filename_ext=None,
            weights_filename=None,
            plot_args={}):

        # Import external modules
        import pyfits as fits
        from mycoords import make_velocity_axis
        from os import system,path

        # Define local variables
        self.av_filename = av_filename
        self.hi_filename = hi_filename
        self.av_error_filename = av_error_filename
        self.hi_error_filename = hi_error_filename
        self.av_error = av_error
        self.dgr_grid = dgr_grid
        self.width_grid = width_grid
        self.intercept_grid = intercept_grid
        self.THRESHOLD_DELTA_DGR = threshold_delta_dgr
        self.RESIDUAL_WIDTH_SCALE = residual_width_scale
        self.verbose = verbose
        self.clobber_likelihoods = clobber_likelihoods
        self.likelihood_filename = likelihood_filename
        self.hi_noise_vel_range = hi_noise_vel_range
        self.binsize = binsize
        self.use_bin_weights=use_bin_weights
        self.use_only_binned_data=use_only_binned_data
        self.Tsys = 30.0 # K
        self.VEL_RANGE_DIFF_THRES = vel_range_diff_thres
        self.init_vel_range = init_vel_range
        self.plot_args = plot_args
        self.diagnostic_filename = diagnostic_filename
        self.PIXEL_MASK_INCREASE_FRACTION = pixel_mask_increase_fraction
        self.iter_vars = {}
        self.binned_data_filename_ext = binned_data_filename_ext
        self.weights_filename = weights_filename

        # Initialize empty variables
        self.av_data_bin = None
        self.av_error_data_bin = None
        self.hi_data_bin = None
        self.hi_error_data_bin = None

        # Load data
        self.av_data, self.av_header = fits.getdata(av_filename, header=True)
        self.hi_data, self.hi_header = fits.getdata(hi_filename, header=True)
        if av_error_filename is not None:
            self.av_error_data, self.av_error_header = \
                    fits.getdata(av_error_filename, header=True)
        elif av_error is not None:
            self.av_error_data, self.av_error_header = \
                    np.ones(self.av_data.shape) * av_error, self.av_header
        else:
            self.av_error_data = None
        if hi_error_filename is not None and path.isfile(hi_error_filename):
            self.hi_error_data, self.hi_error_header = \
                fits.getdata(hi_error_filename, header=True)
        else:
            self.hi_error_data = None
        if cloud_prop_filename is not None:
            self.load_cloud_properties(cloud_prop_filename)

        # Use only binned data throughout the analysis?
        # Either load the binned images or create new images
        if self.use_only_binned_data:
            filenames = self._prepare_bin_filenames()
            if binned_data_filename_ext is not None:
                load_data = self._check_files(filenames)
                write_data = True
            else:
                load_data = False
                write_data = False

            if not load_data:
                self._bin_data(write_data=write_data)
            else:
                self._load_bin_data()

            self.av_data = self.av_data_bin
            self.av_error_data = self.av_error_data_bin
            self.hi_data = self.hi_data_bin
            self.av_header = self.av_header_bin.copy()
            self.av_error_header = self.av_error_header_bin.copy()
            self.hi_header = self.hi_header_bin.copy()

        if av_background is not None:
            self.av_background = av_background
            self.av_data = self.av_data - av_background

        # Make velocity axis for HI cube
        self.hi_vel_axis = make_velocity_axis(self.hi_header)

        # Derive velocity center
        self._calc_vel_center(self.hi_data)

        if self.diagnostic_filename is not None:
            import sys
            import os

            print('\nWriting output to \n' + self.diagnostic_filename)
            os.system('rm -rf ' + self.diagnostic_filename)

            self._orig_stdout = sys.stdout
            self._diagnostics = file(self.diagnostic_filename, 'a')
            sys.stdout = self._diagnostics

        else:
            self._diagnostics = None

    def _load_bin_data(self,):

        import pyfits as fits

        print('\n\tLoading binned data...')

        self.av_data_bin, self.av_header_bin = \
                fits.getdata(self.av_filename_bin, header=True)
        self.av_error_data_bin, self.av_error_header_bin = \
                fits.getdata(self.av_error_filename_bin, header=True)
        self.hi_data_bin, self.hi_header_bin = \
                fits.getdata(self.hi_filename_bin, header=True)
        if self.use_bin_weights:
            self.bin_weights = \
                    fits.getdata(self.hi_filename_bin, header=False)
        else:
            self.bin_weights = None

    def _check_files(self, filenames):

        from os import path

        if self.use_bin_weights:
            filenames.append(self.weights_filename)

        existing = np.zeros(len(filenames))

        for i, filename in enumerate(filenames):
            existing[i] = path.isfile(filename)

        return np.alltrue(existing)

    def _prepare_bin_filenames(self,):

        ext = self.binned_data_filename_ext

        self.av_filename_bin = self.av_filename.replace('.fits', ext + '.fits')
        self.hi_filename_bin = self.hi_filename.replace('.fits', ext + '.fits')
        if self.av_error_filename is None:
            self.av_error_filename_bin = \
                    self.av_filename.replace('.fits', '_error' + ext + '.fits')
        else:
            self.av_error_filename_bin = \
                    self.av_filename.replace('.fits', ext + '.fits')

        return [self.av_filename_bin, self.av_error_filename_bin, \
                self.hi_filename_bin,]

    def _calc_vel_center(self, hi_data, single_vel_center=True):

        import numpy as np

        if single_vel_center:
            self.hi_spectrum = np.nansum(hi_data, axis=(1,2))
            vel_center = np.array((np.average(self.hi_vel_axis,
                                   weights=self.hi_spectrum**2),))[0]
        else:
            vel_center = np.zeros(hi_data.shape[1:])
            for i in xrange(0, hi_data.shape[1]):
                for j in xrange(0, hi_data.shape[2]):
                    self.hi_spectrum = hi_data[:, i, j]
                    self.hi_spectrum[np.isnan(self.hi_spectrum)] = 0.0
                    if np.nansum(self.hi_spectrum) > 0:
                        vel_center[i,j] = \
                                np.array((np.average(self.hi_vel_axis,
                                            weights=self.hi_spectrum**2),))[0]
                    else:
                        vel_center[i,j] = np.nan

            if mask is not None:
                vel_center = vel_center[~mask]

        self.vel_center = vel_center

    def _derive_region_mask(self, binned=False):

        import mygeometry as myg

        if not binned:
            # Derive relevant region
            region_vertices = \
                self.props['regions'][self.region]['poly_verts']['pixel']

            # block off region
            region_mask = np.logical_not(myg.get_polygon_mask(self.av_data,
                                                              region_vertices))

            self.region_mask = region_mask
        else:
            # Derive relevant region
            region_vertices = \
                self.props['regions'][self.region]['poly_verts_bin']['pixel']

            # block off region
            region_mask = np.logical_not(myg.get_polygon_mask(self.av_data_bin,
                                                              region_vertices))

            self.region_mask = region_mask

    def _prep_mask(self, binned=False):

        ''' Derives mask for NaNs and 0 errors.
        '''

        mask = (np.isnan(self.av_data) | \
                np.isnan(self.av_error_data) | \
               (self.av_error_data == 0) | \
                np.isnan(self.nhi_image) | \
                np.isnan(self.nhi_error_image) | \
               (self.nhi_error_image == 0))

        if 0:
            if self.av_data_bin is not None:
                mask_bin = (np.isnan(self.av_data_bin) | \
                            np.isnan(self.av_error_data_bin) | \
                           (self.av_error_data_bin == 0) | \
                            np.isnan(self.nhi_image_bin) | \
                            np.isnan(self.nhi_error_image_bin) | \
                           (self.nhi_error_image_bin == 0))

        if 0:
            if not binned:
                return mask
            else:
                return mask_bin

        return mask

    def _derive_bin_mask(self,):

        self.props['region_limit_bin'] = self.props['region_limit'].copy()
        self.props['plot_limit_bin'] = self.props['plot_limit'].copy()
        self.props = convert_limit_coordinates(self.props,
                                               header=self.av_header_bin,
                                               coords=('region_limit_bin',
                                                       'plot_limit_bin'))

        region_vertices = \
                self.props['regions'][self.region]['poly_verts']['pixel']

        # block off region
        region_mask = np.logical_not(myg.get_polygon_mask(av_data,
                                                          region_vertices))

        mask_bin = self._prep_mask(binned=True)

    def _get_pix_coords(self, ra=None, dec=None, header=None):

        ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec), or Ra in degrees
        and dec in degrees.
        '''

        import pywcsgrid2 as wcs
        import pywcs

        # convert to degrees if ra and dec are array-like
        try:
            if len(ra) == 3 and len(dec) == 3:
                ra_deg, dec_deg = self._hrs2degs(ra=ra, dec=dec)
            else:
                raise ValueError('RA and Dec must be in (hrs,min,sec) and' + \
                        ' (deg,arcmin,arcsec) or in degrees.')
        except TypeError:
            ra_deg, dec_deg = ra, dec

        wcs_header = pywcs.WCS(header)
        pix_coords = wcs_header.wcs_sky2pix([[ra_deg, dec_deg, 0]], 0)[0]

        return pix_coords

    def _hrs2degs(self, ra=None, dec=None):
        ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
        '''

        ra_deg = 15*(ra[0] + ra[1]/60. + ra[2]/3600.)
        dec_deg = dec[0] + dec[1]/60. + dec[2]/3600.

        return (ra_deg, dec_deg)

    def _convert_coordinates(self,
            coords=('region_limit','co_noise_limits','plot_limit'),
            header=None):

        ''' Converts WCS coordinates to pixel coordinates for a few sky
        positions.

        Parameters
        ----------
        header : fits.header
            If None, then uses self.av_header'''

        if header is None:
            header = self.av_header

        # Initialize pixel keys
        for coord in coords:
            self.props[coord].update({'pixel': []})

            if coord in ('region_limit',
                         'plot_limit',
                         'region_limit_bin',
                         'plot_limit_bin'):
                limit_wcs = self.props[coord]['wcs']

                for limits in limit_wcs:
                    # convert centers to pixel coords
                    limit_pixels = self._get_pix_coords(ra=limits[0],
                                                 dec=limits[1],
                                                 header=header)[:2].tolist()

                    self.props[coord]['pixel'].append(limit_pixels[0])
                    self.props[coord]['pixel'].append(limit_pixels[1])
            elif coord == 'co_noise_limits':
                region_limits = self.props[coord]['wcs']

                # Cycle through each region, convert WCS limits to pixels
                for region in region_limits:
                    region_pixels = []
                    for limits in region:
                        # convert centers to pixel coords
                        limit_pixels = self._get_pix_coords(ra=limits[0],
                                                      dec=limits[1],
                                                      header=header)[:2].tolist()
                        region_pixels.append(limit_pixels)

                    # Append individual regions back to CO noise
                    self.props[coord]['pixel'].append(region_pixels)

    def _calc_model_error(self, vel_range_max, dgr_max, intercept_max):

        from myimage_analysis import calculate_nhi

        # Calulate chi^2 for best fit models
        # ----------------------------------
        nhi_image_temp = calculate_nhi(cube=self.hi_data_bin,
                                       velocity_axis=self.hi_vel_axis,
                                       velocity_range=vel_range_max,
                                       return_nhi_error=False)

        self.av_image_model_bin = nhi_image_temp * dgr_max + intercept_max
        numerator = (self.av_data_bin - self.av_image_model_bin)**2
        denominator = np.sum(~np.isnan(self.av_data_bin)) - 3
        std = np.sqrt(np.nansum(numerator) / denominator)

        #std = np.nanstd(self.av_data_bin - self.av_image_model_bin)

        # Derive new av data error
        self.av_error_data_bin = std * np.ones(self.av_error_data_bin.shape)

        self.iter_vars[self.iter_step]['init_likelihood_results']['std'] = \
                std

    def _bin_data(self, write_data=False):

        from myimage_analysis import calculate_nhi, calculate_noise_cube, \
                bin_image
        import pyfits as fits

        if self.verbose:
            print('\n\tBinning data...')

        # Bin the images, retain only one bin_weight image since they are all
        # the same
        # -------------------------------------------------------------------
        binsize = self.binsize

        # Av image
        self.av_data_bin, self.av_header_bin, self.bin_weights = \
                bin_image(self.av_data,
                          binsize=(binsize, binsize),
                          header=self.av_header,
                          statistic=np.nanmean,
                          return_weights=True)

        # Av image error
        # Errors add in square
        # mean = sum(a_i) / n
        # error on mean = sqrt(sum(a_i**2 / n**2))
        noise_func = lambda x: np.nansum(x**2)**0.5 / \
                                     np.sum(~np.isnan(x))

        self.av_error_data_bin, self.av_error_header_bin = \
                bin_image(self.av_error_data,
                          binsize=(binsize, binsize),
                          header=self.av_error_header,
                          statistic=noise_func,)

        # Hi image
        self.hi_data_bin, self.hi_header_bin = \
                bin_image(self.hi_data,
                          binsize=(1, binsize, binsize),
                          header=self.hi_header,
                          statistic=np.nanmean)


        if write_data:
            fits.writeto(self.hi_filename_bin, self.hi_data_bin,
                         self.hi_header_bin, clobber=True)
            fits.writeto(self.av_filename_bin, self.av_data_bin,
                         self.av_header_bin, clobber=True)
            fits.writeto(self.av_error_filename_bin, self.av_error_data_bin,
                         self.av_error_header_bin, clobber=True)
            if self.use_bin_weights:
                fits.writeto(self.weights_filename, self.bin_weights,
                             self.av_header_bin, clobber=True)

    def _iterate_mle_calc(self, vel_range=None,):

        # Import module
        # --------------
        from myimage_analysis import calculate_nhi, calculate_noise_cube,\
                                     bin_image
        from os import system,path

        # Prep data
        # ---------
        self.iter_vars[self.iter_step] = {}
        self.iter_vars[self.iter_step]['vel_range'] = vel_range

        # derive the region mask
        self._derive_region_mask()

        # check if HI error cube present
        if self.hi_error_data is None:
            self.hi_error_data = \
                    calculate_noise_cube(cube=self.hi_data,
                                    velocity_axis=self.hi_vel_axis,
                                velocity_noise_range=self.hi_noise_vel_range,
                                    header=self.hi_header,
                                    Tsys=self.Tsys,
                                    filename=self.hi_error_filename)

        # Derive mask by excluding correlated residuals
        # ---------------------------------------------
        # Derive initial N(HI) image
        self.nhi_image, self.nhi_error_image = \
                calculate_nhi(cube=self.hi_data,
                              velocity_axis=self.hi_vel_axis,
                              velocity_range=vel_range,
                              noise_cube=self.hi_error_data,
                              velocity_noise_range=self.hi_noise_vel_range,
                              Tsys=self.Tsys,
                              return_nhi_error=True,
                              )

        # Write iteration step variables
        self.iter_vars[self.iter_step]['nhi_image'] = self.nhi_image
        self.iter_vars[self.iter_step]['nhi_image_error'] = \
                self.nhi_error_image

        # Finally, derive the mask to be used in calculation
        self._iterate_residual_masking()

        # delete next two lines
        if 0:
            mask = (self.av_data < 1.0)
            self.iter_vars[self.iter_step]['mask'] = mask
            self.plot_args['iter_ext'] = '0_0'

        # Apply mask to data, then bin data to avoid correlated pixels
        # ------------------------------------------------------------
        mask = self.iter_vars[self.iter_step]['mask']
        self.av_data[mask] = np.nan
        self.av_error_data[mask] = np.nan
        self.hi_data[:, mask] = np.nan

        # Bin the data
        if not self.use_only_binned_data:
            self._bin_data()

            # Write binned filenames
            self.av_bin_filename = \
                    self.av_filename.replace('.fits', '_bin.fits')
            self.av_error_bin_filename = \
                    self.hi_error_filename.replace('.fits', '_bin.fits')
            self.hi_bin_filename = \
                    self.hi_filename.replace('.fits', '_bin.fits')
            self.hi_error_bin_filename = \
                    self.hi_error_filename.replace('.fits', '_bin.fits')

            # Remove binned data
            _check_file(self.av_bin_filename,
                        clobber=True)
            _check_file(self.av_error_bin_filename,
                        clobber=True)
            _check_file(self.hi_bin_filename,
                        clobber=True)
            _check_file(self.hi_error_bin_filename,
                        clobber=True)
        else:
            self.av_data_bin = self.av_data
            self.av_error_data_bin = self.av_error_data
            self.hi_data_bin = self.hi_data

        self.nhi_image_bin = \
                calculate_nhi(cube=self.hi_data_bin,
                              velocity_axis=self.hi_vel_axis,
                              velocity_range=vel_range,
                              )

        if 0:
            import matplotlib.pyplot as plt
            plt.clf(); plt.close()
            plt.imshow(self.nhi_image_bin, origin='lower')
            plt.colorbar()
            plt.title(self.region)
            plt.savefig('/usr/users/ezbc/Desktop/' + self.region + '.png')

        if 'av_bin_map_filename_base' in self.plot_args:
            filename = self.plot_args['av_bin_map_filename_base'] + \
                       str(self.iter_step) + '.png'
            plot_av_bin_map(self.av_data,
                            self.av_data_bin,
                            av_header=self.av_header,
                            av_header_bin=self.av_header_bin,
                            filename=filename)

        # Calculate MLE for parameters with unbinned data
        # -----------------------------------------------
        #self._derive_bin_mask()
        if self.verbose:
            print('\n\n\tCalculating MLE parameters with initial errors...')

        results = \
            _calc_likelihoods(
                              hi_cube=self.hi_data_bin,
                              av_image=self.av_data_bin,
                              av_image_error=self.av_error_data_bin,
                              hi_vel_axis=self.hi_vel_axis,
                              width_grid=self.width_grid,
                              dgr_grid=self.dgr_grid,
                              intercept_grid=self.intercept_grid,
                              vel_center=self.vel_center,
                              bin_weights=self.bin_weights,
                              use_bin_weights=self.use_bin_weights,
                              likelihood_filename=self.likelihood_filename,
                              clobber=self.clobber_likelihoods,
                              verbose=self.verbose,
                              )

        self.iter_vars[self.iter_step]['init_likelihood_results'] = results

        if 1:
            self.iter_vars[self.iter_step]['scaled_likelihood_results'] = results
            likelihood_filename_base = \
                    '/d/bip3/ezbc/perseus/figures/likelihood/' + \
                    'perseus_likelihood_lee12_init_'
            self._write_final_params()

            plot_likelihoods_hist(cloud=self,
                          plot_axes=('widths', 'dgrs'),
                          show=0,
                          returnimage=False,
                          filename=likelihood_filename_base + 'wd.png',
                          limits=[0, 60, 0.0, 0.2],
                          )
            plot_likelihoods_hist(cloud=self,
                          plot_axes=('widths', 'intercepts'),
                          show=0,
                          returnimage=False,
                          filename=likelihood_filename_base + 'wi.png',
                          limits=[0, 60, -1.0, 2.0],
                          )

        # Unpack output of likelihood calculation
        dgr_max = results['dgr_max']
        intercept_max = results['intercept_max']
        vel_range_max = results['vel_range_max']

        # Rerun analysis with new error calculated Error should be calculated
        # across entire image, not just atomic regions, in order to understand
        # variation in
        # DGR
        # ---------------------------------------------------------------------
        # Calculate new standard deviation, set global variable npix - 2 is the
        # number of degrees of freedom see equation 15.1.6 in Numerical Recipes
        self._calc_model_error(vel_range_max, dgr_max, intercept_max)

        if self.verbose:
            print('\n\n\tCalculating MLE parameters with revised errors...')

        results = \
            _calc_likelihoods(
                              hi_cube=self.hi_data_bin,
                              av_image=self.av_data_bin,
                              av_image_error=self.av_error_data_bin,
                              hi_vel_axis=self.hi_vel_axis,
                              vel_center=self.vel_center,
                              bin_weights=self.bin_weights,
                              use_bin_weights=self.use_bin_weights,
                              width_grid=self.width_grid,
                              dgr_grid=self.dgr_grid,
                              intercept_grid=self.intercept_grid,
                              likelihood_filename=self.likelihood_filename,
                              clobber=self.clobber_likelihoods,
                              verbose=self.verbose,
                              )

        self.iter_vars[self.iter_step]['scaled_likelihood_results'] = results

        return results['vel_range_max']

    def _get_faint_mask(self, init_fraction=0.1):

        self._av_data_sorted = np.sort(self.av_data, axis=None)
        self._mask_pixel_fraction = init_fraction

        self._mask_pixel_number = int(self.av_data.size * init_fraction)

        self._mask_threshold = self._av_data_sorted[self._mask_pixel_number]

        mask = (self.av_data > self._mask_threshold)

        return mask

    def _add_pixels_to_mask(self, additional_fraction=0.01):

        self._mask_pixel_fraction += additional_fraction

        self._mask_pixel_number = \
                int(self.av_data.size * self._mask_pixel_fraction)

        self._mask_threshold = self._av_data_sorted[self._mask_pixel_number]

        mask = self.av_data > self._mask_threshold

        return mask

    def _trim_mask(self, mask_to_trim, trim_mask, abs_mask=None):

        mask = mask_to_trim + trim_mask

        mask[trim_mask == 0] = 0

        if abs_mask is not None:
            mask += abs_mask

        return mask

    def _combine_masks(self, mask_parent, mask_child, child_action='add'):

        ''' Mask will be wherever the parent is masked. The child will mask
         unmasked elements from the parent

        Parameters
        ----------
        child_action : str
            if 'add', then child is added to parent
            if 'subtract' then child is subtracted from parent
        '''

        mask = np.copy(mask_parent)

        if child_action == 'add':
            mask[mask_child] = 1
        elif child_action == 'subtract':
            mask[~mask_child] = 0

        return mask

    def _iterate_residual_masking__(self,):

        '''

        Returns
        -------
        av_model : numpy array
        mask : numpy array
        dgr : float

        '''

        import numpy as np

        # Mask out nans
        mask = self._prep_mask()

        # Apply initial mask to exclude throughout process
        mask_init = self.region_mask
        if mask_init is not None:
            mask += mask_init

        mask_faint = self._get_faint_mask()

        mask = self._combine_masks(mask_init, mask_faint)

        mask_residuals = np.zeros(mask.shape, dtype=bool)

        if 1:
            import matplotlib.pyplot as plt
            plt.clf(); plt.close()
            av_image = self.av_data.copy()
            av_image[mask] = np.nan
            plt.imshow(av_image, origin="lower", aspect="auto")
            plt.savefig('/usr/users/ezbc/Desktop/av_init.png')

        # solve for DGR using linear least squares
        if self.verbose:
            print('\n\tBeginning iterative DGR calculations + masking...')
            print('\n\tUsing faint masking method')

        # Iterate masking pixels which are correlated and rederiving a linear
        # least squares solution for the DGR
        #----------------------------------------------------------------------
        masking_results = {}
        delta_dgr = 1e10
        dgr = 1e10
        iteration = 0

        while delta_dgr > self.THRESHOLD_DELTA_DGR:


            if 0:
                import matplotlib.pyplot as plt
                plt.clf(); plt.close()
                av_image = self.av_data.copy()
                av_image[mask] = np.nan
                plt.imshow(av_image, origin="lower", aspect="auto")
                plt.savefig('/usr/users/ezbc/Desktop/av_iter' + \
                            str(self.iter_step) + '_' + \
                            str(iteration) + '.png')
                plt.show()

            if self.verbose:
                print('\t\tIteration {0:.0f} results:'.format(iteration))


            if 0:
                results = \
                    _calc_likelihoods(
                                      nhi_image=self.nhi_image[~mask],
                                      av_image=self.av_data[~mask],
                                      av_image_error=self.av_error_data[~mask],
                                      #bin_weights=bin_weights[~mask],
                                      #vel_center=vel_center_masked,
                                      #width_grid=np.arange(0,1,1),
                                      dgr_grid=self.dgr_grid,
                                      intercept_grid=self.intercept_grid,
                                      vel_center=self.vel_center,
                                      #results_filename='',
                                      #return_likelihoods=True,
                                      likelihood_filename=\
                                              self.likelihood_filename,
                                      clobber=self.clobber_likelihoods,
                                      verbose=self.verbose,
                                      )
            else:
                results = _fit_params(
                                      nhi_image=self.nhi_image[~mask],
                                      av_image=self.av_data[~mask],
                                      av_image_error=self.av_error_data[~mask],
                                      verbose=self.verbose,
                                      )

            # add new pixels
            mask_faint = \
                self._add_pixels_to_mask(additional_fraction=\
                                            self.PIXEL_MASK_INCREASE_FRACTION)
            mask = self._combine_masks(mask, mask_faint,
                                       child_action='subtract')
            mask = self._combine_masks(mask_init, mask)

            # Create model with the DGR
            dgr_new = results['dgr_max']
            intercept = results['intercept_max']
            av_image_model = self.nhi_image * dgr_new + intercept
            residuals = self.av_data - av_image_model
            residuals[mask] = np.nan

            # plot progress
            if 0:
                if plot_args is not None:
                    plot_residual_map(residuals,
                                      header=av_header,
                                      dgr=dgr_new)

            if 0:
                if iteration == 0:
                    plot_filename = self.plot_args['results_filename']
                else:
                    plot_filename = None

            # Include only residuals which are white noise
            self.plot_args['iter_ext'] = str(self.iter_step) + '_' + \
                                         str(iteration)
            self.plot_args['av_header'] = self.av_header
            mask_residuals_new, residual_threshold, fit_params = \
                _get_residual_mask(residuals,
                               residual_width_scale=self.RESIDUAL_WIDTH_SCALE,
                               plot_args=self.plot_args
                               )

            #intercept_grid = np.linspace(intercept, intercept + 1.0, 1.0)

            # Mask non-white noise, i.e. correlated residuals.
            mask_residuals = \
                    self._combine_masks(mask_residuals, mask_residuals_new)

            if 0:
                import matplotlib.pyplot as plt
                plt.clf(); plt.close()
                av_image = self.av_data.copy()
                av_image[mask] = np.nan
                plt.imshow(av_image, origin="lower", aspect="auto")
                plt.savefig('/usr/users/ezbc/Desktop/av_residmask_iter' + \
                            str(self.iter_step) + '_' + \
                            str(iteration) + '.png')
                plt.show()

            # Record variables
            masking_results[iteration] = {}
            masking_results[iteration]['residuals'] = residuals
            masking_results[iteration]['dgr'] = dgr_new
            #masking_results[iteration]['width'] = width
            masking_results[iteration]['intercept'] = intercept
            masking_results[iteration]['fit_params'] = fit_params
            masking_results[iteration]['residual_threshold'] = \
                    residual_threshold
            masking_results[iteration]['mask'] = mask

            #if dgr == 1e10:
            #    mask[self.av_data > 1] = 1

            # Reset while loop conditions
            #delta_dgr = np.abs(1 - dgr / dgr_new)
            delta_dgr = np.abs(dgr - dgr_new)
            dgr = dgr_new
            iteration += 1

            # Derive new mask
            mask = self._combine_masks(mask_residuals, mask)

            if self.verbose:
                npix = mask.size - np.sum(mask)
                print('\t\t\tNumber of non-masked pixels = {0:.0f}'.format(npix))

        # Plot results
        if 0:
            mask_new = get_residual_mask(residuals,
                                      residual_width_scale=self.RESIDUAL_WIDTH_SCALE,
                                      plot_progress=plot_progress,
                                      results_filename=results_filename)

        # Create model of Av
        self.av_model = dgr * self.nhi_image + intercept
        self.mask = mask
        self.av_model[self.mask] = np.nan
        self.iter_vars[self.iter_step]['masking_var'] = masking_results
        self.iter_vars[self.iter_step]['mask'] = mask

    def _iterate_residual_masking(self,):

        '''

        Returns
        -------
        av_model : numpy array
        mask : numpy array
        dgr : float

        '''

        import numpy as np

        # Mask out nans
        mask = self._prep_mask()

        # Apply initial mask to exclude throughout process
        mask_init = self.region_mask
        if mask_init is not None:
            mask += mask_init

        mask_faint = self._get_faint_mask()

        mask = self._combine_masks(mask_init, mask_faint)

        mask_residuals = np.zeros(mask.shape, dtype=bool)

        if 1:
            import matplotlib.pyplot as plt
            plt.clf(); plt.close()
            av_image = self.av_data.copy()
            av_image[mask] = np.nan
            plt.imshow(av_image, origin="lower", aspect="auto")
            plt.savefig('/usr/users/ezbc/Desktop/av_init.png')

        # solve for DGR using linear least squares
        if self.verbose:
            print('\n\tBeginning iterative DGR calculations + masking...')
            print('\n\tUsing faint masking method')

        # Iterate masking pixels which are correlated and rederiving a linear
        # least squares solution for the DGR
        #----------------------------------------------------------------------
        masking_results = {}
        delta_dgr = 1e10
        dgr = 1e10
        iteration = 0

        while delta_dgr > self.THRESHOLD_DELTA_DGR:


            if 0:
                import matplotlib.pyplot as plt
                plt.clf(); plt.close()
                av_image = self.av_data.copy()
                av_image[mask] = np.nan
                plt.imshow(av_image, origin="lower", aspect="auto")
                plt.savefig('/usr/users/ezbc/Desktop/av_iter' + \
                            str(self.iter_step) + '_' + \
                            str(iteration) + '.png')
                plt.show()

            if self.verbose:
                print('\t\tIteration {0:.0f} results:'.format(iteration))


            if 0:
                results = \
                    _calc_likelihoods(
                                      nhi_image=self.nhi_image[~mask],
                                      av_image=self.av_data[~mask],
                                      av_image_error=self.av_error_data[~mask],
                                      #bin_weights=bin_weights[~mask],
                                      #vel_center=vel_center_masked,
                                      #width_grid=np.arange(0,1,1),
                                      dgr_grid=self.dgr_grid,
                                      intercept_grid=self.intercept_grid,
                                      vel_center=self.vel_center,
                                      #results_filename='',
                                      #return_likelihoods=True,
                                      likelihood_filename=\
                                              self.likelihood_filename,
                                      clobber=self.clobber_likelihoods,
                                      verbose=self.verbose,
                                      )
            else:
                results = _fit_params(
                                      nhi_image=self.nhi_image[~mask],
                                      av_image=self.av_data[~mask],
                                      av_image_error=self.av_error_data[~mask],
                                      verbose=self.verbose,
                                      )

            # add new pixels
            mask_faint = \
                self._add_pixels_to_mask(additional_fraction=\
                                            self.PIXEL_MASK_INCREASE_FRACTION)
            mask = self._combine_masks(mask, mask_faint,
                                       child_action='subtract')
            mask = self._combine_masks(mask_init, mask)

            # Create model with the DGR
            dgr_new = results['dgr_max']
            intercept = results['intercept_max']
            av_image_model = self.nhi_image * dgr_new + intercept
            residuals = self.av_data - av_image_model
            residuals[mask] = np.nan

            # plot progress
            if 0:
                if plot_args is not None:
                    plot_residual_map(residuals,
                                      header=av_header,
                                      dgr=dgr_new)

            if 0:
                if iteration == 0:
                    plot_filename = self.plot_args['results_filename']
                else:
                    plot_filename = None

            # Include only residuals which are white noise
            self.plot_args['iter_ext'] = str(self.iter_step) + '_' + \
                                         str(iteration)
            self.plot_args['av_header'] = self.av_header
            mask_residuals_new, residual_threshold = \
                _get_residual_mask(residuals,
                               residual_width_scale=self.RESIDUAL_WIDTH_SCALE,
                               plot_args=self.plot_args
                               )

            #intercept_grid = np.linspace(intercept, intercept + 1.0, 1.0)

            # Mask non-white noise, i.e. correlated residuals.
            mask_residuals = \
                    self._combine_masks(mask_residuals, mask_residuals_new)

            if 0:
                import matplotlib.pyplot as plt
                plt.clf(); plt.close()
                av_image = self.av_data.copy()
                av_image[mask] = np.nan
                plt.imshow(av_image, origin="lower", aspect="auto")
                plt.savefig('/usr/users/ezbc/Desktop/av_residmask_iter' + \
                            str(self.iter_step) + '_' + \
                            str(iteration) + '.png')
                plt.show()

            # Record variables
            masking_results[iteration] = {}
            masking_results[iteration]['residuals'] = residuals
            masking_results[iteration]['dgr'] = dgr_new
            #masking_results[iteration]['width'] = width
            masking_results[iteration]['intercept'] = intercept
            masking_results[iteration]['mask'] = mask
            masking_results[iteration]['fit_params'] = fit_params
            masking_results[iteration]['residual_threshold'] = \
                    residual_threshold

            #if dgr == 1e10:
            #    mask[self.av_data > 1] = 1

            # Reset while loop conditions
            delta_dgr = np.abs(dgr - dgr_new)
            dgr = dgr_new
            iteration += 1

            # Derive new mask
            mask = self._combine_masks(mask_residuals, mask)

            if self.verbose:
                npix = mask.size - np.sum(mask)
                print('\t\t\tNumber of non-masked pixels = {0:.0f}'.format(npix))

        # Plot results
        if 0:
            mask_new = get_residual_mask(residuals,
                                      residual_width_scale=self.RESIDUAL_WIDTH_SCALE,
                                      plot_progress=plot_progress,
                                      results_filename=results_filename)

        # Create model of Av
        self.av_model = dgr * self.nhi_image + intercept
        self.mask = mask
        self.av_model[self.mask] = np.nan
        self.iter_vars[self.iter_step]['masking_var'] = masking_results
        self.iter_vars[self.iter_step]['mask'] = mask



    def _write_final_params(self,):

        last_iter = max(self.iter_vars.keys())

        props = self.iter_vars[last_iter]['scaled_likelihood_results']

        # Write results to global properties
        self.props['dust2gas_ratio'] = {}
        self.props['dust2gas_ratio_error'] = {}
        self.props['intercept'] = {}
        self.props['intercept_error'] = {}
        self.props['hi_velocity_width'] = {}
        self.props['hi_velocity_width_error'] = {}
        self.props['dust2gas_ratio_max'] = {}
        self.props['intercept_max'] = {}
        self.props['hi_velocity_center'] = {}
        self.props['hi_velocity_center_bin'] = {}
        self.props['hi_velocity_width_max'] = {}
        self.props['hi_velocity_range_max'] =  {}
        self.props['av_threshold'] = {}
        self.props['co_threshold'] = {}
        self.props['hi_velocity_width_max']['value'] = \
                props['width_max']
        self.props['hi_velocity_width_max']['unit'] = 'km/s'
        self.props['hi_velocity_width']['value'] = \
                props['width_confint'][0]
        self.props['hi_velocity_width']['unit'] = 'km/s'
        self.props['hi_velocity_width_error']['value'] = \
                props['width_confint'][1:]
        self.props['hi_velocity_width_error']['unit'] = 'km/s'
        self.props['hi_velocity_range'] = \
                props['vel_range_confint'][0:2]
        self.props['hi_velocity_range_error'] = \
                props['vel_range_confint'][2:]
        self.props['dust2gas_ratio']['value'] = \
                props['dgr_confint'][0]
        self.props['dust2gas_ratio_error']['value'] = \
                props['dgr_confint'][1:]
        self.props['dust2gas_ratio_max']['value'] = \
                props['dgr_max']
        self.props['intercept_max']['value'] = \
                props['intercept_max']
        self.props['intercept']['value'] = \
                props['intercept_confint'][0]
        self.props['intercept_error']['value'] = \
                props['intercept_confint'][1:]
        self.props['hi_velocity_center']['value'] = \
                self.vel_center
        self.props['hi_velocity_center']['unit'] = 'km/s'
        self.props['single_vel_center'] = self.vel_center
        #self.props['hi_velocity_width_max']['value'] = width_max
        #self.props['hi_velocity_range_max']['value'] = vel_range_max
        #self.props['hi_velocity_range_conf'] = self.conf
        self.props['width_likelihood'] = \
                props['width_likelihood']
        self.props['dgr_likelihood'] = \
                props['dgr_likelihood']
        self.props['vel_centers'] = self.vel_center
        self.props['width_grid'] = self.width_grid
        self.props['dgr_grid'] = self.dgr_grid
        self.props['intercept_grid'] = self.intercept_grid
        self.props['likelihoods'] = props['likelihoods']

    def _plot_likelihoods(self,):

        # Plot the likelihoods
        if 'likelihood_filename_base' in self.plot_args:

            import os

            filename_wd = self.plot_args['likelihood_filename_base'] + \
                          self.plot_args['iter_ext'] + '_wd.png'
            filename_wi = self.plot_args['likelihood_filename_base'] + \
                          self.plot_args['iter_ext'] + '_wi.png'


            if self.iter_step == 0:
                os.system('rm -rf ' + \
                        self.plot_args['likelihood_filename_base'] + '*')

            plot_likelihoods_hist(cloud=self,
                                  filename=filename_wd,
                                  plot_axes=('widths', 'dgrs'))
            plot_likelihoods_hist(cloud=self,
                                  filename=filename_wi,
                                  plot_axes=('widths', 'intercepts'))

    def load_cloud_properties(self, prop_filename):

        import json

        # Load global properties
        with open(prop_filename, 'r') as f:
            self.props = json.load(f)

    def load_region(self, region_filename, header=None):

        import pyregion as pyr

        # region[0] in following format:
        # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
        # [ra center, dec center, width, height, rotation angle]

        regions = pyr.open(region_filename)

        self.props['regions'] = {}

        for region in regions:
            # Cores defined in following format: 'tag={L1495A}'
            tag = region.comment
            region_name = tag[tag.find('text={')+6:tag.find('}')].lower()

            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            self.props['regions'][region_name] = {}
            self.props['regions'][region_name]['poly_verts'] = {}
            self.props['regions'][region_name]['poly_verts']['wcs'] = poly_verts

            if self.av_header is not None or header is not None:
                poly_verts_pix = []
                for i in xrange(0, len(poly_verts)):
                    poly_verts_pix.append(\
                            self._get_pix_coords(ra=poly_verts[i][0],
                                                 dec=poly_verts[i][1],
                                                 header=self.av_header)[:-1][::-1].tolist())

                self.props['regions'][region_name]['poly_verts']['pixel'] = \
                    poly_verts_pix

    def run_analysis(self, region_filename=None, region='california'):

        self.iter_step = 0
        self.region = region

        # Change WCS coords to pixel coords of images
        self._convert_coordinates()

        # Load cloud division regions from ds9
        self.load_region(region_filename)

        # Initialize
        vel_range = self.init_vel_range
        vel_range_new = (-1.0, 1.0)
        vel_range_diff = np.sum(np.abs(np.array(vel_range) - \
                                       np.array(vel_range_new)))


        while vel_range_diff > self.VEL_RANGE_DIFF_THRES:
            if self.verbose:
                print('\nIteration ' + str(self.iter_step))

            vel_range_new = \
                    self._iterate_mle_calc(vel_range=vel_range,
                                           )

            # Write parameters
            self._write_final_params()

            # Check to see if likelihoods should be plotted
            self._plot_likelihoods()

            vel_range_diff = np.sum(np.abs(np.array(vel_range) - \
                                           np.array(vel_range_new)))

            if self.verbose:
                print('\t\tVel range old = ' + \
                      '{0:.1f} to {1:.1f}'.format(*vel_range))
                print('\t\tVel range new = ' + \
                      '{0:.1f} to {1:.1f}'.format(*vel_range_new))

            vel_range = vel_range_new

            self.iter_step += 1

        # Write prints to file?
        if self._diagnostics is not None:
            import sys

            sys.stdout = self._orig_stdout

            self._diagnostics.close()

            # Set to none so class can be pickled
            self._diagnostics = None
            self._orig_stdout = None

class KeyboardInterruptError(Exception): pass

'''
Module Functions

'''

def _calc_likelihoods__(
        hi_cube=None,
        hi_vel_axis=None,
        nhi_image=None,
        av_image=None,
        av_image_error=None,
        bin_weights=None,
        use_bin_weights=False,
        vel_center=None,
        width_grid=None,
        dgr_grid=None,
        intercept_grid=None,
        plot_results=False,
        results_filename='',
        return_likelihoods=True,
        likelihood_filename=None,
        clobber=False,
        conf=0.68,
        threshold_delta_dgr=0.0005,
        verbose=False,
        ):

    '''
    Parameters
    ----------

    Returns
    -------
    vel_range : tuple
        Lower and upper bound of HI velocity range in km/s which provides the
        best likelihoodelated N(HI) distribution with Av.
    likelihoods : array-like, optional
        Array of Pearson likelihoodelation coefficients likelihoodesponding to each
        permutation through the velocity centers and velocity widths.

    '''

    import numpy as np
    from myimage_analysis import calculate_nhi
    from os import path
    from astropy.io import fits
    from mystats import calc_symmetric_error, calc_logL
    import multiprocessing as mp

    # Check if likelihood grid should be derived
    if likelihood_filename is not None:
        if not path.isfile(likelihood_filename):
            perform_mle = True
            write_mle = True
        elif clobber:
            perform_mle = True
            write_mle = True
        else:
            perform_mle = False
            write_mle = False
    # If no filename provided, do not read file and do not write file
    else:
        write_mle = False
        perform_mle = True

    if perform_mle:
        # calculate the likelihoodelation coefficient for each velocity
        # range
        if width_grid is None:
            width_grid = np.linspace(0, 1, 1)
        if dgr_grid is None:
            dgr_grid = np.linspace(0, 1, 1)
        if intercept_grid is None:
            intercept_grid = np.linspace(0, 1, 1)

        likelihoods = np.empty((len(width_grid),
                                len(dgr_grid),
                                len(intercept_grid)))

        # Progress bar parameters
        total = float(likelihoods.size)
        count = 0

        if 0:
            import matplotlib.pyplot as plt
            plt.clf(); plt.close()
            plt.imshow(av_image_error,
                       origin='lower', aspect='auto')
            plt.colorbar()
            plt.savefig('/usr/users/ezbc/Desktop/resid.png')


        if nhi_image is not None:
            # Remove nans
            mask = ((av_image != av_image) | \
                    (av_image_error != av_image_error) | \
                    (av_image_error == 0) | \
                    (nhi_image != nhi_image))

            av_image = av_image[~mask]
            av_image_error = av_image_error[~mask]
            nhi_image = nhi_image[~mask]
            if use_bin_weights and bin_weights is not None:
                bin_weights = bin_weights[~mask]

            print_vel_range = False

            # Cycle through DGR to estimate error
            width_grid = np.linspace(0, 1, 1)

            for j, vel_width in enumerate(width_grid):
                for k, dgr in enumerate(dgr_grid):
                    for m, intercept in enumerate(intercept_grid):
                        # Create model of Av with N(HI) and DGR
                        av_image_model = nhi_image * dgr + intercept

                        logL = calc_logL(av_image_model,
                                         av_image,
                                         data_error=av_image_error,
                                         )

                        #print logL
                        likelihoods[j, k, m] = logL
        else:
            print_vel_range = True

            # Remove nans
            if not use_bin_weights:
                hi_cube[hi_cube != hi_cube] = 0.0
                mask = ((av_image != av_image) | \
                        (av_image_error != av_image_error) | \
                        (av_image_error == 0))

                av_image = av_image[~mask]
                av_image_error = av_image_error[~mask]
                hi_cube = hi_cube[:, ~mask]

            # Weight images
            #bin_weights = None
            if use_bin_weights and bin_weights is not None:

                av_image, av_image_error, hi_cube = \
                    _prep_weighted_data(av_image=av_image,
                                        av_image_error=av_image_error,
                                        hi_cube=hi_cube,
                                        weights=bin_weights)

            for j, vel_width in enumerate(width_grid):
                # Construct N(HI) image outside of DGR loop, then apply
                # dgr_grid in loop

                vel_range = np.array((vel_center - vel_width / 2.,
                                      vel_center + vel_width / 2.))

                nhi_image = calculate_nhi(cube=hi_cube,
                                          velocity_axis=hi_vel_axis,
                                          velocity_range=vel_range,
                                          return_nhi_error=False)

                # Cycle through DGR to estimate error
                for k, dgr in enumerate(dgr_grid):
                    for m, intercept in enumerate(intercept_grid):
                        # Create model of Av with N(HI) and DGR
                        av_image_model = nhi_image * dgr + intercept

                        logL = calc_logL(av_image_model,
                                         av_image,
                                         data_error=av_image_error,
                                         )

                        #print logL

                        likelihoods[j, k, m] = logL

    # Load file of likelihoods
    elif not perform_mle:
        print('Reading likelihood grid file:')
        print(likelihood_filename)

        hdu = fits.open(likelihood_filename)
        likelihoods = hdu[0].data

        if len(width_grid) != likelihoods.shape[0] or \
           len(dgr_grid) != likelihoods.shape[1]:
            raise ValueError('Specified parameter grid not the same as in' + \
                    'loaded data likelihoods.')

        likelihoods = np.ma.array(likelihoods,
                mask=(likelihoods != likelihoods))

    # Normalize the log likelihoods
    likelihoods -= np.nanmax(likelihoods)

    # Convert to likelihoods
    likelihoods = np.exp(likelihoods)

    likelihoods[np.isnan(likelihoods)] = 0.0

    # Normalize the likelihoods
    likelihoods = likelihoods / np.nansum(likelihoods)

    #print width_grid
    #print dgr_grid
    #print intercept_grid
    #print likelihoods

    # Derive marginal distributions of both centers and widths
    intercept_likelihood = np.sum(likelihoods, axis=(0, 1)) / \
                                  np.sum(likelihoods)

    width_likelihood = np.sum(likelihoods, axis=(1, 2)) / \
            np.sum(likelihoods)

    dgr_likelihood = np.sum(likelihoods, axis=(0, 2)) / \
            np.sum(likelihoods)

    # Derive confidence intervals of parameters
    width_confint = calc_symmetric_error(width_grid,
                                   width_likelihood,
                                   alpha=1.0 - conf)
    dgr_confint = calc_symmetric_error(dgr_grid,
                                 dgr_likelihood,
                                 alpha=1.0 - conf)
    intercept_confint = calc_symmetric_error(intercept_grid,
                                 intercept_likelihood,
                                 alpha=1.0 - conf)

    # Get values of best-fit model parameters
    max_loc = np.where(likelihoods == np.max(likelihoods))

    width_max = width_grid[max_loc[0][0]]
    dgr_max = dgr_grid[max_loc[1][0]]
    intercept_max = intercept_grid[max_loc[2][0]]

    if 0:
        import matplotlib.pyplot as plt
        plt.imshow(likelihoods[:, :, max_loc[2][0]], origin='lower', aspect='auto',
                   extent=[width_grid[0], width_grid[-1], dgr_grid[0],
                           dgr_grid[-1]])
        plt.colorbar()
        plt.show()
        plt.imshow(likelihoods[:, 0, :], origin='lower', aspect='auto',
                   extent=[width_grid[0], width_grid[-1], intercept_grid[0],
                           intercept_grid[-1]])
        plt.colorbar()
        plt.show()

    if verbose:
        if print_vel_range:
            print('\n\t\t\tWidth confint = ' + \
                    '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                            width_confint[2],
                                                   np.abs(width_confint[1])))
        print('\n\t\t\tDGR confint = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} 10^-20 cm^2 mag'.format(dgr_confint[0],
                                                    dgr_confint[2],
                                                    np.abs(dgr_confint[1])))
        print('\n\t\t\tIntercept confint = ' + \
        '{0:.2f} +{1:.2f}/-{2:.2f} mag'.format(intercept_confint[0],
                                                intercept_confint[2],
                                                np.abs(intercept_confint[1])))

    # Write PDF
    if vel_center is None:
        vel_center = 0.0

    upper_lim = (np.nanmean(vel_center) + width_confint[0]/2.)
    lower_lim = (np.nanmean(vel_center) - width_confint[0]/2.)
    upper_lim_error = width_confint[2]**2
    lower_lim_error = width_confint[1]**2

    vel_range_confint = (lower_lim, upper_lim, lower_lim_error,
                         upper_lim_error)
    vel_range_max = (vel_center - width_max/2.0, vel_center + width_max/2.0)

    results = {
               'vel_range_confint': vel_range_confint,
               'width_confint': width_confint,
               'dgr_confint': dgr_confint,
               'intercept_confint': intercept_confint,
               'likelihoods': likelihoods,
               'width_likelihood': width_likelihood,
               'dgr_likelihood': dgr_likelihood,
               'intercept_likelihood': intercept_likelihood,
               'width_max': width_max,
               'dgr_max': dgr_max,
               'intercept_max': intercept_max,
               'vel_range_max': vel_range_max
               }

    if not return_likelihoods:
        return vel_range_confint, dgr_confint
    else:
        return results

def _calc_likelihoods(
        hi_cube=None,
        hi_vel_axis=None,
        nhi_image=None,
        av_image=None,
        av_image_error=None,
        bin_weights=None,
        use_bin_weights=False,
        vel_center=None,
        width_grid=None,
        dgr_grid=None,
        intercept_grid=None,
        plot_results=False,
        results_filename='',
        return_likelihoods=True,
        likelihood_filename=None,
        clobber=False,
        conf=0.68,
        threshold_delta_dgr=0.0005,
        verbose=False,
        ):

    '''
    Parameters
    ----------

    Returns
    -------
    vel_range : tuple
        Lower and upper bound of HI velocity range in km/s which provides the
        best likelihoodelated N(HI) distribution with Av.
    likelihoods : array-like, optional
        Array of Pearson likelihoodelation coefficients likelihoodesponding to each
        permutation through the velocity centers and velocity widths.

    '''

    import numpy as np
    from myimage_analysis import calculate_nhi
    from os import path
    from astropy.io import fits
    from mystats import calc_symmetric_error, calc_logL
    import multiprocessing as mp

    # Check if likelihood grid should be derived
    if likelihood_filename is not None:
        if not path.isfile(likelihood_filename):
            perform_mle = True
            write_mle = True
        elif clobber:
            perform_mle = True
            write_mle = True
        else:
            perform_mle = False
            write_mle = False
    # If no filename provided, do not read file and do not write file
    else:
        write_mle = False
        perform_mle = True

    if perform_mle:
        # calculate the likelihoodelation coefficient for each velocity
        # range
        if width_grid is None:
            width_grid = np.linspace(0, 1, 1)
        if dgr_grid is None:
            dgr_grid = np.linspace(0, 1, 1)
        if intercept_grid is None:
            intercept_grid = np.linspace(0, 1, 1)

        likelihoods = np.empty((len(width_grid),
                                len(dgr_grid),
                                len(intercept_grid)))

        # Progress bar parameters
        total = float(likelihoods.size)
        count = 0

        if 0:
            import matplotlib.pyplot as plt
            plt.clf(); plt.close()
            plt.imshow(av_image_error,
                       origin='lower', aspect='auto')
            plt.colorbar()
            plt.savefig('/usr/users/ezbc/Desktop/resid.png')


        if nhi_image is not None:
            # Remove nans
            mask = ((av_image != av_image) | \
                    (av_image_error != av_image_error) | \
                    (av_image_error == 0) | \
                    (nhi_image != nhi_image))

            av_image = av_image[~mask]
            av_image_error = av_image_error[~mask]
            nhi_image = nhi_image[~mask]
            if use_bin_weights and bin_weights is not None:
                bin_weights = bin_weights[~mask]

            print_vel_range = False

            # Cycle through DGR to estimate error
            width_grid = np.linspace(0, 1, 1)

            for j, vel_width in enumerate(width_grid):
                for k, dgr in enumerate(dgr_grid):
                    for m, intercept in enumerate(intercept_grid):
                        # Create model of Av with N(HI) and DGR
                        av_image_model = nhi_image * dgr + intercept

                        logL = calc_logL(av_image_model,
                                         av_image,
                                         data_error=av_image_error,
                                         )

                        #print logL
                        likelihoods[j, k, m] = logL
        else:
            print_vel_range = True

            # Remove nans
            if not use_bin_weights:
                hi_cube[hi_cube != hi_cube] = 0.0
                mask = ((av_image != av_image) | \
                        (av_image_error != av_image_error) | \
                        (av_image_error == 0))

                av_image = av_image[~mask]
                av_image_error = av_image_error[~mask]
                hi_cube = hi_cube[:, ~mask]

            # Weight images
            #bin_weights = None
            if use_bin_weights and bin_weights is not None:

                av_image, av_image_error, hi_cube = \
                    _prep_weighted_data(av_image=av_image,
                                        av_image_error=av_image_error,
                                        hi_cube=hi_cube,
                                        weights=bin_weights)

            def worker(args):

                j = args['j']
                vel_width = args['vel_width']
                likelihoods = np.empty((len(dgr_grid), len(intercept_grid)))
                output = args['output']

                try:
                    vel_range = np.array((vel_center - vel_width / 2.,
                                          vel_center + vel_width / 2.))

                    nhi_image = calculate_nhi(cube=hi_cube,
                                              velocity_axis=hi_vel_axis,
                                              velocity_range=vel_range,
                                              return_nhi_error=False)

                    # Cycle through DGR to estimate error
                    for k, dgr in enumerate(dgr_grid):
                        for m, intercept in enumerate(intercept_grid):
                            # Create model of Av with N(HI) and DGR
                            av_image_model = nhi_image * dgr + intercept

                            logL = calc_logL(av_image_model,
                                             av_image,
                                             data_error=av_image_error,
                                             )

                            likelihoods[k, m] = logL

                    output.put([j, likelihoods])

                except KeyboardInterrupt:
                    raise KeyboardInterruptError()

            processes = []
            # use all but one cpu
            output = mp.Queue(mp.cpu_count() - 1)

            args = {'likelihoods': likelihoods,
                    'output': output
                    }

            for j, vel_width in enumerate(width_grid):
                # Construct N(HI) image outside of DGR loop, then apply
                # dgr_grid in loop

                args['j'] = j
                args['vel_width'] = vel_width

                try:
                    processes.append(mp.Process(target=worker,
                                                args=(args,)))

                    processes[j].start()
                except KeyboardInterrupt():
                    p.terminate()

            for p in processes:
                result = output.get()
                likelihoods[result[0], :, :] = result[1]

    # Load file of likelihoods
    elif not perform_mle:
        print('Reading likelihood grid file:')
        print(likelihood_filename)

        hdu = fits.open(likelihood_filename)
        likelihoods = hdu[0].data

        if len(width_grid) != likelihoods.shape[0] or \
           len(dgr_grid) != likelihoods.shape[1]:
            raise ValueError('Specified parameter grid not the same as in' + \
                    'loaded data likelihoods.')

        likelihoods = np.ma.array(likelihoods,
                mask=(likelihoods != likelihoods))

    # Normalize the log likelihoods
    likelihoods -= np.nanmax(likelihoods)

    # Convert to likelihoods
    likelihoods = np.exp(likelihoods)

    likelihoods[np.isnan(likelihoods)] = 0.0

    # Normalize the likelihoods
    likelihoods = likelihoods / np.nansum(likelihoods)

    #print width_grid
    #print dgr_grid
    #print intercept_grid
    #print likelihoods

    # Derive marginal distributions of both centers and widths
    intercept_likelihood = np.sum(likelihoods, axis=(0, 1)) / \
                                  np.sum(likelihoods)

    width_likelihood = np.sum(likelihoods, axis=(1, 2)) / \
            np.sum(likelihoods)

    dgr_likelihood = np.sum(likelihoods, axis=(0, 2)) / \
            np.sum(likelihoods)

    # Derive confidence intervals of parameters
    width_confint = calc_symmetric_error(width_grid,
                                   width_likelihood,
                                   alpha=1.0 - conf)
    dgr_confint = calc_symmetric_error(dgr_grid,
                                 dgr_likelihood,
                                 alpha=1.0 - conf)
    intercept_confint = calc_symmetric_error(intercept_grid,
                                 intercept_likelihood,
                                 alpha=1.0 - conf)

    # Get values of best-fit model parameters
    max_loc = np.where(likelihoods == np.max(likelihoods))

    width_max = width_grid[max_loc[0][0]]
    dgr_max = dgr_grid[max_loc[1][0]]
    intercept_max = intercept_grid[max_loc[2][0]]

    if 0:
        import matplotlib.pyplot as plt
        plt.imshow(likelihoods[:, :, max_loc[2][0]], origin='lower', aspect='auto',
                   extent=[width_grid[0], width_grid[-1], dgr_grid[0],
                           dgr_grid[-1]])
        plt.colorbar()
        plt.show()
        plt.imshow(likelihoods[:, 0, :], origin='lower', aspect='auto',
                   extent=[width_grid[0], width_grid[-1], intercept_grid[0],
                           intercept_grid[-1]])
        plt.colorbar()
        plt.show()

    if verbose:
        if print_vel_range:
            print('\n\t\t\tWidth confint = ' + \
                    '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                            width_confint[2],
                                                   np.abs(width_confint[1])))
        print('\n\t\t\tDGR confint = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} 10^-20 cm^2 mag'.format(dgr_confint[0],
                                                    dgr_confint[2],
                                                    np.abs(dgr_confint[1])))
        print('\n\t\t\tIntercept confint = ' + \
        '{0:.2f} +{1:.2f}/-{2:.2f} mag'.format(intercept_confint[0],
                                                intercept_confint[2],
                                                np.abs(intercept_confint[1])))

    # Write PDF
    if vel_center is None:
        vel_center = 0.0

    upper_lim = (np.nanmean(vel_center) + width_confint[0]/2.)
    lower_lim = (np.nanmean(vel_center) - width_confint[0]/2.)
    upper_lim_error = width_confint[2]**2
    lower_lim_error = width_confint[1]**2

    vel_range_confint = (lower_lim, upper_lim, lower_lim_error,
                         upper_lim_error)
    vel_range_max = (vel_center - width_max/2.0, vel_center + width_max/2.0)

    results = {
               'vel_range_confint': vel_range_confint,
               'width_confint': width_confint,
               'dgr_confint': dgr_confint,
               'intercept_confint': intercept_confint,
               'likelihoods': likelihoods,
               'width_likelihood': width_likelihood,
               'dgr_likelihood': dgr_likelihood,
               'intercept_likelihood': intercept_likelihood,
               'width_max': width_max,
               'dgr_max': dgr_max,
               'intercept_max': intercept_max,
               'vel_range_max': vel_range_max
               }

    if not return_likelihoods:
        return vel_range_confint, dgr_confint
    else:
        return results

def _fit_params(nhi_image=None,
                av_image=None,
                av_image_error=None,
                verbose=False):

    mask = (np.isnan(nhi_image) | \
            np.isnan(av_image) | \
            np.isnan(av_image_error))

    b = av_image[~mask]
    A = np.array([nhi_image[~mask], np.ones(nhi_image[~mask].shape)]).T
    #A = np.array([nhi_image[~mask],]).T

    params = np.dot(np.linalg.pinv(A), b)

    if verbose:
        print('\t\t\tDGR = {0:.2f} [10^-20 cm^2 mag]'.format(params[0]))
        print('\t\t\tIntercept = {0:.2f} [mag]'.format(params[1]))

    results = {}
    results['dgr_max'] = params[0]
    results['intercept_max'] = params[1]
    #results['intercept_max'] = 0

    return results

def _prep_weighted_data(av_image=None,
                        av_image_error=None,
                        hi_cube=None,
                        weights=None):

        hi_cube[hi_cube != hi_cube] = 0.0
        mask = ((av_image != av_image) | \
                (av_image_error != av_image_error) | \
                (av_image_error == 0) | \
                (weights <= 0) | \
                (weights != weights))

        av_image = av_image[~mask]
        av_image_error = av_image_error[~mask]
        hi_cube = hi_cube[:, ~mask]
        weights = weights[~mask]

        # omit 0 weights
        data = av_image
        error = av_image_error
        model = hi_cube

        weights = np.array(weights / np.nanmin(weights), dtype=int)
        av_image = np.zeros(np.sum(weights))
        av_image_error = np.zeros(np.sum(weights))
        hi_cube = np.zeros((hi_cube.shape[0], np.sum(weights)))
        count = 0
        for i in xrange(0, len(weights)):
            av_image[count:count + weights[i]] = data[i]
            av_image_error[count:count + weights[i]] = \
                    error[i]
            for j in xrange(0, hi_cube.shape[0]):
                hi_cube[j, count:count + weights[i]] = model[j, i]
            count += weights[i]


        return (av_image, av_image_error, hi_cube)

def _get_residual_mask(residuals, residual_width_scale=3.0, plot_args={},
        use_GMM=False):

    '''

    '''

    #import numpy as np
    from mystats import gauss
    from scipy.optimize import curve_fit, minimize
    from sklearn.mixture import GMM, DPGMM

    # Fit the rising portion of the residuals
    residuals_crop = residuals[(residuals < 0) & \
                               #(residuals > -1.5) & \
                                ~np.isnan(residuals)]

    #if not use_GMM:
    if not use_GMM:
        #
        #print('histing')
        counts, bin_edges = np.histogram(np.ravel(residuals_crop),
                                         bins=100,
                                         )
        p0=(0.1, np.nanmax(counts), 0)

        from lmfit import minimize, Parameters

        # Set parameter limits and initial guesses
        params = Parameters()
        params.add('width',
                   value=0.1,
                   min=0.05,
                   max=10,
                   )
        params.add('amp',
                   value=np.nanmax(counts),
                   min=0,
                   max=2 * np.nanmax(counts),
                   )
        params.add('x0',
                   value=0,
                   min=-2,
                   max=2,
                   vary=False,
                   )

        #bin_edges = residuals_crop
        #counts = np.ones(residuals_crop.size - 1)

        def norm(params, bin_edges, counts):
            width = params['width'].value
            amp = params['amp'].value
            x0 = params['x0'].value

            model = gauss(bin_edges, width, amp, x0)

            norm = np.sum((counts - model)**2)

            return norm

        #print('fitting')
        # Perform the fit!
        result = minimize(norm,
                          params,
                          args=(bin_edges[:-1], counts),
                          method='nelder'
                          )

        #print('done fitting')
        fit_params = (params['width'].value, params['amp'].value,
                params['x0'].value)

    if use_GMM:
        g = DPGMM()
        g.fit(residuals_crop)
        fit_params = [g.covars_[0,0]**0.5, 1, g.means_[0,0]]
        #print("mean : %f, std : %f" % (g.means_[0, 0], g.covars_[0, 0]**0.5))

        x_fit = np.linspace(np.nanmin(residuals),
                            np.nanmax(residuals),
                            1000)
        #y_fit = np.exp(g.score_samples(x_fit)[0])
        y_fit = gauss(x_fit, *fit_params)

        #fit_params[0] = y_fit[y_fit == 2.35 * max(y_fit)]/2.0

        if 1:
            import matplotlib.pyplot as plt
            plt.clf(); plt.close()
            counts, bin_edges = np.histogram(residuals[~np.isnan(residuals)],
                                             bins=50,
                                             normed=True
                                             )
            residual_thres = residual_width_scale * np.abs(fit_params[0]) + \
                             fit_params[2]
            plt.plot(x_fit, y_fit / y_fit.max() * counts.max())
            plt.plot(bin_edges[:-1], counts)
            plt.title("mean : %f, var : %f" % \
                      (fit_params[2], fit_params[0]))
            plt.axvline(residual_thres, color='k', linewidth=2)
            plt.xlim(-0.3, 0.3)
            plt.savefig('/usr/users/ezbc/Desktop/residual_hist.png')
            #plt.show()


    # Include only residuals within 3 sigma
    #fit_parms[2] = 0
    intercept = fit_params[2]
    residual_thres = residual_width_scale * np.abs(fit_params[0]) + intercept
    mask = residuals > residual_thres

    print('\t\t\tResidual threshold = {0:.2f} [mag]'.format(residual_thres))
    print('\t\t\tResidual intercept = {0:.2f} [mag]'.format(intercept))

    if 'residual_hist_filename_base' in plot_args:
        import os

        x_fit = np.linspace(-10,
                            10,
                            1000)

        y_fit = gauss(x_fit, *fit_params)
        y_fit / np.nanmax(residuals)

        #y_fit = np.exp(g.eval(x_fit)[0])

        filename = plot_args['residual_hist_filename_base'] + \
                    plot_args['iter_ext'] + '.png'

        if plot_args['iter_ext'] == '0_0':
            os.system('rm -rf ' + plot_args['residual_hist_filename_base'] + \
                      '*')

        #print('\nSaving residual mask PDF figure to\n' + results_filename)
        plot_mask_residuals(residuals=residuals,
                            x_fit=x_fit,
                            y_fit=y_fit,
                            residual_thres=residual_thres,
                            filename=filename,
                            show=0)

    if 'residual_map_filename_base' in plot_args:
        import os

        filename = plot_args['residual_map_filename_base'] + \
                    plot_args['iter_ext'] + '.png'

        if plot_args['iter_ext'] == '0_0':
            os.system('rm -rf ' + plot_args['residual_map_filename_base'] + '*')

        #print('\nSaving residual mask PDF figure to\n' + results_filename)
        plot_residual_map(residuals=residuals,
                          header=plot_args['av_header'],
                          filename=filename,
                          show=0)

    return mask, residual_thres, fit_params

def _check_file(filename, clobber=False, verbose=False):

    import os

    exists = False

    if os.path.isfile(filename) or os.path.isdir(filename):
        exists = True
        if verbose:
            print('\tImage {:s} exists'.format(filename))
        if clobber:
            if verbose:
                print('\tDeleting image {:s}'.format(filename))
            os.system('rm -rf {:s}'.format(filename))
            exists = False

    return exists

def save(cloud, filename):

    import pickle

    with open(filename, 'wb') as output:
        pickle.dump(cloud, output)

def load(filename):

    import pickle

    with open(filename, 'rb') as input:
        cloud = pickle.load(input)

    return cloud

def plot_likelihoods_hist(cloud=None, props=None, limits=None,
        returnimage=False, plot_axes=('widths','dgr_grid'), show=0,
        filename='', contour_confs=(0.95,)):

    ''' Plots a heat map of likelihoodelation values as a function of velocity
    width and velocity center.

    Parameters
    ----------
    cloud : cloudpy.Cloud
        If provided, properties taken from cloud.props.


    '''

    # Import external modules
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    if cloud is not None:
        props = cloud.props

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    font_scale = 9
    params = {
              'figure.figsize': (3.6, 3.6),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    fig, ax_image = plt.subplots()

    # Define parameters from cloud
    if cloud is not None:
        props = cloud.props

    if plot_axes[0] == 'widths':
        x_grid = props['width_grid']
        x_confint = (props['hi_velocity_width']['value'],
                     props['hi_velocity_width_error']['value'][0],
                     props['hi_velocity_width_error']['value'][1],
                     )
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'Velocity Width [km/s]')
        x_sum_axes = (1, 2)
        y_pdf_label = 'Width PDF'
        if limits is None:
            x_limits = (x_grid[0], x_grid[-1])
        else:
            x_limits = limits[:2]
    if plot_axes[1] == 'dgrs':
        y_grid = props['dgr_grid']
        y_confint = (props['dust2gas_ratio']['value'],
                     props['dust2gas_ratio_error']['value'][0],
                     props['dust2gas_ratio_error']['value'][1],
                     )
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'DGR [10$^{-20}$ cm$^2$ mag]')
        y_sum_axes = (0, 2)
        x_pdf_label = 'DGR PDF'
        if limits is None:
            y_limits = (y_grid[0], y_grid[-1])
        else:
            y_limits = limits[2:]
    if plot_axes[1] == 'intercepts':
        y_grid = props['intercept_grid']
        y_confint = (props['intercept']['value'],
                     props['intercept_error']['value'][0],
                     props['intercept_error']['value'][1],
                     )
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'Intercept [mag]')
        y_sum_axes = (0, 1)
        x_pdf_label = 'Intercept PDF'
        if limits is None:
            y_limits = (y_grid[0], y_grid[-1])
        else:
            y_limits = limits[2:]


    # Create axes
    sum_axes = np.array((x_sum_axes, y_sum_axes))
    sum_axis = np.argmax(np.bincount(np.ravel(sum_axes)))

    # Mask NaNs
    likelihoods = np.array(props['likelihoods'])
    image = np.ma.array(likelihoods, mask=np.isnan(likelihoods))

    # Create likelihood image
    image = np.sum(likelihoods, axis=sum_axis) / np.sum(likelihoods)

    if 0:
        from mystats import rv3d_discrete
        import triangle
        samples = rv3d_discrete(likelihoods,
                                param_grid1=props['width_grid'],
                                param_grid2=props['dgr_grid'],
                                param_grid3=props['intercept_grid'],
                                L_scalar=100000,
                                ).pdf
        print samples.shape
        fig = triangle.corner(samples, labels=["width", "dgr", "int"],)
        fig.savefig("/usr/users/ezbc/Desktop/triangle.png")

    # Derive marginal distributions of both centers and widths
    x_sum = np.sum(likelihoods, axis=x_sum_axes)
    x_pdf = x_sum / np.sum(x_sum)
    y_sum = np.sum(likelihoods, axis=y_sum_axes)
    y_pdf = y_sum / np.sum(y_sum)

    extent = np.ravel(np.array((x_extent, y_extent)))

    #plt.rc('text', usetex=False)
    im = ax_image.imshow(image.T, interpolation='nearest', origin='lower',
            extent=extent,
            #cmap=plt.cm.gist_stern,
            #cmap=plt.cm.gray,
            cmap=plt.cm.binary,
            #norm=matplotlib.colors.LogNorm(),
            aspect='auto',
            )

    show_pdfs = 1

    if show_pdfs:
        divider = make_axes_locatable(ax_image)
        ax_pdf_x = divider.append_axes("top", 0.6, pad=0.1, sharex=ax_image)
        ax_pdf_y  = divider.append_axes("right", 0.6, pad=0.1,
                sharey=ax_image)

        # make some labels invisible
        plt.setp(ax_pdf_x.get_xticklabels() + \
                 ax_pdf_y.get_yticklabels(),
                 visible=False)

        ax_pdf_x.plot(x_grid,
                      x_pdf,
                      color='k',
                      drawstyle='steps-mid',
                      linewidth=2,
                      )

        ax_pdf_y.plot(y_pdf,
                      y_grid,
                      color='k',
                      drawstyle='steps-mid',
                      linewidth=2,
                      )

        #axHistx.axis["bottom"].major_ticklabels.set_visible(False)

        # Tick marks on the pdf?
        pdf_ticks = False

        for tl in ax_pdf_x.get_xticklabels():
            tl.set_visible(False)

        if pdf_ticks:
            wmax = x_pdf.max()
            ticks = [0, 0.5*wmax, 1.0*wmax]
            tick_labels = ['{0:.1f}'.format(ticks[0]),
                           '{0:.1f}'.format(ticks[1]),
                           '{0:.1f}'.format(ticks[2]),
                            ]
            ax_pdf_x.set_yticks(ticks)
            ax_pdf_x.set_yticklabels(tick_labels)
        else:
            for tl in ax_pdf_x.get_yticklabels():
                tl.set_visible(False)

        ax_pdf_x.set_ylabel(y_pdf_label)

        for tl in ax_pdf_y.get_yticklabels():
            tl.set_visible(False)
        if pdf_ticks:
            cmax = y_pdf.max()
            ticks = [0, 0.5*cmax, 1.0*cmax]
            tick_labels = ['{0:.1f}'.format(ticks[0]),
                           '{0:.1f}'.format(ticks[1]),
                           '{0:.1f}'.format(ticks[2]),
                            ]
            ax_pdf_y.set_xticks(ticks)
            ax_pdf_y.set_xticklabels(tick_labels)
        else:
            for tl in ax_pdf_y.get_xticklabels():
                tl.set_visible(False)

        ax_pdf_y.set_xlabel(x_pdf_label)

        # Show confidence limits
        if y_confint is not None:
            ax_pdf_y.axhspan(y_confint[0] - y_confint[1],
                             y_confint[0] + y_confint[2],
                             color='k',
                             linewidth=1,
                             alpha=0.2)
            ax_pdf_y.axhline(y_confint[0],
                             color='k',
                             linestyle='--',
                             linewidth=3,
                             alpha=1)
        if x_confint is not None:
            ax_pdf_x.axvspan(x_confint[0] - x_confint[1],
                                 x_confint[0] + x_confint[2],
                                  color='k',
                                 linewidth=1,
                                  alpha=0.2)
            ax_pdf_x.axvline(x_confint[0],
                                 color='k',
                                 linestyle='--',
                                 linewidth=3,
                                 alpha=1)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    #cb.set_label_text(r'log L')

    # Plot contours
    if contour_confs is not None:

        fractions = (1.0 - np.asarray(contour_confs))
        levels = (fractions * image.max())

        cs = ax_image.contour(image.T, levels=levels,
                extent=extent,
                colors='k'
                )

        # Define a class that forces representation of float to look a certain
        # way This remove trailing zero so '1.0' becomes '1'
        class nf(float):
             def __repr__(self):
                 str = '%.1f' % (self.__float__(),)
                 if str[-1]=='0':
                     return '%.0f' % self.__float__()
                 else:
                     return '%.1f' % self.__float__()

        # Recast levels to new class
        cs.levels = [nf(val) for val in np.asarray(contour_confs)*100.0]

        #fmt = {}
        #for level, fraction in zip(cs.levels, fractions):
        #    fmt[level] = fraction
        fmt = '%r %%'

        ax_image.clabel(cs, cs.levels, fmt=fmt, fontsize=9, inline=1)

        if limits is None:
            try:
                x_path = cs.collections[0].get_paths()[0].vertices[:,0]
                y_path = cs.collections[0].get_paths()[0].vertices[:,1]

                x_offset = np.max(x_path)*0.25
                y_offset = np.max(y_path)*0.25

                x_min, x_max = np.min(x_path) - x_offset,\
                               np.max(x_path) + x_offset
                y_min, y_max = np.min(y_path) - y_offset,\
                               np.max(y_path) + y_offset

                try:
                    ax_image.set_xlim((x_min, x_max))
                    ax_image.set_ylim((y_min, y_max))
                except UnboundLocalError:
                    pass
            except IndexError:
                pass

    if limits is not None:
        ax_image.set_xlim(limits[0],limits[1])
        ax_image.set_ylim(limits[2],limits[3])

    if 0:
    #if npix is not None or av_threshold is not None:
        text = ''
        if npix is not None:
            text += r'N$_{\rm pix}$ = ' + \
                     '{0:.0f}'.format(npix)
            if av_threshold is not None:
                text += '\n'
        if av_threshold is not None:
            text += r'$A_V$ threshold = {0:.1f} mag'.format(av_threshold)
            text += '\n'
        text += r'DGR = {0:.2f} '.format(y_confint[0]) + \
                r'$\times$ 10$^{-20}$ (cm$^2$ mag$^1$)'
        text += '\n'
        text += r'Velocity width = {0:.2f} '.format(x_confint[0]) + \
                r'km/s'
        ax_image.annotate(text,
                xytext=(0.95, 0.95),
                xy=(0.95, 0.95),
                textcoords='axes fraction',
                xycoords='axes fraction',
                color='k',
                fontsize=font_scale*0.75,
                bbox=dict(boxstyle='round',
                          facecolor='w',
                          alpha=0.3),
                horizontalalignment='right',
                verticalalignment='top',
                )

    if filename is not None:
        plt.draw()
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.draw()
        plt.show()
    if returnimage:
        return likelihoods

def plot_mask_residuals(residuals=None, x_fit=None, y_fit=None,
        fit_params=None, residual_thres=None, filename=None, show=True,
        anno_text=None, return_fig=False, limits=None, counts=None,
        bin_edges=None):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from scipy.integrate import simps as integrate
    from mystats import gauss

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    font_scale = 9
    params = {
              'figure.figsize': (3.6, 3.6),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    ax = fig.add_subplot(111)

    xlim = (-1, 2.0)


    if residuals is not None:
        residuals_nonans = np.ravel(residuals[~np.isnan(residuals)])
        counts, bin_edges = \
            np.histogram(residuals_nonans,
                         bins=20,
                         )
    else:
        assert counts is not None

    bin_edges_ext = np.zeros(len(counts) + 1)
    counts_ext = np.zeros(len(counts) + 1)

    bin_edges_ext[0] = bin_edges[0] - (bin_edges[1] - bin_edges[0])
    bin_edges_ext[1:] = bin_edges[:-1]
    counts_ext[0] = 0
    counts_ext[1:] = counts
    bin_edges_ext = np.hstack([xlim[0], bin_edges_ext, xlim[1]])
    counts_ext = np.hstack([0, counts_ext, 0])

    # Normalize so area = 1
    #counts_ext /= np.nansum(counts_ext) * (bin_edges_ext[2] - bin_edges_ext[1])
    #counts_ext = counts_ext / integrate(counts_ext, x=bin_edges_ext)
    #counts_ext /= counts_ext.max()
    if fit_params is not None:
        x_fit = np.linspace(-3,
                            3,
                            1000)

        y_fit = gauss(x_fit, *fit_params)
        #y_fit / np.nanmax(residuals)
    y_fit /= np.max(y_fit)
    y_fit *= np.max(counts_ext)

    pdf_max = bin_edges_ext[counts_ext == counts_ext.max()][0]
    #print('\t\t\tResidual max = {0:.2f} [mag]'.format(pdf_max))

    ax.plot(bin_edges_ext, counts_ext, drawstyle='steps-pre',
            linewidth=1.5)
    ax.plot(x_fit, y_fit,
            linewidth=2,
            alpha=0.6)


    ax.set_xlim([-0.5,0.5])
    #ax.set_xlim(xlim)
    ax.set_ylim(limits)
    ax.axvline(residual_thres,
               color='k',
               linestyle='--',
               linewidth=1.5)
    ax.set_xlabel(r'Residual $A_V$ [mag]')
    ax.set_ylabel('Counts')

    if anno_text is not None:
        ax.set_title(anno_text)
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=60)
    if show:
        plt.show()
    if return_fig:
        return fig

def plot_residual_map(residuals, header=None, dgr=None, show=False,
        filename=None, anno_text=None, return_fig=False, vlimits=[None,None]):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon

    # Set up plot aesthetics
    plt.clf()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 15
    params = {
              'figure.figsize': (3.6, 3.6),
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    nrows_ncols=(1,1)
    ngrids=1

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # create axes
    ax = imagegrid[0]
    cmap = cm.jet # colormap
    # show the image
    im = ax.imshow(residuals,
            interpolation='nearest',
            origin='lower',
            cmap=cmap,
            #norm=matplotlib.colors.LogNorm()
            vmin=vlimits[0],
            vmax=vlimits[1],
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)

    if anno_text is not None:
        ax.set_title(anno_text)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    limits = None
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Write label to colorbar
    cb.set_label_text(r'$A_V$ [Mag]',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=60)
    if show:
        plt.show()
    if return_fig:
        return fig

def plot_av_bin_map(av_map, av_bin_map, av_header=None, av_header_bin=None,
        filename=None):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon

    # Set up plot aesthetics
    plt.clf()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 15
    params = {
              'figure.figsize': (3.6, 3.6),
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    nrows_ncols=(1,1)
    ngrids=1

    # Original map
    # ------------
    imagegrid = ImageGrid(fig, (2,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=1,
                 axes_class=(wcs.Axes,
                             dict(header=av_header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # create axes
    ax = imagegrid[0]
    cmap = cm.jet # colormap
    # show the image
    im = ax.imshow(av_map,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            #norm=matplotlib.colors.LogNorm()
            #vmin=-1,
            #vmax=1
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    limits = None
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Write label to colorbar
    cb.set_label_text(r'$A_V$ [Mag]',)

    # Binned map
    # ----------
    imagegrid = ImageGrid(fig, (2,1,2),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=1,
                 axes_class=(wcs.Axes,
                             dict(header=av_header_bin)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # create axes
    ax = imagegrid[0]
    cmap = cm.jet # colormap
    # show the image
    im = ax.imshow(av_bin_map,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            #norm=matplotlib.colors.LogNorm()
            #vmin=-1,
            #vmax=1
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    limits = None
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Write label to colorbar
    cb.set_label_text(r'$A_V$ [Mag]',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')

def plot_dgr_intercept_progression(cloud, filename=None, show=True,
        title=None):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from scipy.integrate import simps as integrate

    nrows = len(cloud.iter_vars)

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    font_scale = 9
    params = {
              'figure.figsize': (2 * nrows, 2 * nrows),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig, axes = plt.subplots(nrows=2, ncols=nrows)

    for i in xrange(0, nrows):

        ax1 = axes[0, i]
        ax2 = axes[1, i]

        iter_dict = cloud.iter_vars[i]['masking_var']

        n_points = len(iter_dict)
        dgrs = np.empty(n_points)
        intercepts = np.empty(n_points)
        npix = np.empty(n_points)

        for j, iteration in enumerate(iter_dict):
            dgrs[j] = iter_dict[iteration]['dgr']
            intercepts[j] = iter_dict[iteration]['intercept']
            npix[j] = np.sum(~iter_dict[iteration]['mask'])

        ax1.plot(npix, dgrs,
                 linestyle='',
                 marker='o',
                 )
        ax2.plot(npix, intercepts,
                 linestyle='',
                 marker='o',
                 )

        ax2.set_xlabel(r'Number of Unmasked Pixels')
        if i == 0:
            ax1.set_ylabel(r'DGR [10$^{-20}$ cm$^2$ mag]')
            ax2.set_ylabel(r'Intercept [mag]')

        ax1.locator_params(nbins = 4)
        ax2.locator_params(nbins = 4)

        vel_range = cloud.iter_vars[i]['vel_range']
        width = np.abs(vel_range[1] - vel_range[0])

        ax1.set_title(r'$\Delta_V$ = {0:.0f} km/s'.format(width))


    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        fig.tight_layout()
        plt.savefig(filename, bbox_inches='tight', dpi=600)
    if show:
        plt.show()

def plot_map_movie(cloud, filename=None,):

    import moviepy.editor as mpy
    from moviepy.video.io.bindings import mplfig_to_npimage

    residuals_list = []
    anno_text_list = []

    resid_min_abs = np.Inf
    resid_max_abs = -np.Inf

    for parent_iter in cloud.iter_vars:
        mask_dict = cloud.iter_vars[parent_iter]['masking_var']
        for mask_iter in mask_dict:
            anno_text_list.append('Parent iter = ' + \
                                  '{0:02d}'.format(parent_iter) + \
                                  ', mask iter = {0:02d}'.format(mask_iter))
            # grab residuals, record global limits
            residuals = mask_dict[mask_iter]['residuals']
            residuals_list.append(residuals)
            resid_min, resid_max = np.nanmin(residuals), np.nanmax(residuals)
            if resid_min < resid_min_abs:
                resid_min_abs = resid_min
            if resid_max > resid_max_abs:
                resid_max_abs = resid_max

    def make_frame(t):
        fig = plot_residual_map(residuals_list[int(t)],
                                anno_text=anno_text_list[int(t)],
                                return_fig=True,
                                header=cloud.av_header,
                                vlimits=[resid_min_abs, resid_max_abs])

        return mplfig_to_npimage(fig)

    # Make the movie
    animation = mpy.VideoClip(make_frame, duration=len(residuals_list))
    if filename.split('.',1)[1] == 'mp4':
        animation.write_videofile(filename, fps=2)
    elif filename.split('.',1)[1] == 'gif':
        animation.write_gif(filename, fps=2)

def plot_residual_hist_movie(cloud, filename=None,):

    import moviepy.editor as mpy
    from moviepy.video.io.bindings import mplfig_to_npimage

    residuals_list = []
    residual_thres_list = []
    fit_params_list = []
    anno_text_list = []

    resid_max_abs = -np.Inf

    for parent_iter in cloud.iter_vars:
        mask_dict = cloud.iter_vars[parent_iter]['masking_var']
        for mask_iter in mask_dict:
            anno_text_list.append('Parent iter = ' + \
                                  '{0:02d}'.format(parent_iter) + \
                                  ', mask iter = {0:02d}'.format(mask_iter))

            fit_params_list.append(mask_dict[mask_iter]['fit_params'])
            residual_thres_list.\
                    append(mask_dict[mask_iter]['residual_threshold'])

            # grab residuals, record global limits
            residuals = mask_dict[mask_iter]['residuals']
            counts, bin_edges = \
                np.histogram(residuals[residuals == residuals],
                             bins=20,
                             )
            residuals_list.append([counts, bin_edges])
            resid_max = np.nanmax(counts)
            if resid_max > resid_max_abs:
                resid_max_abs = resid_max

    def make_frame(t):
        fig = plot_mask_residuals(counts=residuals_list[int(t)][0],
                                  bin_edges=residuals_list[int(t)][1],
                                  anno_text=anno_text_list[int(t)],
                                  fit_params=fit_params_list[int(t)],
                                  residual_thres=residual_thres_list[int(t)],
                                  return_fig=True,
                                  limits=[-0.1,
                                          resid_max_abs * 1.1])

        return mplfig_to_npimage(fig)

    # Make the movie
    animation = mpy.VideoClip(make_frame, duration=len(residuals_list))
    if filename.split('.',1)[1] == 'mp4':
        animation.write_videofile(filename, fps=2)
    elif filename.split('.',1)[1] == 'gif':
        animation.write_gif(filename, fps=2)

