#!/usr/bin/python

# Import external modules
import numpy as np

class cloud():

    def __init__(self, av_filename, hi_filename, av_error_filename=None,
            hi_error_filename=None, cloud_prop_filename=None, dgrs=None,
            intercepts=None, hi_widths=None, threshold_delta_dgr=0.0001,
            resid_width_scale=3.0,
            verbose=False,
            plot_vars={'plot_progress':False,'residual_filename':None}):

        # Import external modules
        import pyfits as fits
        from mycoords import make_velocity_axis
        from os import system,path

        # Define local variables
        self.av_filename = av_filename
        self.hi_filename = hi_filename
        self.av_error_filename = av_error_filename
        self.hi_error_filename = hi_error_filename
        self.dgrs = dgrs
        self.hi_widths = hi_widths
        self.intercepts = intercepts
        self.threshold_delta_dgr = threshold_delta_dgr
        self.resid_width_scale = resid_width_scale
        self.plot_vars = plot_vars
        self.verbose = verbose
        self.iter_var = {}

        # Load data
        self.av_data, self.av_header = fits.getdata(av_filename, header=True)
        self.hi_data, self.hi_header = fits.getdata(hi_filename, header=True)
        if av_error_filename is not None:
            self.av_error_data, self.av_error_header = \
                    fits.getdata(av_error_filename, header=True)
        else:
            self.av_error_data = None
        if hi_error_filename is not None and path.isfile(hi_error_filename):
            self.hi_error_data, self.hi_error_header = \
                fits.getdata(hi_error_filename, header=True)
        else:
            self.hi_error_data = None
        if cloud_prop_filename is not None:
            self.load_cloud_properties(cloud_prop_filename)

        # Make velocity axis for HI cube
        self.hi_vel_axis = make_velocity_axis(self.hi_header)

    def _derive_region_mask(self,):

        import mygeometry as myg

        # Derive relevant region
        region_vertices = \
            self.props['regions'][self.region]['poly_verts']['pixel']

        # block off region
        region_mask = np.logical_not(myg.get_polygon_mask(self.av_data,
                                                          region_vertices))

        self.region_mask = region_mask

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

    def _iterate_mle_calc(self, hi_vel_range=None, hi_noise_vel_range=None,):

        from myimage_analysis import calculate_nhi, calculate_noise_cube,\
                                     bin_image
        from os import system,path

        self.iter_var[self.iter_step] = {}
        self.iter_var[self.iter_step]['hi_vel_range'] = hi_vel_range

        # derive the region mask
        self._derive_region_mask()

        # check if HI error cube present
        if self.hi_error_data is None:
            self.hi_error_data = \
                    calculate_noise_cube(cube=self.hi_data,
                                    velocity_axis=self.hi_vel_axis,
                                    velocity_noise_range=hi_noise_vel_range,
                                    header=self.hi_header,
                                    Tsys=30.0,
                                    filename=self.hi_error_filename)

        # Derive mask by excluding correlated residuals
        # ---------------------------------------------
        # Derive initial N(HI) image
        self.nhi_image, self.nhi_error_image = \
                calculate_nhi(cube=self.hi_data,
                              velocity_axis=self.hi_vel_axis,
                              velocity_range=hi_vel_range,
                              noise_cube=self.hi_error_data,
                              velocity_noise_range=hi_noise_vel_range,
                              Tsys=30.0,
                              return_nhi_error=True,
                              )

        # Write iteration step variables
        self.iter_var[self.iter_step]['nhi_image'] = self.nhi_image
        self.iter_var[self.iter_step]['nhi_image_error'] = self.nhi_error_image

        av_model, mask, dgr, intercepts, masking_results = \
            _iterate_residual_masking(nhi_image=self.nhi_image,
                                     nhi_image_error=self.nhi_error_image,
                                     av_data=self.av_data,
                                     av_data_error=self.av_error_data,
                                     vel_range=hi_vel_range,
                                     dgrs=self.dgrs,
                                     intercepts=self.intercepts,
                                     threshold_delta_dgr=\
                                             self.threshold_delta_dgr,
                                     resid_width_scale=self.resid_width_scale,
                                     init_mask=self.region_mask,
                                     verbose=self.verbose,
                                     plot_progress=\
                                             self.plot_vars['plot_progress'],
                                     results_filename=\
                                             self.plot_vars['residual_filename']
                                     )

        # last left on line  of
        # california_analysis_likelihood_iterative_mask_binning.py

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

        done = false
        self.iter_step = 0
        self.region = region

        # Change WCS coords to pixel coords of images
        self._convert_coordinates()

        # Load cloud division regions from ds9
        self.load_region(region_filename)

        while not done:
            self._iterate_mle_calc(region_filename=region_filename,
                                   region=region,
                                   hi_vel_range=hi_vel_range,
                                   hi_noise_vel_range=hi_noise_vel_range,
                                   )

            done = True

def _iterate_residual_masking(
                             nhi_image=None,
                             nhi_image_error=None,
                             av_data=None,
                             av_data_error=None,
                             init_mask=None,
                             vel_range=None,
                             dgrs=None,
                             intercepts=None,
                             threshold_delta_dgr=None,
                             resid_width_scale=3.0,
                             plot_progress=False,
                             results_filename=None,
                             verbose=False,
                             return_masking_results=True,
                             ):

    '''

    Returns
    -------
    av_model : numpy array
    mask : numpy array
    dgr : float

    '''

    import numpy as np

    # Mask out nans
    mask = (np.isnan(av_data) | \
            np.isnan(av_data_error) | \
            (av_data_error == 0) | \
            np.isnan(nhi_image) | \
            np.isnan(nhi_image_error) | \
            (nhi_image_error == 0))

    # Apply initial mask to exclude throughout process
    if init_mask is not None:
        mask += init_mask

    # solve for DGR using linear least squares
    if verbose:
        print('\nBeginning iterative DGR calculations + masking...')

    # Iterate masking pixels which are correlated and rederiving a linear
    # least squares solution for the DGR
    # ----------------------------------------------------------------------
    masking_results = {}
    use_intercept = True
    delta_dgr = 1e10
    dgr = 1e10
    iteration = 0
    dgr_list = []
    width_list = []
    intercept_list = []
    mask_list = []

    while delta_dgr > threshold_delta_dgr:
        results = _calc_likelihoods(
                         nhi_image=nhi_image[~mask],
                         av_image=av_data[~mask],
                         av_image_error=av_data_error[~mask],
                         #image_weights=bin_weights[~mask],
                         #vel_center=vel_center_masked,
                         vel_widths=np.arange(0,1,1),
                         dgrs=dgrs,
                         intercepts=intercepts,
                         results_filename='',
                         return_likelihoods=True,
                         likelihood_filename=None,
                         clobber=False,
                         verbose=False
                         )

        # Unpack output of likelihood calculation
        (vel_range_confint, width_confint, dgr_confint, intercepts_confint,
                likelihoods, width_likelihood, dgr_likelihood,
                intercept_likelihood, width_max, dgr_max, intercept_max,
                vel_range_max) = results

        dgr_new = dgr_max
        intercept = intercept_max

        # Create model with the DGR
        if verbose:
            print('Iteration {0:.0f} results:'.format(iteration))
            print('\tDGR = {0:.2} 10^20 cm^2 mag'.format(dgr_new))
            print('\tIntercept = {0:.2f} mag'.format(intercept))
            print('')

        av_image_model = nhi_image * dgr_new + intercept

        residuals = av_data - av_image_model

        residuals[mask] = np.nan

        # plot progress
        if 0:
            plot_residual_map(residuals,
                              header=av_header,
                              dgr=dgr_new)

        # Include only residuals which are white noise
        if iteration == 0:
            plot_filename = results_filename
        else:
            plot_filename = None

        mask_new, intercept = \
                get_residual_mask(residuals,
                                     resid_width_scale=resid_width_scale,
                                     plot_progress=plot_progress,
                                     results_filename=plot_filename)

        #intercepts = np.linspace(intercept, intercept + 1.0, 1.0)

        # Mask non-white noise, i.e. correlated residuals.
        mask[mask_new] = 1

        if verbose:
            npix = mask.size - np.sum(mask)
            print('\tNumber of non-masked pixels = {0:.0f}'.format(npix))

        dgr_list.append(dgr_max)
        width_list.append(width_max)
        intercept_list.append(intercept)
        mask_list.append(mask.tolist())

        # Reset while loop conditions
        delta_dgr = np.abs(dgr - dgr_new)
        dgr = dgr_new
        iteration += 1

    #np.savetxt('/usr/users/ezbc/Desktop/dgr_data.txt', dgr_list)
    masking_results['dgrs'] = dgr_list
    masking_results['widths'] = width_list
    masking_results['intercepts'] = intercept_list
    masking_results['masks'] = mask_list

    # Plot results
    if 0:
        mask_new = get_residual_mask(residuals,
                                     resid_width_scale=resid_width_scale,
                                     plot_progress=plot_progress,
                                     results_filename=results_filename)

    # Create model of Av
    av_model = dgr * nhi_image
    av_model[mask] = np.nan

    output = [av_model, mask, dgr, intercepts,]
    if return_masking_results:
        output.append(masking_results)

    return output

def _calc_likelihoods(
        hi_cube=None,
        hi_vel_axis=None,
        nhi_image=None,
        av_image=None,
        av_image_error=None,
        image_weights=None,
        vel_center=None,
        vel_widths=None,
        dgrs=None,
        intercepts=None,
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
    hi_vel_range : tuple
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
        likelihoods = np.zeros((len(vel_widths),
                                len(dgrs),
                                len(intercepts)))

        # Progress bar parameters
        total = float(likelihoods.size)
        count = 0

        for j, vel_width in enumerate(vel_widths):
            # Construct N(HI) image outside of DGR loop, then apply
            # DGRs in loop

            # use the hi cube and vel range if no nhi image provided
            if nhi_image is None:
                vel_range = np.array((vel_center - vel_width / 2.,
                                      vel_center + vel_width / 2.))
                nhi_image = calculate_nhi(cube=hi_cube,
                                          velocity_axis=hi_vel_axis,
                                          velocity_range=vel_range,
                                          return_nhi_error=False)

            # Cycle through DGR to estimate error
            for k, dgr in enumerate(dgrs):
                for m, intercept in enumerate(intercepts):
                    # Create model of Av with N(HI) and DGR
                    av_image_model = nhi_image * dgr + intercept

                    logL = calc_logL(av_image_model,
                                     av_image,
                                     data_error=av_image_error,
                                     weights=image_weights)

                    #print logL

                    likelihoods[j, k, m] = logL

                    # Shows progress each 10%
                    if verbose:
                        count += 1
                        abs_step = int((total * 1)/10) or 10
                        if count and not count % abs_step:
                            print "\t{0:.0%} processed".format(count/total)

            nhi_image = None

    # Load file of likelihoods
    elif not perform_mle:
        print('Reading likelihood grid file:')
        print(likelihood_filename)

        hdu = fits.open(likelihood_filename)
        likelihoods = hdu[0].data

        if len(vel_widths) != likelihoods.shape[0] or \
           len(dgrs) != likelihoods.shape[1]:
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

    if 0:
        import matplotlib.pyplot as plt
        plt.imshow(likelihoods[:, :, 0], origin='lower')
        plt.show()

    # Derive marginal distributions of both centers and widths
    intercept_likelihood = np.sum(likelihoods, axis=(0, 1)) / \
                                  np.sum(likelihoods)
    width_likelihood = np.sum(likelihoods, axis=(1, 2)) / \
            np.sum(likelihoods)
    dgr_likelihood = np.sum(likelihoods, axis=(0, 2)) / \
            np.sum(likelihoods)

    # Derive confidence intervals of parameters
    width_confint = calc_symmetric_error(vel_widths,
                                   width_likelihood,
                                   alpha=1.0 - conf)
    dgr_confint = calc_symmetric_error(dgrs,
                                 dgr_likelihood,
                                 alpha=1.0 - conf)
    intercept_confint = calc_symmetric_error(intercepts,
                                 intercept_likelihood,
                                 alpha=1.0 - conf)

    # Get values of best-fit model parameters
    max_loc = np.where(likelihoods == np.max(likelihoods))
    width_max = vel_widths[max_loc[0][0]]
    dgr_max = dgrs[max_loc[1][0]]
    intercept_max = intercepts[max_loc[2][0]]

    if verbose:
        print('\nVelocity widths = ' + \
                '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                                    width_confint[2],
                                                    np.abs(width_confint[1])))
        print('\nDGRs = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} 10^20 cm^2 mag'.format(dgr_confint[0],
                                                    dgr_confint[2],
                                                    np.abs(dgr_confint[1])))
        print('\nIntercepts = ' + \
        '{0:.2f} +{1:.2f}/-{2:.2f} 10^20 cm^2 mag'.format(intercept_confint[0],
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

    if not return_likelihoods:
        return vel_range_confint, dgr_confint
    else:
        return (vel_range_confint, width_confint, dgr_confint,
                intercept_confint, likelihoods,
                width_likelihood, dgr_likelihood,
                intercept_likelihood, width_max, dgr_max, intercept_max,
                vel_range_max)

def _get_residual_mask(residuals, resid_width_scale=3.0, plot_progress=False,
        results_filename=None):

    '''

    '''

    import numpy as np
    from mystats import gauss
    from scipy.optimize import curve_fit, minimize


    # Fit the rising portion of the residuals
    residuals_crop = residuals[(residuals < 0) & \
                               #(residuals > -1.5) & \
                                ~np.isnan(residuals)]

    #print('histing')
    counts, bin_edges = np.histogram(np.ravel(residuals_crop),
                                     bins=100,
                                     )

    p0=(1, np.nanmax(counts), 0)

    if 0:
        fit_params = curve_fit(gauss,
                               bin_edges[:-1],
                               counts,
                               p0=p0,
                               maxfev=1000000,
                               )[0]
    elif 1:
        from lmfit import minimize, Parameters

        # Set parameter limits and initial guesses
        params = Parameters()
        params.add('width',
                   value=p0[0],
                   min=0.1,
                   max=10,
                   )
        params.add('amp',
                   value=p0[1],
                   min=0,
                   max=2 * np.nanmax(counts),
                   )
        params.add('x0',
                   value=p0[2],
                   min=-4,
                   max=4,
                   )

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
                          method='lbfgsb')

        #print('done fitting')
        fit_params = (params['width'].value, params['amp'].value,
                params['x0'].value)
    else:
        bounds = ((0, 10), (0, 5 * np.nanmax(counts)), (-10, 10))
        fit_params = minimize(gauss,
                              counts,
                              method='L-BFGS-B',
                              bounds=bounds,)

    # Include only residuals within 3 sigma
    intercept = fit_params[2]
    residual_thres = resid_width_scale * np.abs(fit_params[0]) + intercept
    mask = residuals > residual_thres

    if 0:
        import matplotlib.pyplot as plt
        plt.clf(); plt.close();
        x_fit = np.linspace(np.nanmin(residuals),
                            np.nanmax(residuals),
                            1000)

        y_fit = gauss(x_fit, *fit_params)
        plt.plot(bin_edges[:-1], counts)
        plt.plot(x_fit, y_fit)
        plt.show()
        plt.savefig('/usr/users/ezbc/Desktop/residuals.png')

    if results_filename is not None:
        x_fit = np.linspace(np.nanmin(residuals),
                            np.nanmax(residuals),
                            1000)

        y_fit = gauss(x_fit, *fit_params)
        y_fit / np.nanmax(residuals)

        print('\nSaving residual mask PDF figure to\n' + results_filename)
        plot_mask_residuals(residuals=residuals,
                            x_fit=x_fit,
                            y_fit=y_fit,
                            residual_thres=residual_thres,
                            filename=results_filename,
                            show=plot_progress)
        plot_mask_residuals(residuals=residuals,
                            x_fit=x_fit,
                            y_fit=y_fit,
                            residual_thres=residual_thres,
                            filename=results_filename.replace('.pdf', '.png'),
                            show=plot_progress)

    return mask, intercept

