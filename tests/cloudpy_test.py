#!/usr/bin/python

import cloudpy
import numpy as np

if 1:
    def Test_get_residual_mask():

        from numpy.testing import assert_almost_equal

        residuals = np.random.normal(loc=-1, scale=1.0, size=10000)
        residuals[:3000] = np.random.uniform(high=3, size=3000)

        mask, intercept = cloudpy._get_residual_mask(residuals, plot_args={})

        print('true = -1', 'found = ', intercept)

        assert_almost_equal(-1, intercept, decimal=0)

if 0:
    class TestCloudpy():

        def setup(self):

            # define directory locations
            # --------------------------
            self.output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
            self.figure_dir = \
                '/d/bip3/ezbc/perseus/figures/'
            self.av_dir = '/d/bip3/ezbc/perseus/data/av/'
            self.hi_dir = '/d/bip3/ezbc/perseus/data/hi/'
            self.co_dir = '/d/bip3/ezbc/perseus/data/co/'
            self.core_dir = \
                    '/d/bip3/ezbc/perseus/data/python_output/core_properties/'
            self.property_dir = '/d/bip3/ezbc/perseus/data/python_output/'
            self.region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
            self.likelihood_dir = \
                    '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'

            # define filenames
            self.prop_filename = self.property_dir + \
                    'perseus_global_properties.txt'
            self.av_filename = self.av_dir + \
                    'perseus_av_planck_tau353_5arcmin.fits'
            self.av_error_filename = self.av_dir + \
                    'perseus_av_error_planck_tau353_5arcmin.fits'
            self.hi_filename = self.hi_dir + \
                    'perseus_hi_galfa_cube_regrid_planckres.fits'
            self.hi_error_filename = self.hi_dir + \
                    'perseus_hi_galfa_cube_regrid_planckres_noise.fits'

            self.region = 'perseus'
            self.region_filename = self.region_dir + 'multicloud_divisions.reg'
            self.cloud_filename = \
                    '/d/bip3/ezbc/multicloud/data/cloud.pickle'

            width_grid = np.arange(1, 75, 6*0.16667)
            dgr_grid = np.arange(0.001, 0.3, 5e-3)
            intercept_grid = np.arange(-2, 2, 0.1)

            self.width_grid = width_grid
            self.dgr_grid = dgr_grid
            self.intercept_grid = intercept_grid

            # Define number of pixels in each bin
            self.binsize = 1.0 * 60.0 / 5.0

            self.cloud = cloudpy.Cloud(self.av_filename,
                                  self.hi_filename,
                                  av_error_filename=self.av_error_filename,
                                  hi_error_filename=self.hi_error_filename,
                                  cloud_prop_filename=self.prop_filename,
                                  dgr_grid=self.dgr_grid,
                                  intercept_grid=self.intercept_grid,
                                  width_grid=self.width_grid,
                                  residual_width_scale=3.0,
                                  threshold_delta_dgr=0.0001,
                                  hi_noise_vel_range=[90,110],
                                  vel_range_diff_thres=10,
                                  init_vel_range=[-50,50],
                                  verbose=True,
                                  clobber_likelihoods=True,
                                  binsize=self.binsize,
                                  )

            if 0:
                cloudpy.save(self.cloud,
                             self.cloud_filename.replace('cloud.pickle',
                                                         'cloud_init.pickle'))

        if 0:
            def test_run_anaylsis(self):

                self.cloud.run_analysis(region_filename=self.region_filename,
                                        region=self.region)

                props, iter_vars = self.cloud.props, self.cloud.iter_vars

                print('\nSaving props...')
                cloudpy.save(props,
                             '/d/bip3/ezbc/multicloud/data/props.pickle')
                cloudpy.save(iter_vars,
                                   '/d/bip3/ezbc/multicloud/data/iter_vars.pickle')

                print('\nSaving cloud...')
                cloudpy.save(self.cloud, self.cloud_filename)


        if 0:
            def test_write_final_params(self,):
                props = \
                    cloudpy.load(\
                        '/d/bip3/ezbc/multicloud/data/props.pickle')
                iter_vars = \
                    cloudpy.load(\
                        '/d/bip3/ezbc/multicloud/data/iter_vars.pickle')

                self.cloud.props = props
                self.cloud.iter_vars = iter_vars

                self.cloud._write_final_params()

        if 1:
            def test_plotting(self):
                #cloud = cloudpy.load(self.cloud_filename)

                #cloud._write_final_params()

                props = \
                    cloudpy.load(\
                        '/d/bip3/ezbc/multicloud/data/props.pickle')

                likelihood_filename = \
                        '/d/bip3/ezbc/perseus/figures/likelihood/' + \
                        'perseus_likelihood_planck_bin_scaled_wd.png'

                cloudpy.plot_likelihoods_hist(props=props,
                                  plot_axes=('widths', 'dgrs'),
                                  show=0,
                                  returnimage=False,
                                  filename=likelihood_filename,
                                  limits=[0, 15, 0.1, 0.2],
                                  )

        if 0:
            def test_iterate_mle_calc(self):

                self.cloud.iter_step = 0
                self.cloud.region = self.region

                # Change WCS coords to pixel coords of images
                self.cloud._convert_coordinates()

                # Load cloud division regions from ds9
                self.cloud.load_region(self.region_filename)

                self.cloud._iterate_mle_calc(hi_vel_range=[-10,10],
                                          )




#Test_get_residual_mask()



