#!/usr/bin/python

import cloudpy
import numpy as np

def Test_get_residual_mask():

    from numpy.testing import assert_almost_equal

    residuals = np.random.normal(loc=-1, scale=1.0, size=1000)

    mask, intercept = cloudpy._get_residual_mask(residuals)

    assert_almost_equal(-1, intercept, decimal=0)

if 1:
    class TestCloudpy():

        def setup(self):

            # define directory locations
            # --------------------------
            self.output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
            self.figure_dir = \
                '/d/bip3/ezbc/california/figures/'
            self.av_dir = '/d/bip3/ezbc/california/data/av/'
            self.hi_dir = '/d/bip3/ezbc/california/data/hi/'
            self.co_dir = '/d/bip3/ezbc/california/data/co/'
            self.core_dir = \
                    '/d/bip3/ezbc/california/data/python_output/core_properties/'
            self.property_dir = '/d/bip3/ezbc/california/data/python_output/'
            self.region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
            self.likelihood_dir = \
                    '/d/bip3/ezbc/california/data/python_output/nhi_av/'

            # define filenames
            self.prop_filename = self.property_dir + \
                    'california_global_properties.txt'
            self.av_filename = self.av_dir + 'california_av_planck_5arcmin.fits'
            self.av_error_filename = self.av_dir + \
                    'california_av_planck_5arcmin.fits'
            self.hi_filename = self.hi_dir + \
                    'california_hi_galfa_cube_regrid_planckres.fits'
            self.hi_error_filename = self.hi_dir + \
                    'california_hi_galfa_cube_regrid_planckres_noise.fits'

            self.region = 'california'
            self.region_filename = self.region_dir + 'multicloud_divisions.reg'

            vel_widths = np.arange(1, 75, 6*0.16667)
            dgrs = np.arange(0.001, 0.3, 5e-3)
            intercepts = np.arange(-1, 1, 0.1)

            self.hi_widths = vel_widths
            self.dgrs = dgrs
            self.intercepts = intercepts

            self.cloud = cloudpy.cloud(self.av_filename,
                                  self.hi_filename,
                                  av_error_filename=self.av_error_filename,
                                  hi_error_filename=self.hi_error_filename,
                                  cloud_prop_filename=self.prop_filename,
                                  dgrs=self.dgrs,
                                  intercepts=self.intercepts,
                                  hi_widths=self.hi_widths
                                  )

        def test_iterate_mle_calc(self):

            self.cloud.iter_step = 0
            self.cloud.region = self.region

            # Change WCS coords to pixel coords of images
            self.cloud._convert_coordinates()

            # Load cloud division regions from ds9
            self.cloud.load_region(self.region_filename)

            self.cloud._iterate_mle_calc(hi_vel_range=[-10,10],
                                      hi_noise_vel_range=[90,110],
                                      )

#Test_get_residual_mask()



