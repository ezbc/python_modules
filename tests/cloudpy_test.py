#!/usr/bin/python

import cloudpy
import numpy as np

if 0:
    def Test_get_residual_mask():

        from numpy.testing import assert_almost_equal

        residuals = np.random.normal(loc=-1, scale=1.0, size=10000)
        residuals[:3000] = np.random.uniform(high=3, size=3000)

        mask, intercept = cloudpy._get_residual_mask(residuals, plot_args={})

        print('true = -1', 'found = ', intercept)

        assert_almost_equal(-1, intercept, decimal=0)

if 1:
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

            # Plot args
            residual_hist_filename_base = self.figure_dir + \
                                          'diagnostics/residuals/' + \
                                          self.region + '_residual_hist'
            residual_map_filename_base = self.figure_dir + 'diagnostics/residuals/' + \
                                          self.region + '_residual_map'
            likelihood_filename_base = self.figure_dir + 'diagnostics/likelihoods/' + \
                                          self.region + '_likelihood'
            av_bin_map_filename_base = self.figure_dir + 'diagnostics/maps/' + \
                                          self.region + '_bin_map'

            plot_args = {
                    'residual_hist_filename_base': residual_hist_filename_base,
                    'residual_map_filename_base': residual_map_filename_base,
                    'likelihood_filename_base': likelihood_filename_base,
                    'av_bin_map_filename_base' : av_bin_map_filename_base,
                    }

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
                                  plot_args=plot_args,
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
            def test_faint_masking(self,):

                from numpy.testing import assert_array_almost_equal
                from numpy.testing import assert_almost_equal

                self.cloud.av_data = np.arange(0,9).reshape(3, 3,)

                mask_faint = self.cloud._get_faint_mask()
                mask_faint_answer = np.array([[0, 1, 1],
                                              [1, 1, 1],
                                              [1, 1, 1]])

                assert_array_almost_equal(mask_faint, mask_faint_answer)

                # Test including more pixels
                mask_new = self.cloud._add_pixels_to_mask()

                mask_new_answer = np.array([[0, 0, 1],
                                            [1, 1, 1],
                                            [1, 1, 1]])

                assert_array_almost_equal(mask_new, mask_new_answer)

                # Test combo of region mask with new mask
                mask_region = np.array([[1, 0, 0],
                                        [0, 0, 0],
                                        [0, 0, 0]])

                mask_combo = self.cloud._combine_masks(mask_region, mask_new)

                mask_combo_answer = np.array([[1, 0, 1],
                                              [1, 1, 1],
                                              [1, 1, 1]])

                assert_array_almost_equal(mask_combo, mask_combo_answer)

                # Test combo subtraction of region mask with new mask
                mask_region = np.array([[1, 0, 0],
                                        [0, 0, 0],
                                        [0, 0, 0]])

                mask_combo = self.cloud._combine_masks(mask_region, mask_new,
                                                    child_action='subtract')

                mask_combo_answer = np.array([[1, 0, 1],
                                              [1, 1, 1],
                                              [1, 1, 1]])

                assert_array_almost_equal(mask_combo, mask_combo_answer)


                mask_new = np.array([[1, 0, 1],
                                     [1, 0, 1],
                                     [1, 1, 1]])

                mask_combo = self.cloud._combine_masks(mask_region, mask_new,
                                                    child_action='subtract')

                mask_combo_answer = np.array([[1, 0, 1],
                                              [1, 1, 1],
                                              [1, 1, 1]])

                assert_array_almost_equal(mask_combo, mask_combo_answer)




        if 0:
            def test_mle_derivation(self,):

                from numpy.testing import assert_array_almost_equal
                from numpy.testing import assert_almost_equal
                from myimage_analysis import calculate_nhi

                dgr = 0.1 # cm^2 10^-20 mag
                intercept = 1 # mag
                width = 20 # km/s

                vel_range = (self.cloud.vel_center - width / 2.0,
                             self.cloud.vel_center + width / 2.0)

                nhi_image = calculate_nhi(cube=self.cloud.hi_data,
                                          velocity_axis=self.cloud.hi_vel_axis,
                                          velocity_range=vel_range)

                # Create mock Av_data
                if 0:
                    av_data_mock = dgr * nhi_image + intercept

                    self.cloud.av_data = av_data_mock

                    self.cloud.run_analysis(region_filename=self.region_filename,
                                            region=self.region)

                    print('\nSaving cloud...')
                    cloudpy.save(self.cloud, self.cloud_filename)
                else:
                    self.cloud = cloudpy.load(self.cloud_filename)

                dgr_mle = self.cloud.props['dust2gas_ratio_max']['value']
                intercept_mle = self.cloud.props['intercept_max']['value']
                width_mle = self.cloud.props['hi_velocity_width_max']['value']

                assert_almost_equal(dgr_mle, dgr, decimal=1)
                assert_almost_equal(intercept_mle, intercept, decimal=1)
                assert_almost_equal(width_mle, width, decimal=-1)

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

        if 0:
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
                                  returnimage=false,
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


if 0:
    def test_calc_likelihoods_1():
        from numpy.testing import assert_array_almost_equal
        from numpy.testing import assert_almost_equal

        av_image = np.array([[0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 0],
                             [0, 1, 2, 1, 0],
                             [0, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0]])

        #av_image_error = np.random.normal(0.1, size=av_image.shape)
        av_image_error = 0.1 * av_image

        #nhi_image = av_image + np.random.normal(0.1, size=av_image.shape)
        nhi_image = 5 * av_image

        # add intercept
        intercept_answer = 0.7
        av_image = av_image + intercept_answer

        width_grid = np.arange(0, 10, 0.1)
        dgr_grid = np.arange(0, 2, 0.1)
        intercept_grid = np.arange(-2, 2, 0.1)
        vel_center = 1

        results = \
            cloudpy._calc_likelihoods(
                              nhi_image=nhi_image,
                              av_image=av_image,
                              av_image_error=av_image_error,
                              dgr_grid=dgr_grid,
                              intercept_grid=intercept_grid,
                              )

        dgr_answer = 1/5.0

        assert_almost_equal(results['intercept_max'], intercept_answer)
        assert_almost_equal(results['dgr_max'], dgr_answer)

if 1:
    def test_calc_likelihoods_2():
        from numpy.testing import assert_array_almost_equal
        from numpy.testing import assert_almost_equal
        from myimage_analysis import calculate_nhi
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from mpl_toolkits.axes_grid1 import ImageGrid

        av_image = np.array([[0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 0],
                             [0, 1, 2, 1, 0],
                             [0, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0]])

        #av_image_error = np.random.normal(0.1, size=av_image.shape)
        av_image_error = 0.1 * np.ones(av_image.shape)

        #nhi_image = av_image + np.random.normal(0.1, size=av_image.shape)
        hi_cube = np.zeros((5, 5, 5))

        # make inner channels correlated with av
        hi_cube[:, :, :] = np.array(
            [
             [[  1., 0., 0., 0., 0.],
              [  0., 0., 0., 0., 0.],
              [  0., 0., 0., 0., 0.],
              [  0., 0., 0., 0., 0.],
              [  1., 0., 0., 0., 10.],],

             [[  0., 0., 0., 0., 0.],
              [  0., 0., 2., 0., 0.],
              [  0., 0., 4., 0., 0.],
              [  0., 0., 2., 0., 0.],
              [  0., 0., 0., 0., 0.],],

             [[  0., 0., 0., 0., 0.],
              [  0., 0., 0., 2., 0.],
              [  0., 0., 0., 2., 0.],
              [  0., 0., 0., 2., 0.],
              [  0., 0., 0., 0., 0.],],

             [[  0., 0., 0., 0., 0.],
              [  0., 2., 0., 0., 0.],
              [  0., 2., 0., 0., 0.],
              [  0., 2., 0., 0., 0.],
              [  0., 0., 0., 0., 0.],],

             [[  0., 0., 0., 0., 0.],
              [  0., 0., 0., 0., 0.],
              [  0., 0., 0., 0., 0.],
              [  0., 0., 0., 0., 0.],
              [  1., 0., 0., 0., 0.2],],
             ]
             )

        if 1:
            fig = plt.figure(figsize=(4,4))
            imagegrid = ImageGrid(fig, (1,1,1),
                         nrows_ncols=(1,5),
                         ngrids=5,
                         cbar_mode="single",
                         cbar_location='top',
                         cbar_pad="2%",
                         cbar_size='3%',
                         axes_pad=0.1,
                         aspect=True,
                         label_mode='L',
                         share_all=True)
            cmap = cm.get_cmap('Greys', 5)
            for i in xrange(5):
                im = imagegrid[i].imshow(hi_cube[i, :, :],
                                         origin='lower',
                                         #aspect='auto',
                                         cmap=cmap,
                                         interpolation='none',
                                         vmin=0,
                                         vmax=4)
            #cb = imagegrid[i].cax.colorbar(im)
            cbar = imagegrid.cbar_axes[0].colorbar(im)
            #plt.title('HI Cube')
            plt.savefig('/usr/users/ezbc/Desktop/hi_cube.png')

        # make edge channels totally uncorrelated
        #hi_cube[(0, 4), :, :] = np.arange(0, 25).reshape(5,5)
        #hi_cube[(0, 4), :, :] = - np.ones((5,5))

        hi_vel_axis = np.arange(0, 5, 1)

        # add intercept
        intercept_answer = 0.9
        av_image = av_image + intercept_answer

        if 1:
            fig = plt.figure(figsize=(4,4))
            params = {
              'figure.figsize': (1, 1),
              #'figure.titlesize': font_scale,
             }
            plt.rcParams.update(params)
            imagegrid = ImageGrid(fig, (1,1,1),
                         nrows_ncols=(1,1),
                         ngrids=1,
                         cbar_mode="single",
                         cbar_location='top',
                         cbar_pad="2%",
                         cbar_size='3%',
                         axes_pad=0.1,
                         aspect=True,
                         label_mode='L',
                         share_all=True)
            cmap = cm.get_cmap('Greys', 5)
            im = imagegrid[0].imshow(av_image,
                                         origin='lower',
                                         #aspect='auto',
                                         cmap=cmap,
                                         interpolation='none',
                                         vmin=0,
                                         vmax=4)
            #cb = imagegrid[i].cax.colorbar(im)
            cbar = imagegrid.cbar_axes[0].colorbar(im)
            #plt.title('HI Cube')
            plt.savefig('/usr/users/ezbc/Desktop/av.png')

        width_grid = np.arange(0, 5, 1)
        dgr_grid = np.arange(0, 1, 0.1)
        intercept_grid = np.arange(-1, 2, 0.1)
        vel_center = 2

        results = \
            cloudpy._calc_likelihoods(
                              hi_cube=hi_cube / 1.832e-2,
                              hi_vel_axis=hi_vel_axis,
                              vel_center=vel_center,
                              av_image=av_image,
                              av_image_error=av_image_error,
                              width_grid=width_grid,
                              dgr_grid=dgr_grid,
                              intercept_grid=intercept_grid,
                              )

        dgr_answer = 1/2.0
        width_answer = 2
        width = results['width_max']
        dgr = results['dgr_max']
        intercept = results['intercept_max']
        print width

        if 0:
            width = width_answer
            intercept = intercept_answer
            dgr = dgr_answer

        vel_range = (vel_center - width / 2.0, vel_center + width / 2.0)

        nhi_image = calculate_nhi(cube=hi_cube,
                                  velocity_axis=hi_vel_axis,
                                  velocity_range=vel_range) / 1.823e-2
        if 1:
            fig = plt.figure(figsize=(4,4))
            imagegrid = ImageGrid(fig, (1,1,1),
                         nrows_ncols=(1,1),
                         ngrids=1,
                         cbar_mode="single",
                         cbar_location='top',
                         cbar_pad="2%",
                         cbar_size='3%',
                         axes_pad=0.1,
                         aspect=True,
                         label_mode='L',
                         share_all=True)
            cmap = cm.get_cmap('Greys', 5)
            im = imagegrid[0].imshow(nhi_image,
                                     origin='lower',
                                     #aspect='auto',
                                     cmap=cmap,
                                     interpolation='none',
                                     vmin=0,
                                     vmax=4)
            #cb = imagegrid[i].cax.colorbar(im)
            cbar = imagegrid.cbar_axes[0].colorbar(im)
            #plt.title('HI Cube')
            plt.savefig('/usr/users/ezbc/Desktop/nhi.png')
        if 1:
            fig = plt.figure(figsize=(4,4))
            imagegrid = ImageGrid(fig, (1,1,1),
                         nrows_ncols=(1,1),
                         ngrids=1,
                         cbar_mode="single",
                         cbar_location='top',
                         cbar_pad="2%",
                         cbar_size='3%',
                         axes_pad=0.1,
                         aspect=True,
                         label_mode='L',
                         share_all=True)
            cmap = cm.get_cmap('Greys', 5)
            im = imagegrid[0].imshow(nhi_image * dgr + intercept,
                                         origin='lower',
                                         #aspect='auto',
                                         cmap=cmap,
                                         interpolation='none',
                                         vmin=0,
                                         vmax=4)
            #cb = imagegrid[i].cax.colorbar(im)
            cbar = imagegrid.cbar_axes[0].colorbar(im)
            #plt.title('HI Cube')
            plt.savefig('/usr/users/ezbc/Desktop/av_model.png')

        print('residuals = ')
        print(av_image - (nhi_image * dgr + intercept))
        print('dgr', dgr)
        print('intercept', intercept)
        print('width', width)

        assert_almost_equal(results['intercept_max'], intercept_answer)
        assert_almost_equal(results['dgr_max'], dgr_answer)
        assert_almost_equal(results['width_max'], width_answer)


