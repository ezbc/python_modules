
if 0:
    def nhi_test():

        '''
        operating on following cube:
        [[[  0.   1.]
          [  2.   3.]
          [  4.   5.]]

         [[  6.   7.]
          [  8.   9.]
          [ 10.  11.]]

         [[ 12.  13.]
          [ 14.  15.]
          [ 16.  17.]]

         [[ 18.  19.]
          [ 20.  21.]
          [ 22.  23.]]]

        and velocity axis as
        [0 0.5 1.0 1.5]


        '''

        import numpy as np
        from numpy.testing import assert_almost_equal
        from myimage_analysis import calculate_nhi

        # Create cube
        cube = np.empty((4,3,2))
        count = 0
        for index, element in np.ndenumerate(cube):
            cube[index] = count
            count += 1

        # create velocity axis in km/s
        delta_v = 0.5
        vel_axis = np.arange(0, cube.shape[0]*delta_v, delta_v)

        # Test with cube in 3 dimensions
        # No error returned
        # One vel range
        # --------------------------------------------------------------------------
        nhi_calc = calculate_nhi(cube=cube,
                                 velocity_axis=vel_axis,
                                 velocity_range=(0, 1.0))

        nhi_answer = np.array([[0 + 6 + 12, 1 + 7 + 13],
                               [2 + 8 + 14, 3 + 9 + 15],
                               [4 + 10 + 16, 5 + 11 + 17]],
                              dtype=float)

        nhi_answer *= 1.823e-2 * delta_v

        assert_almost_equal(nhi_calc, nhi_answer)

        # Test with cube in 3 dimensions
        # Error returned
        # One vel range
        # --------------------------------------------------------------------------

        print('3D, Error, 1 vel range')

        noise_cube = 0.1 * np.copy(cube)

        nhi_calc, nhi_error_calc = calculate_nhi(cube=cube,
                                                 velocity_axis=vel_axis,
                                                 velocity_range=(0, 1.0),
                                                 noise_cube=noise_cube,
                                                 return_nhi_error=True)

        nhi_answer = np.array([[0 + 6 + 12, 1 + 7 + 13],
                               [2 + 8 + 14, 3 + 9 + 15],
                               [4 + 10 + 16, 5 + 11 + 17]],
                              dtype=float)

        nhi_answer *= 1.823e-2 * delta_v

        nhi_error_answer = np.array([[0**2 + 0.6**2 + 0.12**2,
                                      0.1**2 + 0.7**2 + 0.13**2],
                                     [0.2**2 + 0.8**2 + 0.14**2,
                                      0.3**2 + 0.9**2 + 0.15**2],
                                     [0.4**2 + 0.10**2 + 0.16**2,
                                      0.5**2 + 0.11**2 + 0.17**2]],
                              dtype=float)**0.5

        nhi_error_answer *= 1.823e-2 * delta_v

        assert_almost_equal(nhi_calc, nhi_answer)
        #assert_almost_equal(nhi_error_calc, nhi_error_answer)

        # Test with cube in 3 dimensions
        # No error returned
        # Image of vel ranges
        # --------------------------------------------------------------------------
        print('3D, no error, multiple vel range')

        velocity_range = np.array([[[0, 0],
                                    [0, 0],
                                    [0.5, 0.5]],
                                   [[1.0, 1.0],
                                    [1.0, 1.0],
                                    [1.5, 1.5]]])

        nhi_calc = calculate_nhi(cube=cube,
                                 velocity_axis=vel_axis,
                                 velocity_range=velocity_range)

        nhi_answer = np.array([[0 + 6 + 12, 1 + 7 + 13],
                               [2 + 8 + 14, 3 + 9 + 15],
                               [10 + 16 + 22, 11 + 17 + 23]],
                              dtype=float)

        nhi_answer *= 1.823e-2 * delta_v

        assert_almost_equal(nhi_calc, nhi_answer)

        # Test with cube in 3 dimensions
        # Error returned
        # Image of vel ranges
        # --------------------------------------------------------------------------
        print('3D, error, multiple vel range')
        velocity_range = np.array([[[0, 0],
                                    [0, 0],
                                    [0.5, 0.5]],
                                   [[1.0, 1.0],
                                    [1.0, 1.0],
                                    [1.5, 1.5]]])

        noise_cube = 0.1 * np.copy(cube)

        nhi_calc, nhi_error = calculate_nhi(cube=cube,
                                 velocity_axis=vel_axis,
                                 velocity_range=velocity_range,
                                 noise_cube=noise_cube,
                                 return_nhi_error=True)

        nhi_answer = np.array([[0 + 6 + 12, 1 + 7 + 13],
                               [2 + 8 + 14, 3 + 9 + 15],
                               [10 + 16 + 22, 11 + 17 + 23]],
                              dtype=float)

        nhi_answer *= 1.823e-2 * delta_v

        assert_almost_equal(nhi_calc, nhi_answer)


        # Test with 2 dimension
        # ==========================================================================
        print('2D, no error, 1 vel range')

        # Create cube
        cube = np.empty((4,6))
        count = 0
        for index, element in np.ndenumerate(cube):
            cube[index] = count
            count += 1

        # create velocity axis in km/s
        delta_v = 0.5
        vel_axis = np.arange(0, cube.shape[0]*delta_v, delta_v)

        # Test with cube in 3 dimensions
        # No error returned
        # One vel range
        # --------------------------------------------------------------------------
        nhi_calc = calculate_nhi(cube=cube,
                                 velocity_axis=vel_axis,
                                 velocity_range=(0, 1.0))

        nhi_answer = np.array([0 + 6 + 12, 1 + 7 + 13,
                               2 + 8 + 14, 3 + 9 + 15,
                               4 + 10 + 16, 5 + 11 + 17],
                              dtype=float)

        nhi_answer *= 1.823e-2 * delta_v

        assert_almost_equal(nhi_calc, nhi_answer)

        # Test with cube in 2 dimensions
        # No error returned
        # Image of vel ranges
        # --------------------------------------------------------------------------
        print('2D, no error, multiple vel range')

        velocity_range = np.array([[0, 0,
                                    0, 0,
                                    0.5, 0.5],
                                   [1.0, 1.0,
                                    1.0, 1.0,
                                    1.5, 1.5]])

        print velocity_range.shape, cube.shape

        nhi_calc = calculate_nhi(cube=cube,
                                 velocity_axis=vel_axis,
                                 velocity_range=velocity_range)

        nhi_answer = np.array([0 + 6 + 12, 1 + 7 + 13,
                               2 + 8 + 14, 3 + 9 + 15,
                               10 + 16 + 22, 11 + 17 + 23],
                              dtype=float)

        nhi_answer *= 1.823e-2 * delta_v

        assert_almost_equal(nhi_calc, nhi_answer)


def test_bin_image():

    import numpy as np
    from numpy.testing import assert_array_almost_equal
    from myimage_analysis import bin_image

    unbinned = np.arange(0, 12, 1).reshape((3, 4))
    #array([[ 0,  1,  2,  3],
    #       [ 4,  5,  6,  7],
    #       [ 8,  9, 10, 11]])

    binned = bin_image(unbinned, binsize=(2, 2), func=np.nansum)

    answer = np.array([[26, 34],
                       ])

    assert_array_almost_equal(binned, answer)

    # bin whole axis
    binned = bin_image(unbinned, binsize=(3, 1), func=np.nansum)

    answer = np.array([[12, 15, 18, 21],
                       ])

    assert_array_almost_equal(binned, answer)

    # Bin
    binned = bin_image(unbinned, binsize=(1, 2), func=np.nansum)

    answer = np.array([[1, 5],
                       [9, 13],
                       [17, 21]])

    assert_array_almost_equal(binned, answer)

    # Bin
    if 0:
        def func(image, axis=None):
            return np.nansum(image, axis=axis)

        binned = bin_image(unbinned, binsize=(1, 2), func=func)

        answer = np.array([[1, 13],
                           [41, 85],
                           [145, 221]])

        assert_array_almost_equal(binned, answer)

    m = np.arange(0,100,1).reshape((10,10))
    n = bin_image(m, binsize=(2,2))

    answer = np.array([[ 22,  30,  38,  46,  54],
                       [102, 110, 118, 126, 134],
                       [182, 190, 198, 206, 214],
                       [262, 270, 278, 286, 294],
                       [342, 350, 358, 366, 374]])
    assert_array_almost_equal(n, answer)


if 0:
    def main():
        import pyfits as fits
        from myimage_analysis import bin_image
        import matplotlib.pyplot as plt
        import numpy as np

        av_dir = '/d/bip3/ezbc/perseus/data/av/'

        av_data, av_header = fits.getdata(av_dir + \
                                'perseus_av_planck_5arcmin.fits',
                                          header=True)

        av_data[av_data > 1] = np.nan

        av_data_binned, av_header_binned = bin_image(av_data,
                                               binsize=(11, 11),
                                               func=np.nanmean,
                                               header=av_header)

        print av_data.shape, av_data_binned.shape

        fits.writeto(av_dir + 'test.fits', av_data_binned, av_header_binned,
                clobber=True)

        if 0:
            plt.imshow(av_data, origin='lower left')
            plt.show()
            plt.imshow(av_data_binned, origin='lower left')
            plt.show()

        hi_dir = '/d/bip3/ezbc/perseus/data/hi/'

        hi_data, hi_header = fits.getdata(hi_dir + \
                    'perseus_hi_galfa_cube_regrid_planckres.fits',
                header=True)

        hi_data[hi_data > 10] = np.nan

        hi_data_binned, hi_header_binned = bin_image(hi_data,
                                               binsize=(11, 11),
                                               func=np.nanmean,
                                               header=hi_header)
        print hi_data.shape, hi_data_binned.shape
        fits.writeto(hi_dir + 'test.fits', hi_data_binned, hi_header_binned,
                clobber=True)

        plt.imshow(hi_data[500], origin='lower left')
        plt.show()
        plt.imshow(hi_data_binned[500], origin='lower left')
        plt.show()

if __name__ == '__main__':
    main()

