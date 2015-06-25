#!/usr/bin/python

def main():

    import cloudpy

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
    figure_dir = \
        '/d/bip3/ezbc/california/figures/'
    av_dir = '/d/bip3/ezbc/california/data/av/'
    hi_dir = '/d/bip3/ezbc/california/data/hi/'
    co_dir = '/d/bip3/ezbc/california/data/co/'
    core_dir = '/d/bip3/ezbc/california/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/california/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    likelihood_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'

    # define filenames
    prop_filename = property_dir + 'california_global_properties.txt'
    av_filename = av_dir + 'california_av_planck_5arcmin.fits'
    av_error_filename = av_dir + 'california_av_planck_5arcmin.fits'
    hi_filename = hi_dir + 'california_hi_galfa_cube_regrid_planckres.fits'
    hi_error_filename = hi_dir + \
            'california_hi_galfa_cube_regrid_planckres_noise.fits'

    cloud = cloudpy.cloud(av_filename,
                          hi_filename,
                          av_error_filename=av_error_filename,
                          hi_error_filename=hi_error_filename,
                          cloud_prop_filename=prop_filename,
                          )


if __name__ == '__main__':
    main()


