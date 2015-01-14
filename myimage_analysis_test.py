
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

