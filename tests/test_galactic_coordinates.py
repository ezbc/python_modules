#!/usr/bin/python

import planckpy as pl
reload(pl)
from astropy.io import fits
import pyfits as pf
import os

def plot_av_image(av_image=None, header=None, contour_image=None,
        contour_header=None, title=None, limits=None, contours=None,
        boxes=False, savedir='./', filename=None, show=True,
        colorbar_lims=None):

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
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    fontScale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (6, 6),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # Create axes
    grid_helper = wcs.GridHelperSky(wcs=header)

    if 1:
        imagegrid = ImageGrid(fig, (1,1,1),
                     nrows_ncols=(1,1),
                     ngrids=1,
                     cbar_mode="single",
                     cbar_location='right',
                     cbar_pad="2%",
                     cbar_size='3%',
                     axes_pad=0,
                     axes_class=(wcs.Axes,
                                 dict(grid_helper=grid_helper)),
                     aspect=False,
                     label_mode='L',
                     share_all=False)

    # create axes
    #ax = wcs.subplot(111, grid_helper=grid_helper)
    ax = imagegrid[0]
    cmap = cm.winter # colormap
    # show the image
    im = ax.imshow(av_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,)

    # Contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Asthetics
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)
    if title is not None:
        ax.set_title(title)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    if colorbar_lims is not None:
        im.set_clim(colorbar_lims[0], colorbar_lims[1])

    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Plot Av contours
    #if contour_image is not None:
    #    ax.contour(contour_image, levels=contours, colors='r')

    # Write label to colorbar
    cb.set_label_text(r'A$_V$ (mag)',)

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def write_test_data():

    def extract_test_data():

        ''' Creates test data of the perseus molecular cloud
        '''

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

        dec_range = (21.3, 30.3)
        ra_range = (60.0, 73.0)

        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 0,
                dr_version = 2,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)

        return (data, header)

    def ebv2av(data, header):
        # Av=E(B-V) x 3.05263, 152 in Tielens' ISM book.

        data *= 3.05363
        header['BUNIT'] = 'mag'

        return (data, header)

    def tau353_to_ebv(data, header):

        '''

        The dust emission maps are derived from:
        Planck 2013 results. XI. All-sky model of thermal dust emission

        Whereby they fit for the dust optical depth at 353 GHz, and then calibrated
        the conversion of the optical depth to color excess using background
        quasars.

        E(B-V) can be calculated from the dust optical depth at 353 GHz by
        (E(B - V)/tau 353 = 1.49 x 10**4

        '''

        data *= 1.49e4
        header['BUNIT'] = 'mag'

        return (data, header)

    def write_data(data, header, filename):

        fits.writeto(filename,
                data,
                header = header,
                clobber = True,
                output_verify = 'fix')

    def regrid_images(image, template):
        import mirpy

        # Avoiding the miriad errors avoids errors from existing files
        try:
            mirpy.fits(image+'.fits',
                    out=image+'.mir',
                    op='xyin')
        except mirpy.wrapper.MiriadError:
            pass

        try:
            mirpy.fits(template+'.fits',
                    out=template+'.mir',
                    op='xyin')
        except mirpy.wrapper.MiriadError:
            pass

        try:
            mirpy.regrid(image+'.mir', tin=template+'.mir',
                    out=image+'_regrid.mir')
        except mirpy.wrapper.MiriadError:
            pass

        try:
            mirpy.fits(image+'_regrid.mir', out=image+'_regrid.fits',
                    op='xyout')
        except mirpy.wrapper.MiriadError:
            pass

        hdu = pf.open(image+'_regrid.fits')
        data, header = hdu[0].data, hdu[0].header

        return data, header

    test_dir = '/d/bip3/ezbc/perseus/data/planck/test_data/'

    # Write out av data
    (data, header) = extract_test_data()
    (data, header) = tau353_to_ebv(data, header)
    (data, header) = ebv2av(data, header)
    write_data(data, header, test_dir + 'test_perseus_av_planck.fits')

    header_file = pf.Header(header)
    header = pf.Header(header)

    # Write the header
    header_file.tofile(test_dir + 'test_perseus_av_planck.header',
                   clobber=True)

    # Compare map with Jouni's K09
    av_image, av_header = pf.getdata(test_dir + \
            'perseus_av_kainulainen2009_nan.fits', header=True)

    av_hdu = pf.open(test_dir + 'perseus_av_kainulainen2009_nan.fits')
    av_data, av_header = av_hdu[0].data, av_hdu[0].header

    #
    os.chdir(test_dir)

    planck_data, planck_header = regrid_images('test_perseus_av_planck',
            'perseus_av_kainulainen2009_nan')

    plot_av_image(av_image=av_image,
            header=av_header,
            contour_image=planck_data,
            contour_header=planck_header,
            contours=[3,7,11],
            colorbar_lims=(-3,14),
            title='perseus: Kainulainen 09 map with Planck Contours',
            savedir=test_dir,
            filename='perseus_k09_planck_map.png',
            show=False)

def main():
    #write_test_data()
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
    test_dir = '/d/bip3/ezbc/planck/tests/'

    (data, header) = pl.get_data(data_type='857',
                                 data_location=data_location,
                                 coord_type='galactic',
                                 x_range=(280, 300),
                                 y_range=(-40, -20))

    fits.writeto(test_dir + 'LMC_857ghz.fits',
                data,
                header = header,
                clobber = True,
                output_verify = 'fix')

    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    l_range = (155, 165)
    b_range = (-25, -10)

    (data, header) = pl.get_data(data_location = data_location,
            data_type = 'CO-Type3',
            x_range = l_range,
            y_range = b_range,
            coord_type = 'galactic',
            field = 0,
            dr_version =1,
            resolution = 0.1,
            cut_last_pixel = False,
            verbose = True)

    fits.writeto(test_dir + 'perseus_galactic_ebv.fits',
                data,
                header = header,
                clobber = True,
                output_verify = 'fix')



if __name__ == '__main__':
    main()
