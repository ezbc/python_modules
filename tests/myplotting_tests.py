import numpy as np
import matplotlib.pyplot as plt
import myplotting as myplt
from matplotlib.ticker import ScalarFormatter

if 0:
    def test_textoverlap_simple():
        plt.close(); plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(range(10), range(10))
        testing = True
        if testing:
        # Mess things up
            fig.canvas.draw()
            labels = [item.get_text().zfill(10) for item in ax.get_xticklabels()]
            ax.set_xticklabels(labels)

        fig.canvas.draw()

        ax = myplt.delete_overlapping_xlabels3(fig,ax)

        #plt.show()

    def test_textoverlap_minorticks():

        from matplotlib.ticker import MultipleLocator, FormatStrFormatter
        minorLocator   = MultipleLocator(5)



        plt.close(); plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.plot(range(100), range(100))
        testing = True
        # Mess things up
        if testing:
            fig.canvas.draw()
            labels = [item.get_text().zfill(3) for item in \
                    ax.get_xticklabels(minor=True)]
            ax.set_xticklabels(labels, minor=True)
            fig.canvas.draw()
            #ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.set_minor_formatter(FormatStrFormatter('000000'))
            #plt.setp(ax.get_xticklabels(minor=True), visible=True)

        fig.canvas.draw()

        #ax = myplt.delete_overlapping_xlabels3(fig,ax)

        plt.show()

    def test_bbox():

        import numpy as np
        import matplotlib.pyplot as plt


        plt.close(); plt.clf()
        fig, ax = plt.subplots(1)
        ax.plot(range(10), range(10))
        fig.canvas.draw()
        bbox1 = ax.get_xticklabels()[0].get_window_extent()
        bbox2 = ax.get_xticklabels()[1].get_window_extent()
        print('Overlap before nan values', bbox1.overlaps(bbox2))
        overlap = bbox2.set_points(np.array([[np.nan, np.nan], [np.nan, np.nan]]))
        print('Overlap after nan values', bbox1.overlaps(bbox2))
        assert not overlap

        '''
            import numpy as np
            import matplotlib.pyplot as plt


            plt.close(); plt.clf()
            fig, ax = plt.subplots(1)
            ax.plot(range(10), range(10))
            fig.canvas.draw()
            bbox1 = ax.get_xticklabels()[0].get_window_extent()
            bbox2 = ax.get_xticklabels()[1].get_window_extent()
            print('Overlap before nan values', bbox1.overlaps(bbox2))
            overlap = bbox2.set_points(np.array([[np.nan, np.nan], [np.nan, np.nan]]))
            print('Overlap after nan values', bbox1.overlaps(bbox2))
        '''



    def test_triangle_dist():

        import numpy as np
        import myplotting as myplt

        dists = np.zeros((4,5,6))
        dists[0,1,2] = 0.1
        dists[1,1,2] = 0.1
        dists[2,1,2] = 0.1
        dists[3,4,4] = 0.1
        dists[3,3,4] = 0.1

        myplt.corner_plot(dists, filename='test_plots/test_corner_plot.png')

if 1:
    def test_scatter_contour_log():

        import myimage_analysis as myia

        # Parameters
        # ----------
        levels = (0.99, 0.985, 0.7)
        levels = (0.999, 0.998, 0.96, 0.86, 0.58,)
        levels = 7
        levels = np.logspace(np.log10(0.995), np.log10(0.50), 5)
        log_counts = 0
        limits = [1, 10, -3, 30]
        limits = None

        x = np.linspace(1, 100, 10000)
        y = x**3 + np.random.normal(10, size=10000) + 100
        y_real = x**3 + 100

        x_nonans, y_nonans = myia.mask_nans((x,y))
        scale = ['linear', 'log']

        fig, ax = plt.subplots()

        if limits is None:
            xmin = np.min(x_nonans)
            ymin = np.min(y_nonans)
            xmax = np.max(x_nonans)
            ymax = np.max(y_nonans)
            if scale[0] == 'linear':
                xscalar = 0.25 * xmax
                xlims = xmin - xscalar, xmax + xscalar
            elif scale[0] == 'log':
                xscalar = 0.25
                xlims = xmin * xscalar, xmax / xscalar
            if scale[1] == 'linear':
                yscalar = 0.25 * ymax
                ylims = ymin - yscalar, ymax + yscalar
            elif scale[1] == 'log':
                yscalar = 0.25
                ylims = ymin * yscalar, ymax / yscalar

            limits = [xlims[0], xlims[1], ylims[0], ylims[1]]

        if 1:
            bins = [30, np.logspace(np.log10(limits[2]),
                                    np.log10(limits[3]),
                                    30),
                    ]
        else:
            bins = 40

        print bins

        contour_range = ((limits[0], limits[1]),
                         (limits[2], limits[3]))

        cmap = myplt.truncate_colormap(plt.cm.binary, 0.2, 1, 1000)

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        l1 = myplt.scatter_contour(x_nonans.ravel(),
                             y_nonans.ravel(),
                             threshold=3,
                             log_counts=log_counts,
                             levels=levels,
                             ax=ax,
                             extent=limits,
                             histogram2d_args=dict(bins=bins,
                                                   range=contour_range,
                                                   ),
                             plot_args=dict(marker='o',
                                            linestyle='none',
                                            color='black',
                                            alpha=0.3,
                                            markersize=2),
                             contour_args=dict(
                                               #cmap=plt.cm.binary,
                                               cmap=cmap,
                                               #cmap=cmap,
                                               ),
                             )

        ax.plot(x, y_real, color='r', alpha=0.5, linewidth=2)

        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'x')
        ax.set_ylabel(r'y')
        #ax.set_title(core_names[i])
        ax.legend(loc='lower right')

        plt.savefig('test_plots/test_scatter_contour_log.png')

    def test_scatter_contour():

        from astropy.io import fits
        from myimage_analysis import calculate_nhi
        import mygeometry as myg
        from mycoords import make_velocity_axis


        # Parameters
        # ----------
        levels = (0.99, 0.985, 0.7)
        levels = (0.999, 0.998, 0.96, 0.86, 0.58,)
        levels = 7
        levels = np.logspace(np.log10(0.995), np.log10(0.50), 5)
        log_counts = 0
        limits = [1, 10, -3, 30]
        limits = None

        # Begin test
        # ----------
        data_dir = '/d/bip3/ezbc/perseus/data/'
        av = fits.getdata(data_dir + 'av/perseus_av_planck_tau353_5arcmin.fits')
        hi, hi_header = fits.getdata(data_dir + \
                          'hi/perseus_hi_galfa_cube_regrid_planckres.fits',
                          header=True)

        hi_vel_axis = make_velocity_axis(hi_header)

        nhi = calculate_nhi(cube=hi,
                            velocity_axis=hi_vel_axis,
                            velocity_range=[0, 10],
                            )

        # Drop the NaNs from the images
        indices = np.where((av == av) &\
                           (nhi == nhi)
                           )

        av_nonans = av[indices]
        nhi_nonans = nhi[indices]

        fig, ax = plt.subplots()

        if limits is None:
            xmin = np.min(nhi_nonans)
            ymin = np.min(av_nonans)
            xmax = np.max(nhi_nonans)
            ymax = np.max(av_nonans)
            xscalar = 0.25 * xmax
            yscalar = 0.25 * ymax
            limits = [xmin - xscalar, xmax + xscalar,
                      ymin - yscalar, ymax + yscalar]

        contour_range = ((limits[0], limits[1]),
                         (limits[2], limits[3]))

        cmap = myplt.truncate_colormap(plt.cm.binary, 0.2, 1, 1000)

        l1 = myplt.scatter_contour(nhi_nonans.ravel(),
                             av_nonans.ravel(),
                             threshold=3,
                             log_counts=log_counts,
                             levels=levels,
                             ax=ax,
                             histogram2d_args=dict(bins=30,
                                                   range=contour_range),
                             plot_args=dict(marker='o',
                                            linestyle='none',
                                            color='black',
                                            alpha=0.3,
                                            markersize=2),
                             contour_args=dict(
                                               #cmap=plt.cm.binary,
                                               cmap=cmap,
                                               #cmap=cmap,
                                               ),
                             )

        scale = ['linear', 'linear']
        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$N($H$\textsc{i}) \times\,10^{20}$ cm$^{-2}$')
        ax.set_ylabel(r'$A_V$ [mag]')
        #ax.set_title(core_names[i])
        ax.legend(loc='lower right')

        plt.savefig('test_plots/test_scatter_contour.png')

