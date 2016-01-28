#!/usr/bin/python

import numpy as np
import gausspy.gp as gp
import traingp as train
import pickle

def construct_spectrum(fit_params, vel_axis):

    import gausspy.AGD_decomposer as agd

    spectrum = np.zeros(len(vel_axis))
    ncomps = len(fit_params) / 3
    amps = fit_params[0:ncomps]
    fwhms = fit_params[ncomps:ncomps*2]
    means = fit_params[ncomps*2:ncomps*3]
    for j in xrange(ncomps):
        comp_func = agd.gaussian(amps[j], fwhms[j], means[j])

        spectrum += comp_func(vel_axis)

    return spectrum

def decompose_data(filename_data, g_train=None, filename_decomposed=None,
        data_dict=None):

    import gausspy.AGD_decomposer as agd

    g = gp.GaussianDecomposer()

    print('\nDecomposing data...')

    #Two phase
    if g_train is not None:
        g.set('alpha1', g_train.p['alpha1'])
        g.set('alpha2', g_train.p['alpha2'])
        g.set('phase', g_train.p['phase'])
        g.set('SNR_thresh', g_train.p['SNR_thresh'])
        g.set('SNR2_thresh', g_train.p['SNR2_thresh'])
    else:
        g.set('alpha1', 2.5)
        g.set('alpha2', 6)
        g.set('BLFrac', 0.02)
        g.set('phase', 'two')
        g.set('SNR_thresh', 3.)
        g.set('SNR2_thresh', 3.)
    g.set('mode', 'conv')
    g.set('verbose', False)

    if data_dict is None:
        new_data = g.batch_decomposition(filename_data)

        if filename_decomposed is not None:
            pickle.dump(new_data, open(filename_decomposed, 'w'))
    else:
        results_dict = {}
        #results_dict['spectra'] = []
        results_dict['results'] = []

        #if filename_decomposed is not None:
        #    results_dict = pickle.load(open(filename_decomposed, 'r'))

        x_values = data_dict['velocity_axis']

        for i in xrange(len(data_dict['data_list'])):
        #for i in xrange(12274, 12276):
            print('\n\titeration ' + str(i))
            try:
                results = g.decompose(x_values,
                                      data_dict['data_list'][i],
                                      data_dict['errors'])
            except (np.linalg.LinAlgError, ValueError):
                results['N_components'] = 0

            # record location of spectrum
            results['spectrum_number'] = i

            if 0:
                # Construct spectrum
                if results['N_components'] > 0:
                    spectrum = \
                        construct_spectrum(results['best_fit_parameters'],
                                           x_values)
                else:
                    spectrum = np.zeros(len(x_values))

                # Plot scratch plot of fits
                import matplotlib.pyplot as plt
                plt.close(); plt.clf()
                plt.plot(x_values,
                         data_dict['data_list'][i])
                plt.plot(x_values,
                         spectrum,
                         alpha=0.5,
                         linewidth=3)
                plt.savefig('/d/bip3/ezbc/scratch/spectrum_fit_' + \
                            str(i) + '.png')

            #results_dict['spectra'].append(spectrum)
            results_dict['results'].append(results)

        # Add positions to results
        results_dict['positions'] = data_dict['positions'].copy()

        # Add velocity axis to results
        results_dict['velocity_axis'] = data_dict['velocity_axis'].copy()

        if filename_decomposed is not None:
            pickle.dump(results_dict, open(filename_decomposed, 'w'))

        return results_dict

def get_decomposed_data(filename_data, g_train=None,
        filename_decomposed=None, data_dict=None, load=False,):

    import os
    import pickle

    # load decomposed data if exists, else perform decomposition
    if load:
        if os.path.isfile(filename_decomposed):
            data_decomp = pickle.load(open(filename_decomposed, 'r'))
            perform_decomposition = False
        else:
            perform_decomposition = True
    else:
        perform_decomposition = True

    # Run AGD on data?
    if perform_decomposition:
        data_decomp = decompose_data(filename_data,
                                     g_train=g_train,
                                     data_dict=data_dict,
                                     filename_decomposed=filename_decomposed,
                                     )

    return data_decomp

def perform_PCA(data, n_components=3, pca_type='regular'):

    if pca_type == 'regular':
        from sklearn.decomposition import PCA

        pca = PCA(n_components=n_components, whiten=True)

        data_reduced = pca.fit_transform(data)

        #print('\tExplained variance:')
        #print('\t', pca.explained_variance_ratio_)
    elif pca_type == 'linear':
        from sklearn.decomposition import KernelPCA

        pca = KernelPCA(n_components=n_components,
                        kernel='linear',
                        #fit_inverse_transform=True,
                        #eigen_solver='arpack',
                        )

        data_reduced = pca.fit_transform(data)

    return data_reduced

def cluster_data(data, n_clusters=2, method='kmeans'):

    ''' Clusters data.

    Parameters
    ----------
    n_cluster : int
        Number of clusters
    method : str
        Either kmeans or spectral

    '''

    from sklearn.cluster import KMeans, SpectralClustering, DBSCAN

    # Initialize the clustering method instance
    if method == 'kmeans':
        estimator = KMeans(n_clusters=n_clusters)
    elif method == 'spectral':
        estimator = SpectralClustering(n_clusters=n_clusters,
                                       #eigen_solver='arpack',
                                       affinity="nearest_neighbors"
                                       #affinity="rbf"
                                       )
    elif method == 'dbscan':
        estimator = DBSCAN(min_samples=10,
                           eps=10,
                           )

    # Fit the data
    estimator.fit(data)

    labels = estimator.labels_

    colors = labels.astype(np.float)
    colors /= colors.max()

    return colors

def create_synthetic_cube(results_dict=None, cube_data=None,
        pix_positions=None, velocity_axis=None, fit_params_list=None,):

    if results_dict is not None:
        # Create cube based on number of pixels
        xy_pix = results_dict['positions']['pix']
        z_pix = np.arange(0, len(results_dict['velocity_axis'])+1)
        shape = (np.max(z_pix), np.max(xy_pix[:,0]), np.max(xy_pix[:, 1]),)
        cube = np.zeros(shape)

        shape = (np.max(z_pix), cube_data.shape[1], cube_data.shape[2])
        cube = np.zeros(shape)

        n_spectra = len(results_dict['results'])

        # Add a spectrum to each xy pixel in cube
        for i in xrange(n_spectra):
            if results_dict['results'][i]['N_components'] > 0:
                result = results_dict['results'][i]
                fit_params = result['best_fit_parameters']

                # Get pixel positions
                x, y = results_dict['positions']['pix'][i]

                # If any gaussians present add them to the cube
                if result['N_components'] > 0:
                    spectrum = construct_spectrum(result['best_fit_parameters'],
                                                  results_dict['velocity_axis'])

                    cube[:, int(x), int(y)] = spectrum
    else:
        # Create cube based on number of pixels
        xy_pix = pix_positions
        z_pix = np.arange(0, len(velocity_axis)+1)

        shape = (np.max(z_pix), cube_data.shape[1], cube_data.shape[2])
        cube = np.zeros(shape)

        n_spectra = len(fit_params_list)

        # Add a spectrum to each xy pixel in cube
        for i in xrange(n_spectra):
            # Get pixel positions
            x, y = xy_pix[i]

            fit_params = fit_params_list[i]

            # add Gaussians to cube
            spectrum = construct_spectrum(fit_params,
                                          velocity_axis)

            cube[:, int(x), int(y)] = spectrum

    return cube

def plot_nhi_maps(nhi_1, nhi_2, header=None, contour_image=None,
        limits=None, contours=None, filename=None, show=False,
        nhi_1_vlimits=[None,None], nhi_2_vlimits=[None,None], vscale='linear'):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    import astropy.io as fits
    import matplotlib.pyplot as plt
    import myplotting as myplt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1 import AxesGrid

    # Import external modules
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import AxesGrid
    import pywcsgrid2 as wcs
    from matplotlib.patches import Polygon
    import matplotlib.patheffects as PathEffects
    import myplotting as myplt

    # Set up plot aesthetics
    plt.clf(); plt.close()

    # Color map
    cmap = plt.cm.gnuplot
    #cmap = myplt.reverse_colormap(plt.cm.copper)
    cmap = plt.cm.copper_r
    cmap = plt.cm.gist_heat_r

    # Create figure instance
    fig = plt.figure(figsize=(3, 4))

    if nhi_2 is not None:
        nrows_ncols=(1,2)
        ngrids=2
    else:
        nrows_ncols=(1,1)
        ngrids=1

    #grid_helper = wcs.GridHelper(wcs=header)
    axes = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="3%",
                 cbar_size='7%',
                 axes_pad=0.1,
                 axes_class=(wcs.Axes,
                             #dict(grid_helper=grid_helper)),
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    if vscale == 'log':
        norm = matplotlib.colors.LogNorm()
        nhi_1[nhi_1 <= 0] = np.nan
        nhi_2[nhi_2 <= 0] = np.nan
        if 0 in nhi_1_vlimits:
            nhi_1_vlimits = [None, None]
        if 0 in nhi_2_vlimits:
            nhi_2_vlimits = [None, None]
    else:
        norm = None

    # ------------------
    # NHI image
    # ------------------
    # create axes
    ax = axes[0]

    # show the image
    im = ax.imshow(nhi_1,
            interpolation='none',
            origin='lower',
            cmap=cmap,
            vmin=nhi_1_vlimits[0],
            vmax=nhi_1_vlimits[1],
            norm=norm,
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("gal")
    #ax.set_ticklabel_type("hms", "dms")
    ax.set_ticklabel_type("absdeg", "absdeg")

    #ax.set_xlabel('Right Ascension [J2000]',)
    #ax.set_ylabel('Declination [J2000]',)
    ax.set_xlabel('Glong',)
    ax.set_ylabel('Glat',)


    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    #ax.tick_params(colors='w')
    ax.locator_params(nbins=4)
    #ax.tick_params(colors='w')

    # plot limits
    if limits is not None:
        limits_pix = myplt.convert_wcs_limits(limits,
                                              header,
                                              frame='galactic')
        ax.set_xlim(limits_pix[0],limits_pix[1])
        ax.set_ylim(limits_pix[2],limits_pix[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Write label to colorbar
    cb.set_label_text(r'$N$(H\textsc{i}) [10$^{20}$ cm$^{-2}$]',)

    # ------------------
    # Av image
    # ------------------
    if nhi_2 is not None:
        # create axes
        ax = axes[1]
        # show the image
        im = ax.imshow(nhi_2,
                interpolation='none',
                origin='lower',
                cmap=cmap,
                vmin=nhi_2_vlimits[0],
                vmax=nhi_2_vlimits[1],
                norm=norm
                )

        # Asthetics
        ax.set_display_coord_system("gal")
        #ax.set_ticklabel_type("hms", "dms")
        ax.set_ticklabel_type("absdeg", "absdeg")

        #ax.set_xlabel('Right Ascension [J2000]',)
        #ax.set_ylabel('Declination [J2000]',)
        ax.set_xlabel('Glong',)
        ax.set_ylabel('Glat',)

        ax.locator_params(nbins=4)
        #ax.tick_params(colors='w')

        # colorbar
        cb = ax.cax.colorbar(im)
        cmap.set_bad(color='w')

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits,
                                                  header,
                                                  frame='galactic')
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        ax.tick_params(axis='xy', which='major', colors='w')

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        cb.set_label_text(r'$N$(H\textsc{i}) [10$^{20}$ cm$^{-2}$]',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

def plot_nhi_map_panels(nhi_list, header=None, contour_image=None,
        limits=None, contours=None, filename=None, show=False,
        nhi_vlimits=[None,None], vscale='linear'):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1 import AxesGrid

    # Set up plot aesthetics
    plt.clf(); plt.close()

    # Color map
    cmap = plt.cm.gnuplot
    #cmap = myplt.reverse_colormap(plt.cm.copper)
    cmap = plt.cm.copper
    cmap = plt.cm.gist_heat_r

    #
    ngrids = len(nhi_list)
    nrows, ncols = myplt.get_square_grid_sides(ngrids)
    nrows, ncols = 1, len(nhi_list)
    nrows_ncols = (nrows, ncols)

    # Create figure instance
    fig = plt.figure(figsize=(3 * ncols + 1, 3 * nrows + 1))

    #grid_helper = wcs.GridHelper(wcs=header)
    axes = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 #cbar_mode="each",
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='7%',
                 axes_pad=0.3,
                 #axes_class=(wcs.Axes,
                 #            dict(grid_helper=grid_helper)),
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)


    for i in xrange(ngrids):
        # create axes
        ax = axes[i]

        if vscale == 'log':
            norm = matplotlib.colors.LogNorm()
            nhi_list[i][nhi_list[i] <= 0] = np.nan
            if 0 in nhi_vlimits:
                nhi_vlimits = [None, None]
        else:
            norm = None

        # show the image
        im = ax.imshow(nhi_list[i],
                       #interpolation='nearest',
                       origin='lower',
                       cmap=cmap,
                       vmin=nhi_vlimits[0],
                       vmax=nhi_vlimits[1],
                       norm=norm,
                       #norm=matplotlib.colors.LogNorm()
                       )

        # Asthetics
        ax.set_display_coord_system("gal")
        #ax.set_ticklabel_type("hms", "dms")
        ax.set_ticklabel_type("absdeg", "absdeg")

        #ax.set_xlabel('Right Ascension [J2000]',)
        #ax.set_ylabel('Declination [J2000]',)
        ax.set_xlabel(r'$l$ [deg]',)
        ax.set_ylabel(r'$b$ [deg]',)
        ax.tick_params(axis='xy', which='major', colors='w')
        ax.locator_params(nbins=4)
        #ax.tick_params(colors='w')

        # colorbar
        cb = ax.cax.colorbar(im)
        cmap.set_bad(color='w')

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits,
                                                  header,
                                                  frame='galactic')
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        cb.set_label_text(r'$N$(H\textsc{i}) [10$^{20}$ cm$^{-2}$]',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

def plot_vel_map_panels(vel_list, header=None, contour_image=None,
        limits=None, contours=None, filename=None, show=False,
        vel_vlimits=[None,None], vscale='linear'):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1 import AxesGrid

    # Set up plot aesthetics
    plt.clf(); plt.close()

    # Color map
    cmap = plt.cm.gnuplot
    #cmap = myplt.reverse_colormap(plt.cm.copper)
    cmap = plt.cm.BrBG
    cmap = plt.cm.winter

    #
    ngrids = len(vel_list)
    nrows, ncols = myplt.get_square_grid_sides(ngrids)
    nrows, ncols = 1, len(vel_list)
    nrows_ncols = (nrows, ncols)

    # Create figure instance
    fig = plt.figure(figsize=(3 * ncols + 1, 3 * nrows + 1))

    #grid_helper = wcs.GridHelper(wcs=header)
    axes = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='7%',
                 axes_pad=0.3,
                 #axes_class=(wcs.Axes,
                 #            dict(grid_helper=grid_helper)),
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)


    for i in xrange(ngrids):
        # create axes
        ax = axes[i]

        if vscale == 'log':
            norm = matplotlib.colors.LogNorm()
            vel_list[i][vel_list[i] <= 0] = np.nan
            if 0 in vel_vlimits:
                vel_vlimits = [None, None]
        else:
            norm = None

        # show the image
        im = ax.imshow(vel_list[i],
                       #interpolation='nearest',
                       origin='lower',
                       cmap=cmap,
                       vmin=vel_vlimits[0],
                       vmax=vel_vlimits[1],
                       norm=norm,
                       #norm=matplotlib.colors.LogNorm()
                       )

        # Asthetics
        ax.set_display_coord_system("gal")
        #ax.set_ticklabel_type("hms", "dms")
        ax.set_ticklabel_type("absdeg", "absdeg")

        #ax.set_xlabel('Right Ascension [J2000]',)
        #ax.set_ylabel('Declination [J2000]',)
        ax.set_xlabel(r'$l$ [deg]',)
        ax.set_ylabel(r'$b$ [deg]',)
        ax.tick_params(axis='xy', which='major', colors='w')
        ax.locator_params(nbins=4)
        #ax.tick_params(colors='w')

        # colorbar
        cb = ax.cax.colorbar(im)
        cmap.set_bad(color='w')

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits,
                                                  header,
                                                  frame='galactic')
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        cb.set_label_text(r'Velocity [km/s]',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

