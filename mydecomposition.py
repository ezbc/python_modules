#!/usr/bin/python

import numpy as np

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
    g.set('verbose', True)

    if data_dict is None:
        new_data = g.batch_decomposition(filename_data)

        if filename_decomposed is not None:
            pickle.dump(new_data, open(filename_decomposed, 'w'))
    else:
        results_dict = {}
        #results_dict['spectra'] = []
        results_dict['results'] = []

        if filename_decomposed is not None:
            results_dict = pickle.load(open(filename_decomposed, 'r'))

        for i in xrange(len(data_dict['data_list'])):
            print '\n\titeration ' + str(i)
            try:
                results = g.decompose(data_dict['x_values'][0],
                                      data_dict['data_list'][i],
                                      data_dict['errors'][i])
            except np.linalg.LinAlgError:
                results['N_components'] = 0

            # record location of spectrum
            results['spectrum_number'] = i

            # Construct spectrum
            if results['N_components'] > 0:
                spectrum = construct_spectrum(results['best_fit_parameters'],
                                              data_dict['x_values'][0])
            else:
                spectrum = np.zeros(len(data_dict['x_values'][0]))

            # Plot scratch plot of fits
            if 0:
                import matplotlib.pyplot as plt
                plt.close(); plt.clf()
                plt.plot(data_dict['x_values'][0],
                         data_dict['data_list'][i])
                plt.plot(data_dict['x_values'][0],
                         spectrum,
                         alpha=0.5,
                         linewidth=3)
                plt.savefig('/d/bip3/ezbc/scratch/spectrum_fit_' + \
                            str(i) + '.png')

            #results_dict['spectra'].append(spectrum)
            results_dict['results'].append(results)

        # Add positions to results
        results_dict['positions'] = data_dict['positions'].copy()

        if filename_decomposed is not None:
            pickle.dump(results_dict, open(filename_decomposed, 'w'))

        return result_dict

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

def perform_PCA(data, cluster=False):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn import datasets
    from sklearn.decomposition import PCA

    for i in xrange(data.shape[1]):
        print np.sort(data[:, i])

    n_components = 2
    pca = PCA(n_components=n_components)

    X_reduced = pca.fit_transform(data)
    #X_reduced = pca.fit(data)
    #print pca.score_samples(data)

    if cluster:
        colors = cluster_data(X_reduced, n_clusters=5)

    print colors

    if n_components <= 100:
        from myplotting import scatter_contour
        fig = plt.figure(1, figsize=(4,4))
        ax = fig.add_subplot(111)
        #X_reduced[:,0] += np.abs(np.min(X_reduced[:,0]) + 1.0)

        ax.scatter(X_reduced[:, 0], X_reduced[:, 1],
                   color=colors,
                   alpha=0.1)
        #ax.set_xscale('log')
        if 0:
            scatter_contour(X_reduced[:, 0],
                                    X_reduced[:,1],
                                 threshold=2,
                                 log_counts=0,
                                 levels=5,
                                 ax=ax,
                                 histogram2d_args=dict(bins=50,),
                                 plot_args=dict(marker='o',
                                                linestyle='none',
                                                markeredgewidth=0,
                                                color='black',
                                                alpha=0.4,
                                                markersize=2.5,
                                                ),
                                 contour_args=dict(cmap=plt.cm.binary,)
                                 )

    elif n_components > 3:
        plt.close(); plt.clf()
        if 0:
            fig = plt.figure(1, figsize=(4,4))

            # plot the first three PCA dimensions
            ax = Axes3D(fig, elev=-150, azim=50)
            ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2],
                       cmap=plt.cm.Paired)
            ax.set_title("First three PCA directions")
            ax.set_xlabel("1st eigenvector")
            ax.w_xaxis.set_ticklabels([])
            ax.set_ylabel("2nd eigenvector")
            ax.w_yaxis.set_ticklabels([])
            ax.set_zlabel("3rd eigenvector")
            ax.w_zaxis.set_ticklabels([])

        # Clustered data
        fig = plt.figure(1, figsize=(4,4))
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
        ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=colors)
        #ax.w_xaxis.set_ticklabels([])
        #ax.w_yaxis.set_ticklabels([])
        #ax.w_zaxis.set_ticklabels([])
        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        ax.set_zlabel('PC 3')
        ax.set_xlim([1e10, 0.5e10])
        ax.set_ylim([5e3, 10e3])
        ax.set_zlim([4e3, -10e3])

    plt.savefig('/d/bip3/ezbc/scratch/pca.png')
    plt.show()

def cluster_data(data, n_clusters=2):

    from sklearn.cluster import KMeans
    from sklearn import datasets

    # Initialize the kmeans instance
    estimator = KMeans(n_clusters=n_clusters)

    # Fit the data
    estimator.fit(data)

    labels = estimator.labels_

    import matplotlib as mpl
    colors = labels.astype(np.float)
    colors /= colors.max()

    colors = mpl.cm.rainbow(colors)

    return colors




