import numpy as np
import matplotlib.pyplot as plt

def scatter_contour(x, y,
                    levels=10,
                    fractional_levels=True,
                    threshold=100,
                    log_counts=False,
                    histogram2d_args={},
                    plot_args={},
                    errors=None,
                    contour_args={},
                    ax=None):

    """Scatter plot with contour over dense regions

    Parameters
    ----------
    x, y : arrays
        x and y data for the contour plot
    levels : integer or array (optional, default=10)
        number of contour levels, or array of contour levels
    fractional_levels : bool (optional, default=True)
        If True, then levels argument (if array-like) is interpreted as levels
        within which contain the fraction of data points specified. The
        threshold and log_counts arguments will not have an effect if True.
    threshold : float (default=100)
        number of points per 2D bin at which to begin drawing contours
    log_counts :boolean (optional)
        if True, contour levels are the base-10 logarithm of bin counts.
    histogram2d_args : dict
        keyword arguments passed to numpy.histogram2d
        see doc string of numpy.histogram2d for more information
    plot_args : dict
        keyword arguments passed to pylab.scatter
        see doc string of pylab.scatter for more information
    contourf_args : dict
        keyword arguments passed to pylab.contourf
        see doc string of pylab.contourf for more information
    ax : pylab.Axes instance
        the axes on which to plot.  If not specified, the current
        axes will be used

    Returns
    -------
    level_fractions : numpy.array
        Data fraction levels for each contour.

    """

    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()

    H, xbins, ybins = np.histogram2d(x, y, **histogram2d_args)

    Nx = len(xbins)
    Ny = len(ybins)

    # Get level fractions
    if isinstance(levels, int):
        if not log_counts:
            level_fractions = 1 - np.linspace(threshold, H.max(), levels) /\
                    (H.max() + threshold)
        else:
            H = np.log10(1 + H)
            threshold = np.log10(1 + threshold)
            level_fractions = 1 - np.logspace(np.log10(threshold + 1),
                                              np.log10(H.max() + 1),
                                              levels) / \
                              (H.max() + threshold + 2)

        levels = np.linspace(threshold, H.max(), levels)
    else:
        levels = np.asarray(levels)
        level_fractions = levels
        if fractional_levels:
            levels = np.sort((1.0 - levels) * np.max(H))
            levels = np.hstack((levels, np.max(H)))

    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]

    # draw a zero-width line: this gives us the outer polygon to
    # reduce the number of points we draw
    # somewhat hackish... we could probably get the same info from
    # the filled contour below.
    # make outline color same as edge of filled contour
    outline_colors = (plt.get_cmap()(0),)
    outline = ax.contour(H.T, levels=(np.min(levels),),
                         linewidths=0, extent=extent,
                         colors=outline_colors, alpha=1)

    outer_polys = outline.allsegs[0]

    #print outer_polys
    #print len(outer_polys)

    # Create filled contour
    C = ax.contourf(H.T, levels, extent=extent, **contour_args)

    # Make list of point coordinates, N x 2 in shape
    X = np.hstack([x[:, None], y[:, None]])

    points_inside = np.zeros(x.shape, dtype=bool)
    try:
        # this works in newer matplotlib versions
        from matplotlib.path import Path
        for outer_poly in outer_polys:
            points_inside[Path(outer_poly).contains_points(X)] = True
    except ImportError:
        # this works in older matplotlib versions
        import matplotlib.nxutils as nx
        points_inside = nx.points_inside_poly(X, outer_poly)

    Xplot = X[~points_inside]

    if errors is None:
        ax.plot(Xplot[:, 0], Xplot[:, 1], zorder=0, **plot_args)
    else:
        Xplot_error = errors[~points_inside]
        ax.errorbar(Xplot[:, 0], Xplot[:, 1], yerr=Xplot_error, zorder=0,
                    **plot_args)

    return level_fractions[:-1]

def plot_cdf_confint(data, data_error=0, ax=None, plot_kwargs_line={},
        plot_kwargs_fill_between={}, return_axis=False, nsim=100, nbins=20,
        bin_limits=None):

    ''' Performs Monte Carlo simulation with data error to calculate the CDF
    point-wise confidence interval.

    Parameters
    ----------
    data : array-like
        Distribution of data.
    data_error : float, array-like
        Normal error on data. If an array, must have same dimensions as data.
    ax : matplotlib.pyplot.axis, optional
        If provided, adds plot to axis object. Else plots matplotlib.pyplot.
    plot_kwargs_line : dict, optional
        Kwargs to provide to matplotlib.pyplot.plot for median CDF.
    plot_kwargs_fill_between : dict, optional
        Kwargs to provide to matplotlib.pyplot.fill_between for CDF confidence
        interval.
    return_axis : bool, optional
        Return the CDF bin values of the data?
    nsim : int, optional
        Number of Monte Carlo simulations to run.
    nbins : int, optional
        Number of bins with which to sample the simulated data.
    bin_limits : array-like, optional
        Lower and upper bound of bins with which to calculate point-wise
        confidence intervals.

    Returns
    -------
    x : array-like, optional
        If return_axis is True, then the CDF sample locations for the original
        dataset is returned.

    '''

    import mystats

    # Initialize CDF array
    cdfs = np.empty((nsim, data.size))
    xs = np.empty((nsim, data.size))

    # simulate different CDFs in monte carlo
    for i in xrange(nsim):
        data_sim = data + np.random.normal(scale=data_error)
        cdfs[i], xs[i] = mystats.calc_cdf(data_sim, return_axis=True)
        #ax.plot(xs[i], cdfs[i], alpha=0.05, color='k')

    # initialize new plotted x-values / bins for confidence interval
    if bin_limits is None:
        x_fit = np.linspace(np.min(xs), np.max(xs), nbins, endpoint=True)
    else:
        if bin_limits[0] is not None:
            lim_low = bin_limits[0]
        else:
            lim_low = np.min(xs)
        if bin_limits[1] is not None:
            lim_high = bin_limits[1]
        else:
            lim_high = np.max(xs)
        x_fit = np.linspace(lim_low, lim_high, nbins, endpoint=True)

    # initialize empty array for confidence interval
    cdf_confint = np.ones((3, x_fit.size))

    # Calculate the median and uncertainty on the median in each bin given the
    # simulation results
    for i in xrange(x_fit.size - 1):
        cdf_bin = cdfs[(xs >= x_fit[i]) & (xs < x_fit[i+1])]
        median, conf_err = mystats.calc_cdf_error(cdf_bin)
        cdf_confint[1, i] = median
        cdf_confint[0, i] = median - conf_err[0]
        cdf_confint[2, i] = median + conf_err[1]

    # Use center of bins to plot
    x_plot = x_fit + (x_fit[1] - x_fit[0]) / 2.0

    # eliminate nans
    nan_mask = (np.isnan(cdf_confint[0]) | \
                np.isnan(cdf_confint[1]) | \
                np.isnan(cdf_confint[2]))
    cdf_confint = cdf_confint[:, ~nan_mask]
    x_plot = x_plot[~nan_mask]

    # Plot the results with the median estimate and the confidence interval
    if ax is None:
        plt.plot(x_plot, cdf_confint[1], **plot_kwargs_line)
        plt.fill_between(x_plot, cdf_confint[0], cdf_confint[2],
                **plot_kwargs_fill_between)
    else:
        ax.plot(x_plot, cdf_confint[1], **plot_kwargs_line)
        ax.fill_between(x_plot, cdf_confint[0], cdf_confint[2],
                        **plot_kwargs_fill_between)

    # Return the original data x axis?
    if return_axis:
        return x

def set_color_cycle(num_colors=4, cmap=plt.cm.copper, cmap_limits=[0, 0.8]):

    # color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(cmap_limits[0],
                                                cmap_limits[1],
                                                num_colors)]

    params = {'axes.color_cycle': color_cycle, # colors of different plots
             }

    plt.rcParams.update(params)

    return color_cycle

def get_color_cycle(num_colors=4, cmap=plt.cm.copper, cmap_limits=[0, 0.8]):

    # color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(cmap_limits[0],
                                                cmap_limits[1],
                                                num_colors)]
    return color_cycle

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):

    ''' Truncates a matplolib.colors colormap to a smaller range.

    Parameters
    ----------
    cmap : matplotlib.pyplot.cm
        Colormap
    minval : float
        Lower value to truncate.
    maxval : float
        Upper value to truncate
    n : int
        Number of discrete samples of colormap between minval and maxval.

    Returns
    -------
    new_cmap : matplotlib.pyplot.cm
        Truncated colormap

    '''

    import matplotlib.colors as colors

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))

    return new_cmap

def corner_plot(distributions, plot_grids=None, labels=None,
        filename=None, logscale=False, confidence_intervals=None):

    '''

    Parameters
    ----------
    distributions : array-like
        N-dimensional array
    plot_grids : list
        List of arrays, where each array corresponds to the parameter values of
        the corresponding axis of the distributions.
    labels : list
        List of parameter strings.
    filename : string
        Filename to save to.
    logscale : bool
        Plot both 1D and 2D PDFs in logscale?
    confidence_intervals : array-like
        Confidence intervals for 1D distributions. N-dimesional array, with
        lower, best estimate and upper bounds.

    '''

    # Import external modules
    import numpy as np
    import matplotlib
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import myplotting as myplt

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    params = {
              'figure.figsize': (6, 6),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    nside = distributions.ndim

    fig, axes = plt.subplots(nside, nside)
    #fig = plt.figure()

    # create extents of each plot
    if plot_grids is None:
        plot_grids = []
        for i in xrange(distributions.ndim):
            plot_grids.append(np.arange(distributions.shape[i]))
    if labels is None:
        labels = [r'$\theta$' + str(i) for i in xrange(distributions.ndim)]

    # Cycle through the plots, labeling and plotting
    for x_i in xrange(nside):
        for y_j in xrange(nside):
            ax = axes[y_j, x_i]

            # determine which type of plot to use
            if x_i == y_j:
                if confidence_intervals is not None:
                    plot_confs = confidence_intervals[x_i]
                else:
                    plot_confs = None

                _plot_corner_1dhist(ax,
                                    distributions,
                                    x_i,
                                    plot_grids,
                                    logscale=logscale,
                                    plot_confs=plot_confs,
                                    )
                if not logscale:
                    # reduce number of tick
                    ax.locator_params(nbins = 4,)
                else:
                    ax.locator_params(nbins = 4, axis='x')
            elif x_i > y_j:
                fig.delaxes(ax)
            else:
                _plot_corner_2dhist(ax,
                                    distributions,
                                    (x_i, y_j),
                                    plot_grids,
                                    logscale=logscale,
                                    )
                # reduce number of tick
                ax.locator_params(nbins = 4)

            # turn labels on or off
            if y_j > 0 and x_i == 0:
                ax.set_ylabel(labels[y_j])
            else:
                ax.yaxis.set_ticklabels([])

            if y_j == nside - 1:
                ax.set_xlabel(labels[x_i])
            else:
                #no label
                ax.xaxis.set_ticklabels([])


    # Adjust asthetics
    #ax.set_xlabel(axis_label)
    #ax.set_ylabel('Counts')
    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')

def _plot_corner_1dhist(ax, distributions, plot_axis, plot_grids,
        logscale=False, plot_confs=None):

    # Get the axes to marginalize the distributions over
    dist_axes = np.arange(distributions.ndim)
    marg_axes = dist_axes[dist_axes != plot_axis]
    if not isinstance(marg_axes, np.int64):
        marg_axes = tuple(marg_axes)


    # Derive the histogram
    hist = np.nansum(distributions, axis=marg_axes)

    #hist = plot_grids[plot_axis]
    plot_grid = plot_grids[plot_axis]

    lines = ax.plot(plot_grid,
                    hist,
                    color='k',
                    drawstyle='steps-mid',
                    linewidth=2,
                    )

    if plot_confs is not None:
        ax.axvline(plot_confs[0],
                        color='k',
                        alpha=0.5,
                        linestyle='--',
                        linewidth=1,
                        )

        ax.axvline(plot_confs[1],
                   alpha=0.5,
                   linewidth=2,
                   color='k',
                   linestyle='--',
                   )

        ax.axvline(plot_confs[2],
                        color='k',
                        alpha=0.5,
                        linestyle='--',
                        linewidth=1,
                        )

    if logscale:
        ax.set_yscale('log')

    ax.set_ylim([0, 1.1 * np.max(hist)])
    ax.set_xlim([plot_grid[0], plot_grid[-1]])

def _plot_corner_2dhist(ax, distributions, plot_axes, plot_grids,
        logscale=False):

    from matplotlib.colors import LogNorm

    # Get the axes to marginalize the distributions over
    dist_axes = np.arange(distributions.ndim)

    marg_axes = \
            dist_axes[(dist_axes != plot_axes[0]) & (dist_axes!= plot_axes[1])]

    if not isinstance(marg_axes, np.int64):
        marg_axes = tuple(marg_axes)

    # Derive the histogram
    hist_2d = np.nansum(distributions, axis=marg_axes)

    # get extent
    xgrid = plot_grids[plot_axes[0]]
    ygrid = plot_grids[plot_axes[1]]

    if logscale:
        norm = LogNorm(vmin=np.min(hist_2d[hist_2d != 0]),
                       vmax=np.max(hist_2d))
    else:
        norm = None

    ax.imshow(hist_2d.T,
            extent=[xgrid[0], xgrid[-1], ygrid[0], ygrid[-1]],
            cmap=plt.cm.binary,
            origin='lower',
            aspect='auto',
            interpolation='nearest',
            norm=norm,
            )

def delete_overlapping_xlabels(fig, ax):

    ''' Deletes overlapping xtick labels in axis.

    Parameters
    ----------
    fig : matplotlib.pyplot.figure
        Figure object.

    ax : matplotlib.pyplot.axis
        Axis object.

    Returns
    -------
    ax : matplotlib.pyplot.axis
        Axis object without overlapping xtick labels.

    '''

    import numpy as np

    fig.canvas.draw()
    major_labels = ax.get_xticklabels(minor=False)
    new_major_text = [item.get_text() for item in major_labels]

    old_box = None
    for i, label in enumerate(major_labels):
        new_box = label.get_window_extent().get_points()[:, 0]

        if old_box is not None:
            if ((new_box[0] < old_box[1]) & \
                (new_box[0] > old_box[0])):
                new_major_text[i] = ''
                new_box = [np.nan, np.nan]

        old_box = new_box

    ax.set_xticklabels(new_major_text)

    # minor ticks
    fig.canvas.draw()
    minor_labels = ax.get_xticklabels(minor=True)
    new_minor_text = [item.get_text() for item in minor_labels]

    old_box = None
    for i, label in enumerate(minor_labels):
        new_box = label.get_window_extent().get_points()[:, 0]

        if old_box is not None:
            if (((new_box[0] > old_box[0]) & \
                 (new_box[0] < old_box[1])) | \
                ((new_box[1] > old_box[0]) & \
                 (new_box[1] < old_box[1]))):
                new_minor_text[i] = ''
                new_box = [np.nan, np.nan]

        old_box = new_box

    ax.set_xticklabels(new_minor_text, minor=True)

    # minor ticks against major
    fig.canvas.draw()
    minor_labels = ax.get_xticklabels(minor=True)
    major_labels = ax.get_xticklabels(minor=False)
    new_minor_text = [item.get_text() for item in minor_labels]

    for i, minor_label in enumerate(minor_labels):
        box_minor = minor_label.get_window_extent().get_points()[:, 0]

        for major_label in major_labels:
            box_major = \
                    major_label.get_window_extent().get_points()[:, 0]

            if (((box_minor[0] > box_major[0]) & \
                 (box_minor[0] < box_major[1])) | \
                ((box_minor[1] > box_major[0]) & \
                 (box_minor[1] < box_major[1]))):
                new_minor_text[i] = ''

    ax.set_xticklabels(new_minor_text, minor=True)

    return ax

def check_overlaps(bboxes):

    '''
    Returns a list of how many times each bbox overlaps with other bboxes.

    Paramters
    ---------
    bboxes : list
        List of matplotlib.transforms.Bbox objects.

    Returns
    -------
    overlaps : list
        List of number of overlaps the bbox has with neighbors.

    '''

    overlaps = [0] * len(bboxes)
    for i, box in enumerate(bboxes):
        for other_box in bboxes:
            if (box != other_box):
                #overlaps[i] += box.overlaps(other_box)
                overlaps[i] += check_overlaps(box, other_box)
    return overlaps

def check_overlaps(bbox1, bbox2):
    """
    Returns True if this bounding box overlaps with the given
    bounding box *other*. OR IF NANA!
    """
    ax1, ay1, ax2, ay2 = bbox1._get_extents()
    bx1, by1, bx2, by2 = bbox2._get_extents()

    if ax2 < ax1:
        ax2, ax1 = ax1, ax2
    if ay2 < ay1:
        ay2, ay1 = ay1, ay2
    if bx2 < bx1:
        bx2, bx1 = bx1, bx2
    if by2 < by1:
        by2, by1 = by1, by2

    if np.isnan(np.sum(extents1)) or np.isnan(np.sum(extents2)):
        return False

    return not ((bx2 < ax1) or
                (by2 < ay1) or
                (bx1 > ax2) or
                (by1 > ay2))

def reverse_colormap(cmap):

    ''' Reverses a matplolib.colors colormap.

    Parameters
    ----------
    cmap : matplotlib.pyplot.cm
        Colormap

    Returns
    -------
    cmap_r : matplotlib.pyplot.cm
        Reversed colormap

    '''

    from matplotlib.colors import ListedColormap

    cmap_r = ListedColormap(cmap.colors[::-1])

    return cmap_r

def convert_wcs_limits(limits, header, frame='fk5'):

    '''
    Parameters
    ----------
    limits: array-like
        (ra_min, ra_max, dec_min, dec_max) all in degrees

    header: astropy.io.fits.header

    Returns
    -------
        (ra_min, ra_max, dec_min, dec_max) all in pixels

    '''

    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from astropy.io import fits
    from astropy.wcs import WCS

    # write limits as SkyCoords instances
    coords_wcs_a = SkyCoord(limits[0], limits[2], unit='deg', frame=frame)
    coords_wcs_b = SkyCoord(limits[1], limits[3], unit='deg', frame=frame)

    # Convert limits from WCS to pixels
    wcs_header = WCS(header)
    coords_pix_a = coords_wcs_a.to_pixel(wcs_header)
    coords_pix_b = coords_wcs_b.to_pixel(wcs_header)

    limits = coords_pix_a[0], coords_pix_b[0], coords_pix_a[1], coords_pix_b[1]

    return limits

def plot_cdf(data, ax=None, plot_kwargs={}, return_axis=False):

    import mystats

    cdf, x = mystats.calc_cdf(data, return_axis=True)

    if ax is None:
        plt.plot(x, cdf, **plot_kwargs)
    else:
        ax.plot(x, cdf, **plot_kwargs)

    if return_axis:
        return x

def get_square_grid_sides(ngrids):

    ''' Calculates the number of rows and columns needed to encompass N
    subplots.

    '''

    n = int(np.ceil(ngrids**0.5))

    if n**2 - n >= ngrids:
        nrows = n - 1
        ncols = n
    else:
        nrows, ncols = n, n

    nrows_ncols = nrows, ncols

    return nrows_ncols



