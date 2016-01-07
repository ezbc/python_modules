import numpy as np
import matplotlib.pyplot as plt


def delete_overlapping_xlabels(fig, ax):

    ''' Deletes overlapping xtick labels'''

    import numpy as np

    fig.canvas.draw()
    major_labels = ax.get_xticklabels(minor=False)
    new_major_text = [item.get_text() for item in major_labels]

    print new_major_text

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

    print new_minor_text

    old_box = None
    for i, label in enumerate(minor_labels):
        new_box = label.get_window_extent().get_points()[:, 0]

        print new_box, new_minor_text[i]

        if old_box is not None:
            if (((new_box[0] > old_box[0]) & \
                 (new_box[0] < old_box[1])) | \
                ((new_box[1] > old_box[0]) & \
                 (new_box[1] < old_box[1]))):
                print 'bounded'
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


def delete_overlapping_xlabels2(fig, ax):

    ''' Deletes overlapping xtick labels'''

    import numpy as np

    fig.canvas.draw()
    major_labels = ax.get_xticklabels(minor=False)
    new_major_text = [item.get_text() for item in major_labels]

    print new_major_text

    old_box = None
    for i in xrange(len(major_labels)):
        label = major_labels[i]
        new_box = label.get_window_extent().get_points()[:, 0]

        if old_box is not None:
            for j in xrange(i):
                old_box = major_labels[j].get_window_extent().get_points()[:, 0]
                if new_box[1] > old_box[0]:
                    new_major_text[i] = ''
                    new_box = [np.nan, np.nan]

        old_box = new_box

    print new_major_text
    ax.set_xticklabels(new_major_text)

    # minor ticks
    fig.canvas.draw()
    minor_labels = ax.get_xticklabels(minor=True)
    new_minor_text = [item.get_text() for item in minor_labels]

    print new_minor_text

    old_box = None
    for i in xrange(len(minor_labels)):
        label = minor_labels[i]
        new_box = label.get_window_extent().get_points()[:, 0]

        if old_box is not None:
            # check
            for j in xrange(i):
                old_box = minor_labels[j].get_window_extent().get_points()[:, 0]
                if new_box[1] > old_box[0]:
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


def delete_overlapping_xlabels3(fig, ax):

    ''' Deletes overlapping xtick labels'''

    import numpy as np

    fig.canvas.draw()
    major_labels = ax.get_xticklabels(minor=False)
    new_major_text = [item.get_text() for item in major_labels]
    bboxes = [label.get_window_extent() for label in major_labels]

    bbox_overlaps = check_overlaps(bboxes)
    overlaps = any(bbox_overlaps)

    while overlaps:
        i = np.argmax(bbox_overlaps)
        new_major_text[i] = ''
        bboxes[i].set_points(np.array([[np.nan, np.nan], [np.nan, np.nan]]))
        bbox_overlaps = check_overlaps(bboxes)
        overlaps = any(bbox_overlaps)

    ax.set_xticklabels(new_major_text)


def check_overlaps(bboxes):

    '''
    takes list of bboxes
    returns a list of how many times each bbox overlaps with other bboxes
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

def set_color_cycle(num_colors=4, cmap=plt.cm.copper, cmap_limits=[0, 0.8]):

    # color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(cmap_limits[0],
                                                cmap_limits[1],
                                                num_colors)]

    params = {'axes.color_cycle': color_cycle, # colors of different plots
             }

    plt.rcParams.update(params)

    return color_cycle

def corner_plot(distributions, plot_grids=None, labels=None,
        filename=None, logscale=False):

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
                _plot_corner_1dhist(ax,
                                    distributions,
                                    x_i,
                                    plot_grids,
                                    logscale=logscale,
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
        logscale=False, plot_confs=True):

    # Get the axes to marginalize the distributions over
    dist_axes = np.arange(distributions.ndim)
    marg_axes = dist_axes[dist_axes != plot_axis]
    if not isinstance(marg_axes, np.int64):
        marg_axes = tuple(marg_axes)

    # Derive the histogram
    hist = np.sum(distributions, axis=marg_axes)


    #hist = plot_grids[plot_axis]
    plot_grid = plot_grids[plot_axis]

    lines = ax.plot(plot_grid,
                    hist,
                    color='k',
                    drawstyle='steps-mid',
                    linewidth=2,
                    )

    if plot_confs:
        # get plotted hist
        #xdata = lines[0].get_path().vertices[:,0]
        #ydata = lines[0].get_path().vertices[:,1]
        xdata = plot_grid
        ydata = hist
        cdf = np.cumsum(ydata) / sum(ydata)

        # get the left and right boundary of the interval that contains 95% of
        # the probability mass
        right = np.argmax(cdf > 0.975)
        left = np.argmax(cdf > 0.025)
        center = np.argmax(cdf > 0.5)

        conf_interval = np.array([xdata[left],
                                  xdata[center],
                                  xdata[right],]
                                  )
        if 0:
            ax.fill_between(xdata[left:right],
                            ydata[left:right],
                            facecolor='k',
                            alpha=0.5)
        ax.axvline(xdata[center],
                   alpha=0.5,
                   linewidth=2,
                   color='k',
                   linestyle='--',
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
    hist_2d = np.sum(distributions, axis=marg_axes)

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


