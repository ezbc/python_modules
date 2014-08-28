#!/usr/bin/python

def find_clumps(image, threshold):
    '''
    Finds regions with pixels above threshold.

    Returns
    -------
    label_im : array-like
        Array the shape of the image including region locations.

    '''

    from scipy import ndimage
    from sklearn.feature_extraction import image as sklearn_image
    from sklearn.cluster import spectral_clustering

    mask = image > threshold

    label_im, nb_labels = ndimage.label(mask)

    return label_im

def cut_clump_sizes(clump_image, size_threshold):

    ''' Cuts regions in clump image with fewer than the threshold number of
    pixels.

    Returns
    -------
    clean_clump_image : array-like
        Array the shape of the image including region locations.

    '''

    from scipy import ndimage
    from sklearn.feature_extraction import image as sklearn_image
    from sklearn.cluster import spectral_clustering

    # Select larger regions
    sizes = ndimage.sum(mask, clump_image, range(nb_labels + 1))
    mask_size = sizes < size_threshold
    remove_pixel = mask_size[clump_image]
    clump_image[remove_pixel] = 0
    labels = np.unique(clump_image)
    clean_clump_image = np.searchsorted(labels, clump_image)

    return clean_clump_image

def plot_clumps(image, header=None):

    # Simple plot
    if header is None:
        import matplotlib.pyplot as plt

        fig, ax0 = plt.subplots(ncols=1, figsize=(12, 12))

        ax0.imshow(image, cmap=plt.cm.spectral, origin='lower')

        ax0.axis('off')

        fig.show()
    elif header is not None: # Plot with RA and Dec
        # Import external modules
        import matplotlib.pyplot as plt
        import numpy as np
        from mpl_toolkits.axes_grid1 import ImageGrid
        import pyfits as pf
        import matplotlib.pyplot as plt
        import pywcsgrid2 as wcs
        import pywcs
        from pylab import cm # colormaps

        fig = plt.figure(figsize=(8,8))

        imagegrid = ImageGrid(fig, (1,1,1),
                     nrows_ncols=(1,1),
                     ngrids=1,
                     cbar_mode="single",
                     cbar_location='right',
                     cbar_pad="2%",
                     cbar_size='3%',
                     axes_pad=0,
                     axes_class=(wcs.Axes,
                                 dict(header=header)),
                     aspect=False,
                     label_mode='L',
                     share_all=False)

        ax = imagegrid[0]
        im = ax.imshow(image,
                interpolation='nearest',origin='lower',
                cmap=plt.cm.spectral)

        # Asthetics
        ax.set_display_coord_system("fk4")
        ax.set_ticklabel_type("hms", "dms")
        ax.set_xlabel('Right Ascension (J2000)',
                  size = 'small',
                  family='serif')
        ax.set_ylabel('Declination (J2000)',
                  size = 'small',
                  family='serif')
        # colorbar
        cb = ax.cax.colorbar(im)

        fig.show()

def choose_region(image, header=None, return_xy=False,
        limits=None):

    import matplotlib.pyplot as plt

    if header is None:
        fig, ax0 = plt.subplots(ncols=1, figsize=(12, 12))

        ax0.imshow(image, cmap=plt.cm.spectral, origin='lower')

        if limits is not None:
            ax0.set_xlim([limits[0], limits[1]])
            ax0.set_ylim([limits[2], limits[3]])
        ax0.axis('off')
        fig.show()

    elif header is not None: # Plot with RA and Dec
        # Import external modules
        import matplotlib.pyplot as plt
        import numpy as np
        from mpl_toolkits.axes_grid1 import ImageGrid
        import pyfits as pf
        import matplotlib.pyplot as plt
        import pywcsgrid2 as wcs
        import pywcs
        from pylab import cm # colormaps

        fig = plt.figure(figsize=(8,8))

        imagegrid = ImageGrid(fig, (1,1,1),
                     nrows_ncols=(1,1),
                     ngrids=1,
                     cbar_mode="single",
                     cbar_location='right',
                     cbar_pad="2%",
                     cbar_size='3%',
                     axes_pad=0,
                     axes_class=(wcs.Axes,
                                 dict(header=header)),
                     aspect=False,
                     label_mode='L',
                     share_all=False)

        ax = imagegrid[0]
        im = ax.imshow(image,
                interpolation='nearest',origin='lower',
                cmap=plt.cm.spectral)

        # Asthetics
        ax.set_display_coord_system("fk4")
        ax.set_ticklabel_type("hms", "dms")
        ax.set_xlabel('Right Ascension (J2000)',
                  size = 'small',
                  family='serif')
        ax.set_ylabel('Declination (J2000)',
                  size = 'small',
                  family='serif')
        # colorbar
        cb = ax.cax.colorbar(im)
        # limits
        if limits is not None:
            ax.set_xlim([limits[0], limits[1]])
            ax.set_ylim([limits[2], limits[3]])

        fig.show()

    image_regions = clump_image(image, fig)

    if not return_xy:
    	return image_regions.get_region_mask()

    elif return_xy:
        return image_regions.get_xy()

class clump_image():
    ''' Canvas image for choosing a region.
    '''

    x = None
    y = None
    clump_image = None
    cid = None
    fig = None
    region_locs = None

    def __init__(self, clump_image, fig):

        self.clump_image = clump_image
        self.fig = fig
        canvas = self.fig.canvas

        # add a callback that triggers when the text is clicked
        self.cid = canvas.mpl_connect('button_press_event',self.on_click)

        # start a blocking event loop
        canvas.start_event_loop(timeout=-1)

    def on_click(self, event):
        import matplotlib.pyplot as plt
        import numpy as np

        self.set_xy(event.ydata, event.xdata)

        if self.x is not None:
            # exit the blocking event loop
            self.fig.canvas.stop_event_loop()
            plt.close()

            # create mask with region_locs
            self.region_locs = np.where(self.clump_image == \
                    self.clump_image[self.x, self.y])

    def get_xy(self,):
        return self.x, self.y

    def get_image(self,):
        return self.clump_image

    def get_region_mask(self,):
        import numpy as np
        mask = np.ones(self.clump_image.shape)
        mask[self.region_locs] = 0
        return mask

    def set_xy(self, x, y):
        self.x, self.y = int(x), int(y)


