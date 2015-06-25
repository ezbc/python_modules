#!/usr/bin/python

class cloud():

    # Import external modules
    import numpy as np

    def __init__(self, av_filename, hi_filename, av_error_filename=None,
            hi_error_filename=None, cloud_prop_filename=None):

        # Import external modules
        import pyfits as fits
        from mycoords import make_velocity_axis

        # Define local variables
        self.av_filename = av_filename
        self.hi_filename = hi_filename
        self.av_error_filename = av_error_filename
        self.hi_error_filename = hi_error_filename

        # Load data
        self.av_data, self.av_header = fits.getdata(av_filename, header=True)
        self.hi_data, self.hi_header = fits.getdata(hi_filename, header=True)
        if av_error_filename is not None:
            self.av_error_data, self.av_error_header = \
                    fits.getdata(av_error_filename, header=True)
        if hi_error_filename is not None:
            self.hi_error_data, self.hi_error_header = \
                fits.getdata(hi_error_filename, header=True)
        if cloud_prop_filename is not None:
            self.load_cloud_properties(cloud_prop_filename)

        # Make velocity axis for HI cube
        self.hi_vel_axis = make_velocity_axis(self.hi_header)

    def _derive_region_mask(self,):

        # Derive relevant region
        region_vertices = \
            self.props['regions'][self.region_name]['poly_verts']['pixel']

        # block off region
        region_mask = np.logical_not(myg.get_polygon_mask(self.av_data,
                                                          region_vertices))

        self.region_mask = region_mask

    def _convert_coordinates(self,
            coords=('region_limit','co_noise_limits','plot_limit'),
            header=self.av_header):

        ''' Converts WCS coordinates to pixel coordinates for a few sky
        positions.  '''

        # Initialize pixel keys
        for coord in coords:
            self.props[coord].update({'pixel': []})

            if coord in ('region_limit',
                         'plot_limit',
                         'region_limit_bin',
                         'plot_limit_bin'):
                limit_wcs = self.props[coord]['wcs']

                for limits in limit_wcs:
                    # convert centers to pixel coords
                    limit_pixels = get_pix_coords(ra=limits[0],
                                                 dec=limits[1],
                                                 header=header)[:2].tolist()

                    self.props[coord]['pixel'].append(limit_pixels[0])
                    self.props[coord]['pixel'].append(limit_pixels[1])
            elif coord == 'co_noise_limits':
                region_limits = self.props[coord]['wcs']

                # Cycle through each region, convert WCS limits to pixels
                for region in region_limits:
                    region_pixels = []
                    for limits in region:
                        # convert centers to pixel coords
                        limit_pixels = get_pix_coords(ra=limits[0],
                                                      dec=limits[1],
                                                      header=header)[:2].tolist()
                        region_pixels.append(limit_pixels)

                    # Append individual regions back to CO noise
                    self.props[coord]['pixel'].append(region_pixels)

    def load_cloud_properties(self, prop_filename):

        import json

        # Load global properties
        with open(prop_filename, 'r') as f:
            self.props = json.load(f)

    def load_region(self, region_filiename):

        import pyregion as pyr

        # region[0] in following format:
        # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
        # [ra center, dec center, width, height, rotation angle]

        regions = pyr.open(filename)

        self.props['regions'] = {}


        for region in regions:
            # Cores defined in following format: 'tag={L1495A}'
            tag = region.comment
            region_name = tag[tag.find('text={')+6:tag.find('}')].lower()

            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            poly_verts_pix = []
            for i in xrange(0, len(poly_verts)):
                poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                                dec=poly_verts[i][1],
                                                header=header)[:-1][::-1].tolist())

            self.props['regions'][region_name] = {}
            self.props['regions'][region_name]['poly_verts'] = {}
            self.props['regions'][region_name]['poly_verts']['wcs'] = poly_verts
            self.props['regions'][region_name]['poly_verts']['pixel'] = poly_verts_pix

    def run_analysis(self, region_filename=None, region='california'):

        self.region = region

        # Change WCS coords to pixel coords of images
        self._convert_coordinates()

        # Load cloud division regions from ds9
        self.load_region(region_filename)

        # derive the region mask
        self._derive_region_mask()



        # last left on line 1791 of
        # california_analysis_likelihood_iterative_mask_binning.py



