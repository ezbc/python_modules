#!/usr/bin/python


def get_rect(x, y, width, height, angle,):

    ''' Returns four points of a rotated rectangle.

    Author: http://stackoverflow.com/questions/12638790/drawing-a-rectangle-inside-a-2d-numpy-array

    Parameters
    ----------
    x, y : int
        x and y positions of rectangle
    width : float
        Width of rectangle.
    height : float
        Height of rectangle.
    angle : float
        Angle of rotation in degrees counter clockwise from East.
    '''

    import numpy as np

    # Create simple rectangle
    rect = np.array([(-width/2., -height/2.),
                     (width/2., -height/2.),
                     (width/2., height/2.),
                     (-width/2., height/2.),
                     (-width/2., -height/2.)])
    theta = (np.pi / 180.0) * (angle )

    # Define four corners of rotated rectangle
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    offset = np.array([x, y])
    transformed_rect = np.dot(rect, R) + offset

    return transformed_rect[:-1]

def get_rectangular_mask(image, x, y, width=None, height=None, angle=0.0,
        return_indices=False):

    ''' Returns a mask excluding pixels outside of a rotated rectangular region.

    Parameters
    ----------
    image : array-like
        2D array to be masked.
    x, y : int
        x and y positions of rectangle
    width : float
        Width of rectangle.
    height : float
        Height of rectangle.
    angle : float
        Angle of rotation in degrees counter-clockwise from East.

    Returns
    -------
    masked_image : array-like
        Masked numpy array.
    '''

    # Import external modules
    import numpy as np

    # Define four corners of rectangle
    rectangle = get_rect(x, y, width, height, angle)

    if not return_indices:
        mask = get_polygon_mask(image, rectangle)

        masked_image = np.ma.array(image, mask=mask)
        return mask
    else:
        indices = get_polygon_mask(image, rectangle, return_indices=True)
        return indices

def get_polygon_mask(image, polygon, return_indices=False):

    '''
    Parameters
    ----------
    image : array-like
        N x M array.
    polygon : array-like
        K x 2 array of vertices.

    Returns
    -------
    mask : array-like
        N x M boolean array. True for any point inside of polygon.

    '''

    import numpy as np
    import numpy
    import matplotlib.path as Path
    from skimage.draw import polygon as get_poly_indices

    if type(polygon) is not numpy.array:
        polygon = np.copy(polygon)

    rr, cc = get_poly_indices(polygon[:, 0], polygon[:, 1], image.shape)

    if return_indices:
    	return (rr,cc)
    else:
        try:
            mask = np.zeros(image.shape)
            mask[rr,cc] = 1
            return mask
        except IndexError:
            raise IndexError('Polygon contains no pixels')

def point_in_polygon(xy, poly):

    ''' Tests whether a point is in a polygon.

    Parameters
    ----------
    xy : tuple
        (x,y) positions of point to test in polygon.
    poly : list
        List of tuples comprising the polygon.

    Returns
    -------
    inside : bool
        Whether or not the point is inside the polygon.

    '''

    x, y = xy

    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside

def rotate_polygon(polygon, anchor, angle):

    '''
    Parameters
    ----------
    polygon : array-like
        N x 2 array with x coordinates in column 0 and y coordinates in
        column 1.
    anchor : tuple
        x and y coordinates of pivot point.
    angle : float
        Angle to rotate polygon clockwise from North.

    Returns
    -------
    polygon_translate : array-like
        N x 2 array of rotated x and y coordinates.

    '''

    import numpy
    import numpy as np

    angle = np.deg2rad(angle)

    # If not arrays, make them
    if type(polygon) is not numpy.array:
        polygon = np.copy(polygon)
    if type(anchor) is not numpy.array:
        anchor = np.copy(anchor)

    # Center the polygon to the anchor point
    polygon_zero = polygon - anchor

    # Rotate the polygon
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle), np.cos(angle)]])

    polygon_rotated = np.dot(polygon_zero, rotation_matrix)

    # Translate polygon to original position
    polygon_translate = polygon_rotated + anchor

    return polygon_translate

def create_wedge(center_pos, radius, angle, center_rel_pos=0.1, width=None):

    '''
    Parameters
    ----------
    center_pos : array-like
        x and y pixel coordinates of core
    height : int, float
        Height of box along y axis.
    center_rel_pos : float, optional
        Core position in box along y-axis as a fraction of the height with the
        origin in the south.
    width : int, float
        Width of box along x axis.

    Returns
    -------
    wedge_vertices : numpy.array
        4 x 2 array with box pixel vertices, starting from lower-left going
        clockwise.

    '''

    from matplotlib.patches import Circle, Wedge, Polygon

    center_pos = (center_pos[0] - center_rel_pos * radius,
                  center_pos[1])

    wedge_vertices = Wedge(center_pos, radius, -angle/2., angle/2.,
            width=width).get_verts()

    return wedge_vertices

def rotate_wedge(wedge_vertices, anchor, angle):

    '''
    Parameters
    ----------
    wedge_vertices : numpy.array
        4 x 2 array with box pixel vertices, starting from lower-left going
        clockwise.
    anchor : tuple
        x and y coordinates of pivot point.
    angle : float
        Angle to rotate polygon clockwise from North.

    Returns
    -------
    wedge_vertices_rotated : numpy.array
        4 x 2 array with rotated box pixel vertices.
    '''

    from mygeometry import rotate_polygon

    wedge_vertices_rotated = rotate_polygon(wedge_vertices, anchor, angle)

    return wedge_vertices_rotated


