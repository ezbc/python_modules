#!/usr/bin/python


def get_rect(x, y, width, height, angle):

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
        Angle of rotation in degrees clockwise from East.
    '''

    import numpy as np

    # Create simple rectangle
    rect = np.array([(0, 0), (width, 0), (width, height), (0, height), (0, 0)])
    theta = (np.pi / 180.0) * angle

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
        Angle of rotation in degrees clockwise from East.

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

    import numpy as np
    import matplotlib.path as Path

    '''
    import numpy
    from PIL import Image, ImageDraw

    img = Image.new('L', (image.shape[1], image.shape[0]))
    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
    mask = numpy.array(img)
    '''

    from skimage.draw import polygon as get_poly_indices

    rr, cc = get_poly_indices(polygon[:, 0], polygon[:, 1], image.shape)

    if return_indices:
    	return (rr,cc)
    else:
        mask = np.zeros(image.shape)

        mask[rr,cc] = 1

        return mask

def point_in_polygon(target, poly):

    ''' Tests whether a point is in a polygon.

    Parameters
    ----------
    target : tuple
        (x,y) positions of point to test in polygon.
    poly : list
        List of tuples comprising the polygon.

    Returns
    -------
    inside : bool
        Whether or not the point is inside the polygon.

    '''

    from collections import namedtuple

    point = namedtuple("Point", ("x", "y"))
    line = namedtuple("Line", ("p1", "p2"))
    target = point(*target)

    inside = False
    # Build list of coordinate pairs
    # First, turn it into named tuples

    poly = map(lambda p: point(*p), poly)

    # Make two lists, with list2 shifted forward by one and wrapped around
    list1 = poly
    list2 = poly[1:] + [poly[0]]
    poly = map(line, list1, list2)

    for l in poly:
        p1 = l.p1
        p2 = l.p2

        if p1.y == p2.y:
            # This line is horizontal and thus not relevant.
            continue
        if max(p1.y, p2.y) < target.y <= min(p1.y, p2.y):
            # This line is too high or low
            continue
        if target.x < max(p1.x, p2.x):
            # Ignore this line because it's to the right of our point
            continue
        # Now, the line still might be to the right of our target point, but
        # still to the right of one of the line endpoints.
        rise = p1.y - p2.y
        run =  p1.x - p2.x
        try:
            slope = rise/float(run)
        except ZeroDivisionError:
            slope = float('inf')

        # Find the x-intercept, that is, the place where the line we are
        # testing equals the y value of our target point.

        # Pick one of the line points, and figure out what the run between it
        # and the target point is.
        run_to_intercept = target.x - p1.x
        x_intercept = p1.x + run_to_intercept / slope
        if target.x < x_intercept:
            # We almost crossed the line.
            continue

        inside = not inside

    return inside

