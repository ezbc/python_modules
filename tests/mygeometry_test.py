#!/usr/bin/python

import mygeometry as myg
import numpy as np

def test_point_in_polygon():

    point_test = (50, 50)

    polygon = np.array([[10, 10,],
                        [100, 10],
                        [100, 100],
                        [50, 150],
                        [10, 100]])

    print myg.point_in_polygon(point_test, polygon)
    assert myg.point_in_polygon(point_test, polygon)
