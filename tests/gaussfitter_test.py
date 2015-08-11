#!/usr/bin/python

def test_fit_gaussians():

    import numpy as np
    import matplotlib.pyplot as plt
    from myfitting import fit_gaussians

    def gauss(x, amp, x0, width):
        return amp * np.exp(-(x-x0)**2 / (2 * width**2))

    x = np.linspace(-20,15,1000)
    y = np.exp(-(x-3)**2 / (2 * 1**2)) + \
        3 * np.exp(-(x-9)**2 / (2 * 0.5**2)) + \
        5 * np.exp(-(x+0.5)**2 / (2 * 5**2))

    # amp, offset, width
    params = (1, 3, 1,
              2, 9, 1,
              0.5, -1, 1)

    fit = fit_gaussians(x, y, ncomps=3, guesses=params)

    assert len(fit) == 3
    assert len(fit[1]) == 3

