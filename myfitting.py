#!/usr/bin/python

def fit_gaussians(x, y, guesses=None, ncomps=1):

    """ Fits multiple Gaussians to y(x).

    Parameters
    ----------
    x, y : array-like
        Data.
    guesses : array-like, optional
        Initial guesses for fitting: [amp, x0, width] for each componenet. If
        ncomps > len(guesses) / 3  then the ncomps will be lowered to
        len(guesses) / 3.
    ncomps : int
        Number of Gaussians to fit

    Returns
    -------
    model : array-like
    model_comps : list
        Model components
    model_params : list
        list of model_params

    """

    import numpy as np
    from gaussfitter import multigaussfit

    def gauss(x, amp, x0, width):
        return amp * np.exp(-(x-x0)**2 / (2 * width**2))

    fit = multigaussfit(x, y, ngauss=int(ncomps), params=np.asarray(guesses))

    model = fit[1]
    model_params, model_comps = [], []

    for i in xrange(len(fit[0]) / 3):
        model_param = fit[0][3*i:3*i+3]
        model_params.append(model_param)
        model_comps.append(gauss(x, *model_param))

    return model, model_comps, model_params

