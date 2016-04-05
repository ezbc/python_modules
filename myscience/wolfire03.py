#!/usr/bin/python

'''
Module for using the model of Wolfire et al. (2003)
'''

def calc_T(n, Z=1.0, n_error=(0.0,0.0), Z_error=(1.0,1.0), G_0=1.0,
        G_0_error=(0.0,0.0), use_fsolve=False,):

    ''' Calculates T from equation 18 in Krumholz et al. (2009) and equation C16
    in Wolfire et al. (2003)

    Parameters
    ----------
    Z : float, array-like
        Metallicity normalized by solar value.
    G0 : float, array-like
        FUV radiation field normalized by solar value.

    Returns
    -------
    T : float, array-like
        Temperature of gas in K.

    '''

    import numpy as np
    from scipy.optimize import fsolve

    # change tuples or lists to arrays
    n = np.array(n)
    Z =np.array(Z)
    G_0 =np.array(G_0)
    n_error =np.array(n_error)
    Z_error =np.array(Z_error)
    G_0_error =np.array(G_0_error)


    if use_fsolve:
        # write n as a function of T
        calc_resid = lambda T2: n - 20.0 * G_0 * T2**-0.2 * np.exp(1.5 / T2) / \
                                (1 + 2.6 * (T2**0.5 * Z)**0.365)

        # numerically solve for T
        initial_guess = 187.0 / 100.0 * (n / 10.0)
        T2 = fsolve(calc_resid, initial_guess, maxfev=10000, xtol=1e-16)

        T = T2 * 100.0
    else:

        T2 = np.logspace(0, 6, 10000) / 100.0

        # write n as a function of T
        n_interp = 20.0 * G_0 * T2**-0.2 * np.exp(1.5 / T2) / \
                   (1 + 2.6 * (T2**0.5 * Z)**0.365)

        T2_interp = np.interp(n,
                           n_interp,
                           T2,
                           )

        if np.size(n) > 1:
            T2_solution = np.empty(np.shape(n))
            for element in n:
                indices = np.where(n == element)[0]
                T2_solution[indices] = T2[np.argmin(np.abs(n_interp - element))]
        else:
            T2_solution = T2[np.argmin(np.abs(n_interp - n))]

        T = T2_solution * 100.0

    return T


