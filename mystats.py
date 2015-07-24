#!/usr/bin/python -W ignore::DeprecationWarning


import numpy as np
import matplotlib.pyplot as plt
import random
import warnings
warnings.filterwarnings('ignore')

class rv2d_discrete(object):

    ''' A generic discrete 2D random variable class meant for subclassing.
    Similar to scipy.stats.rv_discrete.

    Parameters
    ----------
    likelihoods : array-like
        N x M array of relative likelihoods corresponding to parameter 1 and 2
        values.
    param_grid1, param_grid2 : array-like
        Parameter values corresponding to element positions of likelihoods.
        The lengths of param_grid1 and param_grid2 must be N and M respectively.
    param_name1, param_name2 : str
        Names of parameters 1 and 2.
    L_scalar : float
        Inverse scalar of the likelihoods to calculate the 2D PDF. The
        likelihoods are divided by the minimum non-zero likelihood otherwise.

    Examples
    --------


    '''

    import numpy as np

    def __init__(self, likelihoods=None, param_grid1=None, param_grid2=None,
            param_name1='param1', param_name2='param2', L_scalar=None):

        super(rv2d_discrete, self).__init__()
        self.likelihoods = np.squeeze(likelihoods)
        self.likelihoods[self.likelihoods < 1e-16] = 0.0
        self.pdf = None
        self.param_grid1 = param_grid1
        self.param_grid2 = param_grid2
        self.param_name1 = param_name1
        self.param_name2 = param_name2
        if L_scalar is None:
            self.L_scalar = int(1.0 / np.min(likelihoods[likelihoods > 1e-8]))
        else:
            self.L_scalar = L_scalar

        # Scale likelihoods so that min value is an integer
        likelihoods_scaled = np.floor(self.likelihoods * self.L_scalar)

        # Initialize
        self.pdf = np.empty((np.sum(likelihoods_scaled), 2))
        count = 0

        # Create numbers of parameter pairs proportional to the parameter pair
        # likelihood
        for i, param1 in enumerate(param_grid1):
            for j, param2 in enumerate(param_grid2):
                L = likelihoods_scaled[i, j]
                if L > 0:
                    self.pdf[count:count + L] = (param1, param2)
                    count += L

    def rvs(self,):

        ''' Returns parameters random sample from the pdf.
        '''

        from numpy.random import randint

        # Get a random index of the pdf
        index = randint(0, len(self.pdf[:, 0]))

        # extract the parameters from the pdf
        params = self.pdf[index]

        return params

class rv3d_discrete(object):

    ''' A generic discrete 3D random variable class meant for subclassing.
    Similar to scipy.stats.rv_discrete.

    Parameters
    ----------
    likelihoods : array-like
        N x M X P array of relative likelihoods corresponding to parameter 1
        2, and 3 values.
    param_grid1, param_grid2, param_grid3: array-like
        Parameter values corresponding to element positions of likelihoods.
        The lengths of param_grid1, param_grid2 and param_grid3 must be N and M
        and P respectively.
    param_name1, param_name2, param_name3 : str
        Names of parameters 1, 2 and 3.
    L_scalar : float
        Inverse scalar of the likelihoods to calculate the 3D PDF. The
        likelihoods are divided by the minimum non-zero likelihood otherwise.

    Examples
    --------


    '''

    import numpy as np

    def __init__(self, likelihoods=None, param_grid1=None, param_grid2=None,
            param_grid3=None, param_name1='param1', param_name2='param2',
            param_name3='param3', L_scalar=None):

        super(rv3d_discrete, self).__init__()
        #self.likelihoods = np.squeeze(likelihoods)
        self.likelihoods = likelihoods
        self.likelihoods[self.likelihoods < 1e-16] = 0.0
        self.pdf = None
        self.param_grid1 = param_grid1
        self.param_grid2 = param_grid2
        self.param_grid3 = param_grid3
        self.param_name1 = param_name1
        self.param_name2 = param_name2
        self.param_name3 = param_name3
        if L_scalar is None:
            self.L_scalar = int(1.0 / np.min(likelihoods[likelihoods > 1e-8]))
        else:
            self.L_scalar = L_scalar

        # Scale likelihoods so that min value is an integer
        likelihoods_scaled = np.floor(self.likelihoods * self.L_scalar)

        # Initialize
        self.pdf = np.empty((np.sum(likelihoods_scaled), 3))
        count = 0

        # Create numbers of parameter pairs proportional to the parameter pair
        # likelihood
        for i, param1 in enumerate(param_grid1):
            for j, param2 in enumerate(param_grid2):
                for k, param3 in enumerate(param_grid3):
                    L = likelihoods_scaled[i, j, k]
                    if L > 0:
                        self.pdf[count:count + L] = (param1, param2, param3)
                        count += L

    def rvs(self,):

        ''' Returns parameters random sample from the pdf.
        '''

        from numpy.random import randint

        # Get a random index of the pdf
        index = randint(0, len(self.pdf[:, 0]))

        # extract the parameters from the pdf
        params = self.pdf[index]

        return params

def calc_symmetric_error(x, y=None, alpha=0.05):

    '''

    Parameters
    ----------
    x : array-like

    y : array-like, optional
        If provided, treated as the PDF of x

    '''


    import numpy as np

    from scipy.integrate import simps as integrate


    if len(x) < 4:
        #raise ValueError('x and y must have lengths > 3')
        return x[0], 0, 0

    # Create histogram with bin widths normalized by the density of values
    if y is None:
        x = np.sort(x)
        y = np.ones(x.shape)

    if np.any(y < 0):
        raise ValueError('y values mush be greater than 0')

    confidence = (1.0 - alpha)

    # area under whole function
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",category=DeprecationWarning)
        area = integrate(y, x)

    # Get weighted average of PDF
    mid_pos = np.argmin(np.abs(x - np.average(x, weights=y)))
    #mid_pos = np.argmin(np.abs(x - np.median(x, weights=y)))
    #mid_pos = np.interp

    # If the cum sum had duplicates, then multiple median pos will be derived,
    # take the one in the middle.
    try:
        if len(mid_pos) > 1:
            mid_pos = mid_pos[len(mid_pos) / 2]
    except TypeError:
        pass

    # Lower error
    pos = mid_pos - 1
    low_area = -np.Inf

    #area = integrate(y[0:mid_pos], x[0:mid_pos])

    while low_area <= area * confidence / 2.0 and pos > 0:
        y_clip = y[pos:mid_pos]# + 1]
        x_clip = x[pos:mid_pos]# + 1]

        low_area = integrate(y_clip, x_clip)

        # Catch the error if going to far
        if pos < 0:
            pos = 0
            break

        pos -= 1

    # set result to lower position
    low_pos = pos

    if pos == 0:
        low_pos = np.min(np.where(y != 0))

    # higher error
    pos = mid_pos + 1
    max_pos = len(x)
    high_area = -np.Inf

    #area = integrate(y[mid_pos:-1], x[mid_pos:-1])
    while high_area <= area * confidence / 2.0 and pos < max_pos:
        y_clip = y[mid_pos:pos]
        x_clip = x[mid_pos:pos]

        high_area = integrate(y_clip, x_clip)

        if pos > max_pos:
            pos = max_pos
            break

        pos += 1

    high_pos = pos

    if pos >= max_pos:
        high_pos = np.max(np.where(y != 0))

    median = x[mid_pos]
    low_error = x[mid_pos] - x[low_pos]
    high_error = x[high_pos] - x[mid_pos]

    return median, high_error, low_error

def calc_cdf_error(x, alpha=0.05):

    import numpy as np

    from scipy.integrate import simps as integrate

    if len(x) != len(y):
        raise ValueError('x and y must be the same shape')
    if len(x) < 4:
        raise ValueError('x and y must have lengths > 3')

    cdf = np.cumsum(y) / np.sum(y)

    import matplotlib.pyplot as plt
    plt.plot(x, cdf)
    plt.show()

    mid_pos = np.argmin(np.abs(cdf - 0.5)) + 1
    low_pos = np.argmin(np.abs(cdf - alpha / 2.0))
    high_pos = np.argmin(np.abs(alpha / 2.0 - cdf))

    median = x[mid_pos]
    low_error = x[mid_pos] - x[low_pos]
    high_error = x[high_pos] - x[mid_pos]

    return median, high_error, low_error

# Bootstrapping using medians
def bootstrap(data, num_samples):

    ''' Bootstraps data to determine errors. Resamples the data num_samples
    times. Returns errors of a bootstrap simulation at the 100.*(1 - alpha)
    confidence interval.

    Parameters
    ----------
    data : array-like
        Array of data in the form of an numpy.ndarray
    num_samples : int
        Number of times to resample the data.

    Returns
    -------
    conf_int : tuple, float
        Lower error and upper error at 100*(1-alpha) confidence of the data.
    samples : array-like
        Array of each resampled data. Will have one extra dimension than the
        data of length num_samples, representing each simulation.

    Notes
    -----
    -> arrays can be initialized with numpy.empty
    -> random samples can be retrieved from an array with random.sample

    Examples
    --------
    >>> import scipy
    >>> import numpy as np
    >>> data = scipy.random.f(1, 2, 100)
    >>> data.shape
    (100,)
    >>> samples = bootstrap(data, 50)
    >>> samples.shape
    (50, 100,)
    '''

    samples = np.empty((num_samples, data.size))

    for i in range(num_samples):
        indices = np.random.randint(0, data.size, data.size)
        samples[i,:] = data[indices]


    return samples

def calc_bootstrap_error(samples, alpha):

    ''' Returns errors of a bootstrap simulation at the 100.*(1 - alpha)
    confidence interval. Errors are computed by deriving a cumulative
    distribution function of the medians of the sampled data and determining the
    distance between the median and the value including alpha/2 % of the data,
    and the value including alpha/2 % of the data.

    Parameters
    ----------
    samples : array-like
        Array of each resampled data.

    Returns
    -------
    conf_int : tuple, float
        Median of the data, the lower error and the upper error at 100*(1-alpha)
        confidence of the data.

    Notes
    -----
    -> To find the index in an array closest to a given value, use the
    numpy.argmin function to find the index of the minimum value in an array.
    For example to find the value closest to 11.1 in an array of 10, 11, and 12:
        >>> import numpy as np
        >>> a = np.array([10, 11, 12])
        >>> print(np.argmin(np.abs(a - 11.1)))
        1

    Examples
    --------
    >>> import scipy
    >>> import numpy as np
    >>> data = scipy.random.f(1, 2, 100)
    >>> samples = bootstrap(data, 50)
    >>> errors = calc_bootstrap_error(samples, 0.05)

    '''

    medians, cdf = calc_cdf(samples)
    median = medians[np.argmin(np.abs(cdf - 0.5))]
    error_low = medians[np.argmin(np.abs(cdf - alpha/2.))]
    error_high = medians[np.argmin(np.abs(cdf - (1 - alpha/2.)))]

    return (median, median - error_low, error_high - median)

def calc_cdf(samples):

    ''' Calculates a cumulative distribution function of the medians of each
    instance of resampled data.

    Parameters
    ----------
    samples : array-like
        Array of each resampled data.

    Returns
    -------
    medians : array-like
        Array containing mean values for the cdf.
    cdf : array-like
        Array containing fraction of data below value x.

    '''

    medians = np.sort(np.median(samples, axis=0))
    cdf = np.cumsum(medians) / np.sum(medians)

    return medians, cdf

# Bootstrapping using means
def bootstrap(data, num_samples):

    ''' Bootstraps data to determine errors. Resamples the data num_samples
    times. Returns errors of a bootstrap simulation at the 100.*(1 - alpha)
    confidence interval.

    Parameters
    ----------
    data : array-like
        Array of data in the form of an numpy.ndarray
    num_samples : int
        Number of times to resample the data.

    Returns
    -------
    conf_int : tuple, float
        Lower error and upper error at 100*(1-alpha) confidence of the data.
    samples : array-like
        Array of each resampled data. Will have one extra dimension than the
        data of length num_samples, representing each simulation.

    Notes
    -----
    -> arrays can be initialized with numpy.empty
    -> random samples can be retrieved from an array with random.sample

    Examples
    --------
    >>> import scipy
    >>> import numpy as np
    >>> data = scipy.random.f(1, 2, 100)
    >>> data.shape
    (100,)
    >>> samples = bootstrap(data, 50)
    >>> samples.shape
    (50, 100,)
    '''

    samples = np.empty((num_samples, data.size))

    for i in range(num_samples):
        indices = np.random.randint(0, data.size, data.size)
        samples[i,:] = data[indices]

    return samples

def calc_bootstrap_error(samples, alpha):

    ''' Returns errors of a bootstrap simulation at the 100.*(1 - alpha)
    confidence interval. Errors are computed by deriving a cumulative
    distribution function of the means of the sampled data and determining the
    distance between the mean and the value including alpha/2 % of the data,
    and the value including alpha/2 % of the data.

    Parameters
    ----------
    samples : array-like
        Array of each resampled data.

    Returns
    -------
    conf_int : tuple, float
        Mean of the data, the lower error and the upper error at 100*(1-alpha)
        confidence of the data.

    Notes
    -----
    -> To find the index in an array closest to a given value, use the
    numpy.argmin function to find the index of the minimum value in an array.
    For example to find the value closest to 11.1 in an array of 10, 11, and 12:
        >>> import numpy as np
        >>> a = np.array([10, 11, 12])
        >>> print(np.argmin(np.abs(a - 11.1)))
        1

    Examples
    --------
    >>> import scipy
    >>> import numpy as np
    >>> data = scipy.random.f(1, 2, 100)
    >>> samples = bootstrap(data, 50)
    >>> errors = calc_bootstrap_error(samples, 0.05)

    '''

    means, cdf = calc_cdf(samples)

    import matplotlib.pyplot as plt
    plt.plot(means, cdf)
    plt.show()

    mean = means[np.argmin(np.abs(cdf - 0.5))]
    error_low = means[np.argmin(np.abs(cdf - alpha/2.))]
    error_high = means[np.argmin(np.abs(cdf - (1 - alpha/2.)))]

    return (mean, mean - error_low, error_high - mean)

def calc_cdf(samples):

    ''' Calculates a cumulative distribution function of the means of each
    instance of resampled data.

    Parameters
    ----------
    samples : array-like
        Array of each resampled data.

    Returns
    -------
    means : array-like
        Array containing mean values for the cdf.
    cdf : array-like
        Array containing fraction of data below value x.

    '''

    means = np.sort(np.mean(samples, axis=1))
    cdf = np.cumsum(means) / np.sum(means)

    return means, cdf

def calc_pdf(x, y):

    '''
    Calculates probability density function of the data. Uses a non-parametric
    approach to estimate the PDF.

    '''

    from scipy import interpolate

    inverse_density_function = interpolate.interp1d(x, y)

    return inverse_density_function

def bootstrap_model(data,model,num_samples=100,alpha=0.05,data_error=None,
        sigma=None, verbose=True):

    ''' Bootstraps data with models a given number of times and calculates the
    Goodness of fit for each run. The standard deviation of the Goodness-of-fit
    values is then used to estimate the confidence interval.

    Parameters
    ----------
    data : array_like
        The observed data, must be the same size as the model
    model : array_like
        The model data, must be the same size as the observed data.
    num_samples : int, optional
        Number of runs in the bootstrapping.
    alpha : float, optional
        Significance of confidence interval.
    data_error : float, array_like, optional
        If unset, the error will be the standard deviation of the data. If an
        array, it must have the same dimensions as the observed data.
    sigma : float, optional
        If set, the confidence interval will be calculated using the number of
        standard deviations from the mean.
    verbose : bool, optional
        Print out progress?

    Returns
    -------
    out : list
        A list, [confidence interval, goodness of fit array]

    '''

    import numpy as np
    from scipy.stats import norm

    data_list = data.ravel()
    model_list = model.ravel()
    length = len(data_list)
    #indices = np.arange(0,length,1)

    if data_error is None:
        data_error_list = data.std()
    else:
        data_error_list = data_error.ravel()

    num_samples = int(num_samples)
    gofArray = np.zeros(num_samples)

    if verbose:
        print('Beginning bootstrapping')

    for i in range(num_samples):
        # randomly sample all values of data and model
        indices_sample = np.random.choice(length,size=length,replace=True)
        data_sample = data_list[indices_sample]
        model_sample = model_list[indices_sample]
        gofArray[i] = ((data_sample - model_sample)**2 / \
                data_error_list**2).sum()
        if verbose:
            if i%10 == 0:
                print(str(i) + 'th run complete.')

    mean, std = gofArray.mean(), gofArray.std()
    if sigma is not None:
        alpha = 1 - norm.cdf(sigma)
    confid_int = norm.interval(1 - alpha, loc=mean, sigma=std)

    return (confid_int,gofArray)

def get_rms(x, axis=None):
    ''' Calculates the rms of an array.
    '''

    return np.sqrt(np.mean(x**2, axis=axis))

def fvalue(chi1,chi2,dof1,dof2):
    return (chi1/float(dof1))/(chi2/float(dof2))

def ftest(chi1,chi2,dof1,dof2):
    ''' The function ftest() c omputes the probability for a value drawn
    from the F-distribution to equal or exceed the given value of F.
    This can be used for confidence testing of a measured value obeying
    the F-distribution (i.e., ffor testing the ratio of variances, or
    equivalently for the addition of parameters to a fitted model).

      P_F(X > F; DOF1, DOF2) = PROB

    In specifying the returned probability level the user has three
    choices:

      * return the confidence level when the /CLEVEL keyword is passed;
        OR

      * return the significance level (i.e., 1 - confidence level) when
        the /SLEVEL keyword is passed (default); OR

      * return the "sigma" of the probability (i.e., compute the
        probability based on the normal  distribution) when the /SIGMA
        keyword is passed.

    Note that /SLEVEL, /CLEVEL and /SIGMA are mutually exclusive.

    For the ratio of variance test, the two variances, VAR1 and VAR2,
    should be distributed according to the chi-squared distribution
    with degrees of freedom DOF1 and DOF2 respectively.  The F-value is
    computed as:

       F = (VAR1/DOF1)  / (VAR2/DOF2)

    and then the probability is computed as:

       PROB = MPFTEST(F, DOF1, DOF2, ... )


    For the test of additional parameters in least squares fitting, the
    user should perform two separate fits, and have two chi-squared
    values.  One fit should be the "original" fit with no additional
    parameters, and one fit should be the "new" fit with M additional
    parameters.

      CHI1 - chi-squared value for original fit

      DOF1 - number of degrees of freedom of CHI1 (number of data
             points minus number of original parameters)

      CHI2 - chi-squared value for new fit

      DOF2 - number of degrees of freedom of CHI2

    Note that according to this formalism, the number of degrees of
    freedom in the "new" fit, DOF2, should be less than the number of
    degrees of freedom in the "original" fit, DOF1 (DOF2 < DOF1); and
    also CHI2 < CHI1.

    With the above definition, the F value is computed as:

      F = ( (CHI1-CHI2)/(DOF1-DOF2) )   /  (CHI2/DOF2)

    where DOF1-DOF2 is equal to M, and then the F-test probability is
    computed as:

       PROB = MPFTEST(F, DOF1-DOF2, DOF2, ... )

    Note that this formalism assumes that the addition of the M
    parameters is a small peturbation to the overall fit.  If the
    additional parameters dramatically changes the character of the
    model, then the first "ratio of variance" test is more appropriate,
    where F = (CHI1/DOF1) / (CHI2/DOF2).
    '''

    from scipy.stats import f

    return 1 - f.cdf( (chi1/float(dof1)) / (chi2/float(dof2)), dof1,dof2)

def test_bootstrap():
    import numpy as np
    from scikits.bootstrap import ci

    data = np.random.normal(loc=1, scale=1, size=1000)
    print('std = %.2f' % data.std())

    samples = bootstrap(data, 100)
    boot_error = calc_bootstrap_error(samples, 0.32)

    boot_error_ci = ci(data, np.median, 0.32)

    print('bootstrap error', boot_error)
    print('bootstrap error ci', boot_error_ci)

def main():
    test_bootstrap()

# Likelihood calculations

def calc_logL(model, data, data_error=None, weights=None):

    '''
    Calculates log likelihood

    http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node2.html

    '''

    import numpy as np

    if data_error is None:
        data_error = np.std(data)
    if isinstance(data_error, int):
        data_error = data_error * np.ones(data.shape)

    if weights is None:
        weights = 1.0
        data_weighted = data
        data_error_weighted = data_error
        model_weighted = model
    else:
        weights = weights[weights > 0]
        weights = np.array(weights / np.nanmin(weights), dtype=int)
        data_weighted = np.zeros(np.sum(weights))
        data_error_weighted = np.zeros(np.sum(weights))
        model_weighted = np.zeros(np.sum(weights))
        count = 0
        for i in xrange(0, len(weights)):
            data_weighted[count:count + weights[i]] = data[i]
            data_error_weighted[count:count + weights[i]] = data_error[i]
            model_weighted[count:count + weights[i]] = model[i]
            count += weights[i]

    #logL = -np.nansum((data - model)**2 / (2 * (data_error)**2))
    logL = -np.nansum((data_weighted - model_weighted)**2 / \
            (2 * (data_error_weighted)**2))

    return logL

def gauss(x, width, amp, x0):
    import numpy as np

    return amp * np.exp(-(x - x0)**2 / (2 * width**2))

if __name__ == '__main__':
    main()

