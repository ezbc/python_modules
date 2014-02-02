#!/usr/bin/python

def bootstrap(data,model,num_samples=100,alpha=0.05,data_error=None,
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

