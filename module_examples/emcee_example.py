#!/usr/bin/python

'''Example usage of emcee module.

http://dan.iel.fm/emcee/current/user/line/

'''


import numpy as np

# Choose the "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
y = m_true*x+b_true
y += np.abs(f_true*y) * np.random.randn(N)
y += yerr * np.random.randn(N)

# Likelihood function
def lnlike(theta, x, y, yerr):
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))

    logL = -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
    return logL

# Log likelihood priors
def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf

# log likelihood function
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

# initialize the walkers in a tiny Gaussian ball around the maximum likelihood
# result
ndim, nwalkers = 3, 10
pos = [np.array((-1, 4.2, 1)) + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

pos = [(0, 0, 0) + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# Set up the sampler
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr),
        threads=1)

# Run the sampler
sampler.run_mcmc(pos, 100)

# Flatten the chain
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

# plot
import triangle
fig = triangle.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],
                      truths=[m_true, b_true, np.log(f_true)])
fig.savefig("triangle.png")

import matplotlib.pyplot as pl
pl.clf()
xl = np.array([0, 10])
for m, b, lnf in samples[np.random.randint(len(samples), size=100)]:
    pl.plot(xl, m*xl+b, color="k", alpha=0.1)
pl.plot(xl, m_true*xl+b_true, color="r", lw=2, alpha=0.8)
pl.errorbar(x, y, yerr=yerr, fmt=".k")
#pl.show()

samples[:, 2] = np.exp(samples[:, 2])
print samples.shape
np.save('mcmc_samples.npy', samples)
m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))

print m_mcmc, b_mcmc, f_mcmc



