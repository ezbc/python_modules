#!/usr/bin/python

import numpy
from numpy.ma import median
from numpy import pi
#from scipy import optimize,stats,pi
from agpy.mpfit import mpfit
from gaussfitter import n_gaussian

def multigaussfit(xax, data, ngauss=1, err=None, params=[1,0,1],
        fixed=[False,False,False], limitedmin=[False,False,True],
        limitedmax=[False,False,False], minpars=[0,0,0], maxpars=[0,0,0],
        quiet=True, shh=True, veryverbose=False):
    """
    An improvement on onedgaussfit.  Lets you fit multiple gaussians.

    Inputs:
       xax - x axis
       data - y axis
       ngauss - How many gaussians to fit?  Default 1 (this could supersede onedgaussfit)
       err - error corresponding to data

     These parameters need to have length = 3*ngauss.  If ngauss > 1 and length = 3, they will
     be replicated ngauss times, otherwise they will be reset to defaults:
       params - Fit parameters: [amplitude, offset, width] * ngauss
              If len(params) % 3 == 0, ngauss will be set to len(params) / 3
       fixed - Is parameter fixed?
       limitedmin/minpars - set lower limits on each parameter (default: width>0)
       limitedmax/maxpars - set upper limits on each parameter

       quiet - should MPFIT output each iteration?
       shh - output final parameters?

    Returns:
       Fit parameters
       Model
       Fit errors
       chi2
    """

    if len(params) != ngauss and (len(params) / 3) > ngauss:
        ngauss = len(params) / 3 

    if isinstance(params,numpy.ndarray): params=params.tolist()

    # make sure all various things are the right length; if they're not, fix them using the defaults
    for parlist in (params,fixed,limitedmin,limitedmax,minpars,maxpars):
        if len(parlist) != 3*ngauss:
            # if you leave the defaults, or enter something that can be multiplied by 3 to get to the
            # right number of gaussians, it will just replicate
            if len(parlist) == 3: 
                parlist *= ngauss 
            elif parlist==params:
                parlist[:] = [1,0,1] * ngauss
            elif parlist==fixed or parlist==limitedmax:
                parlist[:] = [False,False,False] * ngauss
            elif parlist==limitedmin:
                parlist[:] = [False,False,True] * ngauss
            elif parlist==minpars or parlist==maxpars:
                parlist[:] = [0,0,0] * ngauss

    def mpfitfun(x,y,err):
        if err is None:
            def f(p,fjac=None): return [0,(y-n_gaussian(pars=p)(x))]
        else:
            def f(p,fjac=None): return [0,(y-n_gaussian(pars=p)(x))/err]
        return f

    if xax == None:
        xax = numpy.arange(len(data))

    parnames = {0:"AMPLITUDE",1:"SHIFT",2:"WIDTH"}

    parinfo = [ {'n':ii, 'value':params[ii],
        'limits':[minpars[ii],maxpars[ii]],
        'limited':[limitedmin[ii],limitedmax[ii]], 'fixed':fixed[ii],
        'parname':parnames[ii%3]+str(ii%3), 'error':ii} 
        for ii in xrange(len(params)) ]

    if veryverbose:
        print "GUESSES: "
        print "\n".join(["%s: %s" % (p['parname'],p['value']) for p in parinfo])

    mp = mpfit(mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet)
    mpp = mp.params
    mpperr = mp.perror
    chi2 = mp.fnorm

    if mp.status == 0:
        raise Exception(mp.errmsg)

    if not shh:
        print "Final fit values: "
        for i,p in enumerate(mpp):
            parinfo[i]['value'] = p
            print parinfo[i]['parname'],p," +/- ",mpperr[i]
        print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

    return mpp,n_gaussian(pars=mpp)(xax),mpperr,chi2


import numpy as np
import multiGaussFit
import matplotlib.pyplot as plt
import createGaussian

def norm(x, mean, sd):
  norm = []
  for i in range(x.size):
    norm += [1.0/(sd*np.sqrt(2*np.pi))*np.exp(-(x[i] - mean)**2/(2*sd**2))]
  return np.array(norm)

mean1, mean2, mean3 = 0, -2, 2
std1, std2, std3 = 0.5, 1 , 1.5

x = np.linspace(-20, 20, 500)
y_real = norm(x, mean1, std1) + norm(x, mean2, std2) + norm(x, mean3, std3)

nGauss=5
fitData = multiGaussFit.multigaussfit(x,y_real,ngauss=nGauss)

plt.plot(x, y_real, label='Real Data')
#plt.plot(x, y_init, 'r.', label='Starting Guess')
#plt.plot(x, y_est, 'g.', label='Fitted')
colors = ['g','r','c','m','y','k','b']
count = 0
compNum = 0
for i in xrange(nGauss):
    comp = modelData(x,
                     fitData[0][i+count],
                     fitData[0][i+count+1],
                     fitData[0][i+count+2])
    if comp.max() != 0:
        print 'yes!'
        plt.plot(x,comp, 
                 color=colors[i],linestyle='dashed',
                 label="Model Gaussian Fit, comp #" + str(compNum))
        compNum = compNum + 1
    count = count + 2

plt.legend()
plt.show()
