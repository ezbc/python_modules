#!/usr/bin/python

def modelData(x,amp,offset,width):

    profile = np.zeros(len(x))

    for i in xrange(len(x)):
        profile[i] = amp*np.exp(-(x[i]-offset)**2/(2*width**2))

    return profile
