# Planckpy Overview 

Extracts a region from the Planck data in HEALPix format, and converts to
galactic coordinates.

Original author: Dr. Robert Benjamin - bobbenjamin@me.com

Edited and ported by: Elijah Bernstein-Cooper - ezbc@astro.wisc.edu

Code hosted at:
git@bitbucket.org:ezbc/planckpy.git

Module requires the following libraries:

+ [numpy](http://www.scipy.org/scipylib/download.html)

+ [pyfits](http://www.stsci.edu/institute/software_hardware/pyfits/Download)

+ [healpy](https://pypi.python.org/pypi/healpy)

### Some Example Code

    >>> import planckpy as pl
    >>> import pyfits as pf
    >>> (data, header) = pl.get_data(data_type = '857', longitude_range =
            (155,165), latitude_range = (-30, -15))
    >>> data.shape
    (151, 101)
    >>> header['TYPE']
    'I_STOKES'
    >>> pf.writeto('planck_region_857GHz.fits', data, header = header)


