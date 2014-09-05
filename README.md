##### Version: 0.0

# Planckpy Overview 

Extracts a region from the Planck Survey data archive in HEALPix format, and
converts to J2000 equatorial coordinates. To access the Planck data archive visit

[Planck archive](http://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/)


Original author: Dr. Robert Benjamin - benjamir@uww.edu

Maintained by: Elijah Bernstein-Cooper - ezbc@astro.wisc.edu

Please send emails to ezbc@astro.wisc.edu with questions, comments, or suggestions.

### Dependencies

The Planckpy module requires the following libraries:

+ [numpy](http://www.scipy.org/scipylib/download.html)

+ [astropy](http://www.astropy.org/)

+ [healpy](https://pypi.python.org/pypi/healpy)


### Installing

To add Planckpy module to your Python library, add the location of the Planckpy
module to your Python path. This line can be added to your .bashrc.

    $ export PYTHONPATH=$PYTHONPATH:<planckpy location>

A Python configuration setup script is in development to make installation
easier.

### Example Code

The following code will extract a region from the 857 GHz DR1 all-sky map
covering the Taurus molecular cloud. The example also shows how to write the
region to a FITS file.
```
>>> import planckpy as pl
>>> from astropy.io import fits
>>> (data, header) = pl.get_data(data_type='857',
                                 x_range=(21.3, 30.3),
                                 y_range=(60.0, 73.0))
>>> data.shape
(131, 91)
>>> header['TYPE']
'I_STOKES'
>>> fits.writeto('planck_region_857GHz.fits', 
                 data,
                 header=header)
```

The next example code will extract a region from the CO-Type1 J3-->2 DR1
all-sky map covering the Taurus molecular cloud.

```
>>> (data, header) = pl.get_data(data_type='CO-Type1', 
                                 x_range=(21.3, 30.3),
                                 y_range=(60.0, 73.0),
                                 field=8)
>>> fits.writeto('co_type1_j32.fits',
                 data,
                 header=header)

```