#!/usr/bin/python

import myplanck as pl
import pyfits as pf

longitude_range = (155, 165)
latitude_range = (-30,-15)


if 0:
    def test_857GHz():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location = data_location,
                data_type = '857', longitude_range = (longitude_range),
                latitude_range = (latitude_range), field = 0, resolution = 0.1,
                cut_last_pixel = False, verbose = True, dr_version = 1)

        # Write the data to FITS format
        pf.writeto(output_dir + '857GHz.fits', data, header = header, clobber
                = True, output_verify = 'fix')

if 0:
    def test_30GHz():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location = data_location,
                data_type = '030', longitude_range = (longitude_range),
                latitude_range = (latitude_range), field = 0, resolution = 0.1,
                cut_last_pixel = False, verbose = True, dr_version = 1)

        # Write the data to FITS format
        pf.writeto(output_dir + '30GHz.fits', data, header = header, clobber
                = True, output_verify = 'fix')

if 0:
    def test_tau353():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location=data_location,
                data_type='Dust Opacity', longitude_range = (longitude_range),
                latitude_range=(latitude_range), field=0, resolution = 0.1,
                cut_last_pixel = False, verbose = True, dr_version = 2)

        # Write the data to FITS format
        pf.writeto(output_dir + 'tau353.fits', data, header = header, clobber
                = True, output_verify = 'fix')

if 0:
    def test_ebv():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity', longitude_range = (145,150),
                latitude_range = (-10,-5), field = 2, resolution = 0.1,
                cut_last_pixel = False, verbose = True, dr_version = 2)

        # Write the data to FITS format
        pf.writeto(output_dir + 'ebv.fits', data, header = header, clobber
                = True, output_verify = 'fix')

if 0:
    def test_CO_type1_j10():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location=data_location,
                data_type='CO-Type1', longitude_range = (longitude_range),
                latitude_range=(latitude_range), field = 0, resolution = 0.1,
                cut_last_pixel = False, verbose = True)

        # Write the data to FITS format
        pf.writeto(output_dir + 'co_type1_j10.fits', data, header = header,
                clobber = True, output_verify = 'fix')

    def test_CO_type1_j32():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location=data_location,
                data_type='CO-Type1', longitude_range = (longitude_range),
                latitude_range=(latitude_range), field = 8, resolution = 0.1,
                cut_last_pixel = False, verbose = True)

        # Write the data to FITS format
        pf.writeto(output_dir + 'co_type1_j32.fits', data, header = header,
                clobber = True, output_verify = 'fix')

if 0:
    def test_CO_type2_j10():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location=data_location,
                data_type='CO-Type2', longitude_range = (longitude_range),
                latitude_range=(latitude_range), field = 4, resolution = 0.1,
                cut_last_pixel = False, verbose = True)

        # Write the data to FITS format
        pf.writeto(output_dir + 'co_type2_j10.fits', data, header = header,
                clobber = True, output_verify = 'fix')

    def test_CO_type2_j21():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location=data_location,
                data_type='CO-Type2', longitude_range = (longitude_range),
                latitude_range=(latitude_range), field = 4, resolution = 0.1,
                cut_last_pixel = False, verbose = True)

        # Write the data to FITS format
        pf.writeto(output_dir + 'co_type2_j21.fits', data, header = header,
                clobber = True, output_verify = 'fix')

if 0:
    def test_CO_type3():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type3', longitude_range = (longitude_range),
                latitude_range = (latitude_range), field = 0, resolution = 0.1,
                cut_last_pixel = False, verbose = True)

        # Write the data to FITS format
        pf.writeto(output_dir + 'co_type3.fits', data, header = header, clobber
                = True, output_verify = 'fix')

if 1:
    def test_CO_type3_no_reverse_xaxis():

        import os
        reload(pl)

        data_location = '/d/bip3/ezbc/planck/planck_raw_data/'
        output_dir = '/d/bip3/ezbc/planck/tests/'

        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type3', longitude_range = (longitude_range),
                latitude_range = (latitude_range), field = 0, resolution = 0.1,
                cut_last_pixel = False, verbose = True, reverse_xaxis=False)

        # Write the data to FITS format
        pf.writeto(output_dir + 'co_type3_no_reverse.fits', data, header =
                header, clobber = True, output_verify = 'fix')




