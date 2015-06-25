#!/usr/bin/python

import grid

if False:
    def test_spectralgrid_init():

        file_dir = '/d/bip3/ezbc/taurus/data/galfa/'
        file_name = 'taurus.galfa.cube.bin.4arcmin.fits'
        box = [300,300,325,325]
        noise_range = [(-110,-90),(90,110)]

        spGrid = grid.SpectralGrid(file_dir + file_name, box = box,
                noiseRange = noise_range,
                basesubtract = True)

        assert spGrid.cube is not None

if True:
    def test_spectralgrid_fit_profiles():

        import os
        reload(grid)
        file_dir = '/d/bip3/ezbc/taurus/data/galfa/'
        file_name = 'taurus.galfa.cube.bin.4arcmin.fits'
        box = [300,300,325,325]
        noise_range = [(-110,-90),(90,110)]

        spGrid = grid.SpectralGrid(file_dir + file_name, box = box,
                noiseRange = noise_range,
                basesubtract = True)

        guesses = [40,5,6]

        spGrid.fit_profiles(
                grow_pos = (312,312),
                guesses=guesses,
                ncomp=len(guesses)/3,
                alpha=0.001,
                coords='image',
                number_of_fits=10,
                verbose=True)

        os.system('rm -rf Test')






