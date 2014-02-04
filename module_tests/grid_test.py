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

        file_dir = '/d/bip3/ezbc/taurus/data/galfa/'
        file_name = 'taurus.galfa.cube.bin.4arcmin.fits'
        box = [300,300,325,325]
        noise_range = [(-110,-90),(90,110)]

        spGrid = grid.SpectralGrid(file_dir + file_name, box = box,
                noiseRange = noise_range,
                basesubtract = True)

        guesses = [40,5,6]

        spGrid.fit_profiles(
                growPos = (312,312),
                guesses=guesses,
                ncomp=len(guesses)/3,
                alpha=0.001,
                coords='image',
                number_of_fits=100,
                verbose=False,
                filename='Test')

        os.system('rm -rf Test')



'''
galfa = grid.SpectralGrid('/d/bip3/ezbc/taurus/data/galfa/' + \
        'taurus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=20.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

guesses = [3,-70,20, 5,-50,10, 5,-25,20, 4,-25,3, 3,-15,3, 40,5,6]

galfa.fit_profiles(
        growPos = (445,386),
        tileSaveFreq=1000,
        threshold = 1,
        filename='galfa.445_386.3CO_widths',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=100,
        COcube=cfa,
        COwidthScale=3.)

'''
