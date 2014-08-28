#!/usr/bin/python
'''
Author: Elijah Bernstein-Cooper; ezbc@astro.wisc.edu

NAME
    hianalyis

FILE
    /home/elijah/python/myModules/hianalysis.py

DESCRIPTION
    HIanalysis
    ==========

    Provides
        1. Tool for blanking HI cubes by hand and creating moment maps.

    Execute in the CASA environment.

'''

def addNoise2Cube(cube, fwhm, std=None):
    '''

    Convolves random Gaussian noise with a psf. The convolved noise is then
    added to the cube.

    Parameters
    ----------
    cube : array_like
        3d array of data. Assumes first dimension is velocity axis
    fwhm : float
        FWHM (pixels) in each dimension.
        Single number to make all the same.
    std : float, optional
        Standard deviation (STD) of Gaussian noise to add to cube. If no value
        provided, STD will be calculated from the input cube.

    Returns
    -------
    out : ndarray
        Input cube with convolved Gaussian noise added.
    Examples
    --------
    '''

    # Import external moduls
    import numpy as np

    # Create PSF image of the cube
    psf = gauss2d(cube.shape[1:], fwhm, normalize=False)

    # Create cube of random noise
    if std is None:
        std = cube.std()
    noise = np.random.normal(0, std, cube.shape)

    # Convolve each plane of the noise cube with the psf
    convolNoise = np.zeros(noise.shape)
    for i in range(noise.shape[0]):
        plane = noise[i,:,:]
        newPlane = convolve(plane, psf)
        convolNoise[i,:,:] = newPlane
        #print convolNoise[i].max()

    # Compute the final noise cube
    def rms(x, axis=None):
        return np.sqrt(np.mean(x**2, axis=axis))

    cubeRms = rms(cube)
    noiseRms = rms(convolNoise)
    finalNoise = cubeRms * noise / noiseRms

    return cube + finalNoise

def blankcube(image,
              dummyMS,
              smooth=True,
              verbose=True,
              region='centerbox[[10h21m45s,18.05.14.9],[15arcmin,15arcmin]]',
              ruthless=False,
              extension='.image',
              beamround='int',
              blankThreshold=2.5,
              moments=[0]):

    '''
    Parameters
    ----------
    image : string
        Base name of image. If target image has extension other than '.image',
        change the extension keyword.

    dummyMS : string
        MS file required for blanking process in CASA

    smooth : bool, optional
        Smooth the image to a circular beam before blanking?
        Default True

    verbose : bool, optional

    region : string, optional
        region parameter featured in CASA

    ruthless : bool, optional
        Delete previous outputs from blankcube

    extension : string, optional
        Extension of the target image. Must include the '.', e.g., '.restored'
        Default '.image'

    beamround : string,float
        P

    blankThreshold : float, optional
        Initial blanking threshold of all pixels scaled by the standard
        deviation times the blankingthreshold
        Default = 2.5 (Walter et al. 2008)

    moments : list, optional
        Moments to calculate from cube. Options are 0,1,2.
        Example: [0,1,2]
        Default: [0]

    Returns
    -------
    out : null

    Examples
    --------

    '''

    from casa import immath,imsmooth,immoments
    from casa import image as ia
    from casa import imager as im
    from casa import quanta as qa
    import os
    import numpy as np

    # Delete files associated with previous runs
    if ruthless:
        os.system('rm -rf ' + image + '.smooth.blk.image')
        os.system('rm -rf ' + image + '.smooth.image')
        os.system('rm -rf ' + image + '.blk.image')
        os.system('rm -rf ' + image + '.mask')

    imageDir = './'

    # Create moment maps
    mom0,mom1,mom2 = False,False,False
    if len(moments) > 0:
        for i, moment in enumerate(moments):
            if moment == 0:
                mom0 = True
            if moment == 1:
                mom1 = True
            if moment == 2:
                mom2 = True
    if mom1 == True or mom2 == True:
        beamScale = 1.01
    elif mom0 == True:
        beamScale = 2.

    # determine beamsize of cube
    ia.open(imageDir + image + extension)
    beamsizes = np.zeros(ia.shape()[2])
    for i in range(ia.shape()[2]):
        beamsizes[i] = ia.restoringbeam(channel=i)['major']['value']
    beamsizeUnit = ia.restoringbeam(channel=i)['major']['unit']
    beamsize = qa.convert(str(beamsizes.max()) + beamsizeUnit,'arcsec')['value']
    if type(beamround) == float:
        beamsize_smooth = beamround*beamsize
    else:
        beamsize_smooth = np.ceil(beamsize)   #beamsize_smooth = 1.01 * beamsize
    ia.close()

    if verbose:
        print 'Max beamsize is ' + str(beamsize) + '"'

    if not os.path.exists(image + '.blk.image'):
        # create cube for blanking
        if smooth: # smooth to a larger beam if the user desires
            if verbose:
                print 'Convolving to ' + str(beamsize_smooth) + '"'

            imsmooth(imagename=image + extension,
                     outfile=image + '.blk.image',
                     major=str(beamsize_smooth) + 'arcsec',
                     minor=str(beamsize_smooth) + 'arcsec',
                     region=region,
                     pa='0deg',
                     targetres=True)
        else: # do no smooth
            immath(imagename=image + extension,
               outfile=image + '.blk.image',
                mode='evalexpr',
                region=region,
                expr='IM0')

    if not os.path.exists(image + '.smooth.image'):
        # convolve cube to 2X beam for blanking
        imsmooth(imagename=image + extension,
             outfile=image + '.smooth.image',
             major=str(beamsize*beamScale) + 'arcsec',
             minor=str(beamsize*beamScale) + 'arcsec',
             pa='0deg',
             region=region,
             targetres=True)

    # determine threshold of cube
    ia.open(image + '.smooth.image')
    threshold = ia.statistics()['sigma'][0] * blankThreshold
    ia.close()

    # blank the cube at threshold*sigma
    ia.open(image + '.smooth.image')
    ia.calcmask(mask=image + '.smooth.image > ' + str(threshold),
             name='mask1')
    wait = 'waits for calcmask to close'
    ia.close()

    # hand blank the cube
    im.open(dummyMS)
    pause = None
    while pause is None:
        im.drawmask(image=image + '.smooth.image',
                           mask=image + '.mask')
        pause = 0
    im.close

    # mask contains values of 0 and 1, change to a mask with only values of 1
    ia.open(image + '.mask')
    ia.calcmask(image + '.mask' + '>0.5')
    ia.close()

    # apply mask on smoothed image
    immath(imagename=image + '.smooth.image',
       outfile=image + '.smooth.blk.image',
       mode='evalexpr',
       mask='mask(' + image + '.mask)',
       expr='IM0')

    # apply mask on image
    ia.open(imageDir + image + '.blk.image')
    ia.maskhandler('copy',[image + '.smooth.blk.image:mask0', 'newmask'])
    ia.maskhandler('set','newmask')
    ia.done()

    cube = '.blk.image' # specify name of cube for moment calculation


    # create moment 0 map
    if mom0:
        if ruthless:
            os.system('rm -rf ' + image + '.mom0.image')
        immoments(imagename=image + cube,
                  moments=[0],
                  axis='spectra',
                  chans='',
                  mask='mask(' + image + cube + ')',
                  outfile=image + '.mom0.image')

    # create moment 1 map
    if mom1:
        if ruthless:
            os.system('rm -rf ' + image + '.mom1.image')
        immoments(imagename=image + cube,
                  moments=[1],
                  axis='spectra',
                  chans='',
                  mask='mask(' + image + cube + ')',
                  outfile=image + '.mom1.image')

    # create moment 2 map
    if mom2:
        if ruthless:
            os.system('rm -rf ' + image + '.mom2.image')
        immoments(imagename=image + cube,
                  moments=[2],
                  axis='spectra',
                  chans='',
                  mask='mask(' + image + cube + ')',
                  outfile=image + '.mom2.image')

    if verbose and mom0:
	from casa import imstat
        flux = imstat(image + '.mom0.image')['flux'][0]
        ia.open(image + '.mom0.image')
        beammaj = ia.restoringbeam(channel=0)['major']['value']
        beammin = ia.restoringbeam(channel=0)['minor']['value']
        beamsizeUnit = ia.restoringbeam(channel=0)['major']['unit']
        ia.close()
        print 'Moment Image: ' + str(image) + '.mom0.image'
        print 'Beamsize: ' + str(beammaj) + '" X ' + str(beammin) + '"'
        print 'Flux: ' + str(flux) + ' Jy km/s'

def convolve(image, psf, ft_psf=None, ft_image=None, no_ft=None,
             correlate=None, auto_correlation=None):
    '''
    NAME:
          CONVOLVE
    PURPOSE:
          Convolution of an image with a Point Spread Function (PSF)
    EXPLANATION:
          The default is to compute the convolution using a product of
          Fourier transforms (for speed).

    CALLING SEQUENCE:

          imconv = convolve( image1, psf, FT_PSF = psf_FT )
     or:
          correl = convolve( image1, image2, /CORREL )
     or:
          correl = convolve( image, /AUTO )

    INPUTS:
          image = 2-D array (matrix) to be convolved with psf
          psf = the Point Spread Function, (size < or = to size of image).

    OPTIONAL INPUT KEYWORDS:

          FT_PSF = passes out/in the Fourier transform of the PSF,
                  (so that it can be re-used the next time function is called).
          FT_IMAGE = passes out/in the Fourier transform of image.

          /CORRELATE uses the conjugate of the Fourier transform of PSF,
                  to compute the cross-correlation of image and PSF,
                  (equivalent to IDL function convol() with NO rotation of PSF)

          /AUTO_CORR computes the auto-correlation function of image using FFT.

          /NO_FT overrides the use of FFT, using IDL function convol() instead.
                  (then PSF is rotated by 180 degrees to give same result)
    METHOD:
          When using FFT, PSF is centered & expanded to size of image.
    HISTORY:
          written, Frank Varosi, NASA/GSFC 1992.
          Appropriate precision type for result depending on input image
                                  Markus Hundertmark February 2006
          Fix the bug causing the recomputation of FFT(psf) and/or FFT(image)
                                  Sergey Koposov     December 2006
    '''
    from numpy import *
    from numpy.fft import fft2, ifft2

    # begin
    n_params = 2
    psf_ft = ft_psf
    imft = ft_image
    noft = no_ft
    auto = auto_correlation

    sp = array(shape(psf_ft))
    sif = array(shape(imft))
    sim = array(shape(image))
    sc = sim / 2
    npix = array(image, copy=0).size

    if image.ndim!=2 or noft!=None:
       if (auto is not None):
          message("auto-correlation only for images with FFT", inf=True)
          return image
       else:
          if (correlate is not None):
             return convol(image, psf)
          else:
             return convol(image, rotate(psf, 2))

    if imft==None or (imft.ndim!=2) or imft.shape!=im.shape: #add the type check
       imft = ifft2(image)

    if (auto is not None):
       return roll(roll(npix * real(fft2(imft * conjugate(imft))), sc[0], 0),sc[1],1)

    if (ft_psf==None or ft_psf.ndim!=2 or ft_psf.shape!=image.shape or
             ft_psf.dtype!=image.dtype):
       sp = array(shape(psf))

       loc = maximum((sc - sp / 2), 0)         #center PSF in new array,
       s = maximum((sp / 2 - sc), 0)        #handle all cases: smaller or bigger
       l = minimum((s + sim - 1), (sp - 1))
       psf_ft = conjugate(image) * 0 #initialise with correct size+type according
       #to logic of conj and set values to 0 (type of ft_psf is conserved)
       psf_ft[loc[1]:loc[1]+l[1]-s[1]+1,loc[0]:loc[0]+l[0]-s[0]+1] = \
                      psf[s[1]:(l[1])+1,s[0]:(l[0])+1]
       psf_ft = ifft2(psf_ft)

    if (correlate is not None):
       conv = npix * real(fft2(imft * conjugate(psf_ft)))
    else:
       conv = npix * real(fft2(imft * psf_ft))

    sc = sc + (sim % 2)   #shift correction for odd size images.

    return roll(roll(conv, sc[0],0), sc[1],1)

def gauss2d(shape, fwhm, normalize=True):
    '''
    Author: P. L. Lim, Space Telescope Science Institute, Baltimore MD

    Parameters
    ----------
    shape : tuple, int
        Number of pixels for first and second dimension.
        Just one number to make all sizes equal.

    fwhm : float
        FWHM (pixels) in each dimension.
        Single number to make all the same.

    normalize : bool, optional
        Normalized so total PSF is 1.

    Returns
    -------
    psf : array_like
        Gaussian point spread function.

    Examples
    --------
    >>> from psf_gaussian import gauss2d
    >>> import matplotlib.pyplot as plt
    >>> psf = gauss2d(50, 2.5)
    >>> plt.imshow(psf, cmap=plt.cm.gray)
    '''

    import numpy as np

    if type(shape) is int or len(shape) == 1:
        xSize,ySize = shape,shape
    else:
        xSize,ySize = shape[0],shape[1]

    # Initialize PSF params
    xCentroid,yCentroid = (xSize - 1.0) * 0.5, (ySize - 1.0) * 0.5
    st_dev = 0.5 * fwhm / np.sqrt(2.0 * np.log(2))

    # Make PSF (Rene Breton 2011-10-20)
    # https://groups.google.com/group/astropy-dev/browse_thread/thread/5ee6cd662236e382
    x = np.indices([xSize, ySize])[0] - xCentroid
    y = np.indices([xSize, ySize])[0] - yCentroid
    psf = np.exp(-0.5 * ((x**2 + y**2) / st_dev**2))

    # Normalize
    if normalize:
        psf /= psf.sum()

    return psf

def flux2Mass(distance=None,distanceErr=None,arcsecPerPix=None,beamsize=None,
        flux=None,fluxErr=None):
    '''
    Distance in units of Mpc
    flux in units of Jy/beam km/s
    beamsize in units of "
    Returns mass in units of Msun
    '''

    import numpy as np

    beamarea = 1.13*beamsize**2/arcsecPerPix**2
    mass = 2.36e5 * distance**2 * flux / beamarea
    massErr = np.sqrt((mass * 2*distanceErr/distance)**2 + (fluxErr/flux)**2)
    return mass, massErr

def mass2SB(mass=None,massErr=None,area=None,arcsecPerPix=None,
        distance=None):
    '''
    mass in units of Msun
    massErr in units of Msun
    area in units of pixels
    arsecPerPix in units of "/pixel
    distance in units of Mpc
    '''
    pcPerArcsec = distance *1e6 / 206265.
    pcPerPixel = arcsecPerPix * pcPerArcsec
    pixSB = mass / area # units in Msun per pixel
    return pixSB / pcPerPixel**2, massErr/area / pcPerPixel**2

def make_SBimage(image='',cellsize=None,beamsize=None,distance=None,
                 extension='.image'):
    from casa import immath,imhead

    beamarea = 1.13*beamsize**2/cellsize**2
    massCoeff = 2.36e5 * distance**2 / beamarea
    pcPerArcsec = distance *1e6 / 206265.
    pcPerPixel = cellsize * pcPerArcsec
    pixSB = massCoeff / pcPerPixel**2 # units in Msun per pixel

    immath(imagename=image + extension,
       outfile=image + '.sb.image',
       mode='evalexpr',
       expr='IM0*' + str(pixSB))
    imhead(imagename=image + '.sb.image',
           mode='put',
           hdkey='bunit',
           hdvalue='',
           hdcomment='Units Msun/pc^-2')

def execute():
    import os
    os.chdir('/d/bip3/ezbc/leop/data/hi/casa/images/gmrt/')
    dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopGMRT.contsub.ms'
    blankcube('./','leop.gmrt.original.4arcsec',dummyMS=dummyMS,pbcor=False,
            ruthless=True,smooth=True)

def printProperties(filenames):

    from casa import imstat
    from casa import image as ia
    from casa import imager as im

    imageDir = './'

    for i, image in enumerate(filenames):
        flux = imstat(imageDir + image+ '.mom0.blk.image')['flux'][0]
        ia.open(imageDir + image + '.mom0.blk.image')
        beammaj = ia.restoringbeam(channel=0)['major']['value']
        beammin = ia.restoringbeam(channel=0)['minor']['value']
        beamsizeUnit = ia.restoringbeam(channel=0)['major']['unit']
        ia.close()

        ia.open(imageDir + image + '.image')
        noise = ia.statistics()['sigma'][0]
        ia.close()

        if True: mosaic,res = 'Multiscale clean', i
        #else: mosaic,res = 'Mosaic', i*3 -3
        print '\nImage: '+mosaic+' resolution #' + str(res) + \
                '\nBeamsize: ' + str(beammaj) + '" X ' + str(beammin) + '"' + \
                '\nStd: ' + str(noise) + ' Jy/Bm' + \
                '\nFlux: ' + str(flux) + ' Jy km/s'




