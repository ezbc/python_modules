
#!/usr/local/casapy-41.0.24668-001-64b/lib64/casapy/bin/python

################################################################################
# Script to image VLA/C, VLA/D, GMRT observations together
################################################################################

# Executed in CASA 3.4.0

# DOES NOT WORK IN CASA 4.0.0 or 4.1.0! Deconvolve components
# .residual and .restored are not written due to a bug

# Imaging directions:

#==============================================================================
#===============================================================================

# 1. Make a mosaics from the GMRT and VLA data separately.

# 2a. Combine the output PSFs with fluxscale as weights.

# 2b. Combine the dirty images (or residuals, if this is the second or
# later iteration) with fluxscale as weights.

# 3a. Deconvolve the combined dirty image for a fraction of the total
# number of desired iterations. This will make a model image
# (incomplete, but some number of components should be there).

# 3 b. If this is the second or later iteration, combine this model
# with the model from the previous loop.

# 4 Use this model to make residual mosaics as in (1) (input to clean
# as "modelimage").  There isn't any need to redo the PSF (2a).

# 5 Go to 2(b) with the residual mosaics.

#===============================================================================
#===============================================================================

################################################################################
# 1. Make a mosaics of each uv dataset separately.

def __clean(uvdir='./',uvdatasets=None,modelImage=None,iteration=0,uvtaper=''):
    from casa import clean

    # Determine whether the mosaic has been gridded yet
    image = []

    # If the file exists, don't clean
    for i, uvdata in enumerate(uvdatasets):
        filename = uvdata + '.' + str(iteration) + '.image'
        try:
            with open(filename + '/table.f0'):
                image.append(False)
                print 'File ' + filename + ' exists'
        except IOError:
            image.append(True)

    for i, uvdata in enumerate(uvdatasets):
        if image[i]:
            if modelImage==None:
                clean(vis=uvdir + uvdata,
                    imagename=uvdata + '.' + str(iteration),
                    mode='velocity',
                    outframe='lsrk',
                    start='190km/s',
                    width='1.8km/s',
                    nchan=100,
                    uvtaper=True,
                    robust=0.5,
                    outertaper=uvtaper,
                    #imagermode='mosaic',
                    #ftmachine='mosaic',
                    interpolation='nearest',
                    imsize=1000,
                    cell='2arcsec',
                    niter=0,
                    restfreq='1.420405752GHz')
            else:
                clean(vis=uvdir + uvdata,
                    imagename=uvdata + '.' + str(iteration),
                    modelimage=modelImage,
                    mode='velocity',
                    outframe='lsrk',
                    start='190km/s',
                    width='1.8km/s',
                    nchan=100,
                    uvtaper=True,
                    robust=0.5,
                    outertaper=uvtaper,
                    #imagermode='mosaic',
                    #ftmachine='mosaic',
                    interpolation='nearest',
                    imsize=1000,
                    cell='2arcsec',
                    niter=0,
                    restfreq='1.420405752GHz')

    # define psf image list
    psfImages = []
    fluxImages = []
    residualImages = []
    dirtyImages = []
    for i, uvdata in enumerate(uvdatasets):
        psfImages.append(uvdata + '.' + str(iteration) + '.psf')
        fluxImages.append(uvdata + '.' + str(iteration) + '.flux')
        residualImages.append(uvdata + '.' + str(iteration) + '.residual')
        dirtyImages.append(uvdata + '.' + str(iteration) + '.image')
    return psfImages,fluxImages,residualImages,dirtyImages

def __combineModels(modelImages,iteration=0):
    from casa import immath

    filename = 'modelCombine.' + str(iteration)
    try:
       with open(filename + '/table.f0'):
           combine = False
           print 'File ' + filename + ' exists'
    except IOError:
        combine = True

    if combine:
        immath(imagename=modelImages,
            outfile=filename,
            expr = 'IM0 + IM1')
    return filename

def __combinePSFs(psfImages,fluxImages):
    from casa import immath

    imageList = []
    num = '('
    den = '('
    for i in range(len(psfImages)):
        imageList.append(psfImages[i])
        imageList.append(fluxImages[i])
        tempNum = 'IM'+str(2*i)+'*IM'+str(2*i+1)
        tempDen = 'IM'+str(2*i+1)
        if i < len(psfImages) -1:
            num = num + tempNum + '+'
            den = den + tempDen + '+'
        else:
            num = num + tempNum + ')'
            den = den + tempDen + ')'
    expr = num + '/' + den

    #print imageList
    #print expr
    filename = 'psfCombine.im'
    try:
       with open(filename + '/table.f0'):
           combine = False
           print 'File ' + filename + ' exists'
    except IOError:
        combine = True

    if combine:
        immath(imagename=imageList,
           outfile=filename,
           expr=expr)
    return filename

def __combineResiduals(residualImages,fluxImages,iteration=0):

    from casa import immath
    imageList = []
    num = '('
    den = '('
    for i in range(len(residualImages)):
        imageList.append(residualImages[i])
        imageList.append(fluxImages[i])
        tempNum = 'IM'+str(2*i)+'*IM'+str(2*i+1)
        tempDen = 'IM'+str(2*i+1)
        if i < len(residualImages) -1:
            num = num + tempNum + '+'
            den = den + tempDen + '+'
        else:
            num = num + tempNum + ')'
            den = den + tempDen + ')'
    expr = num + '/' + den

    filename = 'residualCombine.' + str(iteration) + '.im'

    try:
       with open(filename + '/table.f0'):
           combine = False
           print 'File ' + filename + ' exists'
    except IOError:
        combine = True

    if combine:
        immath(imagename=imageList,
           outfile=filename,
           expr=expr)
    return filename

def __deconvolve(image=None,psf=None,niter=0,iteration=0,threshold='0.0mJy'):

    from casa import deconvolve
    filename = 'model.' + str(iteration)+ '.im'

    try:
       with open(filename + '/table.f0'):
           clean = False
           print 'File ' + filename + ' exists'
    except IOError:
        clean = True

    if clean:
        deconvolve(imagename=image,
                psf=psf,
                model=filename,
                niter=niter,
                threshold=threshold,
                #alg='multiscale',
                #scales=[0,10,50]
                )
    return filename

def heteroClean(uvdir='./',uvdatasets=None,niter=0,majorCycles=3,
        threshold='0.0mJy',imagename=None,uvtaper=''):

    import os
    if os.system('mkdir heteroCleanImages') != 256:
        os.chdir('heteroCleanImages')
    else:
        os.system('mkdir heteroCleanImages')
        os.chdir('heteroCleanImages')

    for i in range(majorCycles):

        if i==0:
            # 1. Make a mosaics of each uv dataset separately.
            psfImages,fluxImages,residualImages,dirtyImages = \
                    __clean(uvdir=uvdir,
                            uvdatasets=uvdatasets,
                            iteration=i,
                            uvtaper=uvtaper)
            # 2a. Combine the output PSFs with fluxscale as weights.
            psfImage = __combinePSFs(psfImages,
                    fluxImages)
        if i>0:
            if i==1:
                modelCombineImage = modelImage

            # 4 Use this model to make residual mosaics as in (1) (input to
            # clean as "modelimage").  There isn't any need to redo the PSF
            # (2a).
            psfImages,fluxImages,residualImages,dirtyImages = \
                __clean(uvdir=uvdir,
                        uvdatasets=uvdatasets,
                        modelImage=modelCombineImage,
                        iteration=i,
                        uvtaper=uvtaper)

        # 2b. Combine the dirty images (or residuals, if this is the second or
        # later iteration) with fluxscale as weights.
        #if i==0: images = dirtyImages
        #if i>0: images = residualImages
        #images = dirtyImages
        residualImage = __combineResiduals(images,
                fluxImages,
                iteration=i)

        # 3a. Deconvolve the combined dirty image for a fraction of the total
        # number of desired iterations. This will make a model image
        # (incomplete, but some number of components should be there).
        modelImage = __deconvolve(image=residualImage,
                 psf=psfImage,
                 niter=niter/majorCycles,
                 iteration=i,
                 threshold=threshold)

        if i>0 and i < majorCycles-1:
            # 3 b. If this is the second or later iteration, combine this model
            # with the model from the previous loop.
            modelCombineImage = __combineModels([modelImage,
                modelCombineImage],
                iteration=i)

        # After finished, move the final image to the original directory
        elif i == majorCycles-1:
            if imagename==None:
                imagename=modelImage
            os.system('mv ' + modelImage + '.restored ../' + imagename)

def execute():
    '''
    '''
    uvdir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
    uvdatasets = ['leopGMRT.contsub.normwts.ms',
                  'leopVLAc.contsub.ms',
                  'leopVLAd.contsub.ms']
    #print uvdatasets

    import os

    #    print 'temp 10'

    os.chdir('../temp12')

    heteroClean(uvdir=uvdir,
            uvdatasets=uvdatasets,
            niter=1000,
            majorCycles=3,
            threshold='2mJy',
            uvtaper='1arcsec',
            imagename='leop.vlacd.gmrt.1arcsec.image')

    os.chdir('../temp13')

    heteroClean(uvdir=uvdir,
            uvdatasets=uvdatasets,
            niter=1000,
            majorCycles=3,
            threshold='2mJy',
            uvtaper='20arcsec',
            imagename='leop.vlacd.gmrt.20arcsec.image')

    os.chdir('../temp14')

    heteroClean(uvdir=uvdir,
            uvdatasets=uvdatasets,
            niter=1000,
            majorCycles=3,
            threshold='2mJy',
            uvtaper='40arcsec',
            imagename='leop.vlacd.gmrt.40arcsec.image')

    # temp4 has 0 uv taper
    # temp5 has 1" uv taper 
    # temp6 has 40" uv taper
    # temp7 has 40" uv taper + mosaic imaging
    # temp8 has 20" uv taper + multiscale

    # Images using residual images
    # temp9 has 20" uv taper
    # temp10 has 1" uv taper
    # temp11 has 40" uv taper

    # Images using dirty images
    # temp12 1" uv taper 
    # temp13 20" uv taper
    # temp14 40" uv taper




