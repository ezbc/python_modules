#!/usr/bin/python

################################################################################
# Classes
################################################################################
import numpy as np

class SpectralGrid(object):
    ''' A class for decomposing a cube into Gaussians. :-)
    '''

    # Class variables
    gridlist = []
    cube = None
    spectralAxis = None
    raAxis = None
    decAxis = None
    gaussAreaList = []
    region = []
    header = None
    fitProfileCount = 0
    invalidProfileCount = 0

    def __init__(self,fitsfile,box=None,noiseScale=None,
            basesubtract=False,noiseRange=None,Tsys=None):
        """ Create an instance of SpectralGrid class. Requires an
        N-dimensional numpy array. Two dimensions of the array must be
        specified to build the Grid on. Creates a new Tile instance for each
        element in the 2-dimensional numpy array. The additional dimensions of
        the input array are stored in the Tile class instances.

        Parameters
        ----------
        fitsFile : string
            Path to fits file to load.
        gridfile : gridlist

        box : list
            [brc_x,brc_y,tlc_x,tlc_y]

        noiseScale : int
            Noise profile multiplied by this scale.

        noiseRange : tuple,list
            2D tuple or list specifying noise range.
            ((velMin1,velMax1),(velMin2,velMax2))

        Tsys : float
            System temperature of the telescope in units of K.
        """

        import pyfits as pf
        # Check optional keywords
        if fitsfile is None:
            print 'Please provide a fits image'
            sys.exit()

        if noiseScale is None:
            noiseScale = 1.

        if noiseRange is None:
            print 'Error: no velocity range selected for noise.'

        # Create a Cube object
        print ' Loading cube...'

        f = pf.open(fitsfile)
        h, self.cube = f[0].header, f[0].data

        self.header = h

        # Define axes arrays
        self.spectralAxis = make_velocityAxis(h)
        self.raAxis,self.decAxis = make_radecAxes(h)

        # Create Tiles for the grid
        tempList = []

        cubeXSize = self.cube.shape[1]
        cubeYSize = self.cube.shape[2]

        if box is None or len(box) != 4:
            box = [0,0,cubeXSize,cubeYSize]
        elif box[0] < 0:
            box[0] = 0
        elif box[1] < 0:
            box[1] = 0
        elif box[2] > cubeXSize:
            box[2] = cubeXSize
        elif box[3] > cubeYSize:
            box[3] = cubeYSize
        self.region = box

        print ' Loading profiles...'

        for i in range(box[0],box[2]):
            for j in range(box[1],box[3]):
                if self.cube[:,i,j].max() != self.cube[:,i,j].max():
                    profile = None
                    self.invalidProfileCount += 1
                else:
                    profile = self.cube[:,i,j]
                tempList.append([SpectralTile(gridXPos = i - box[0],
                                              gridYPos = j - box[1],
                                              imageXPos = i,
                                              imageYPos = j,
                                              box = box,
                                              profile = profile)])
            self.gridlist.append(tempList)
            tempList = []

        noise1 = noiseRange[0]
        noise2 = noiseRange[1]

        # Subtract baselines from each profile
        if basesubtract:
            print ' Subtracting baseline from profiles...'
            for i in range(self.get_xsize()):
                for j in range(self.get_ysize()):
                    self.get_tile(i,j).subtract_baseline(self.spectralAxis,
                                                         noise1,noise2)

        print ' Calculating noise characteristics...'

        # Calculate noises from each profile
        #print noiseScale
        self.calculate_profileNoises(noise1,noise2,noiseScale=noiseScale,
                                     Tsys=Tsys)

    def at_boundary(self,x,y,direction):
        """ Returns True boundary exists in the direction specified.

        Parameters
        ----------
        x : int
            Tile x position
        y : int
            Tile y position
        direction : str
            Direction to check boundary. Options are:
            'north','south','east', or 'west'
        """

        if self.get_neighbor(x,y,direction) is None:
            return True
        else:
            return False

    def calculate_profileNoises(self,velMin,velMax,noiseScale=None,Tsys=None):
        """ Calcuates noise envelopes for each pixel.
        """

        if noiseScale is None:
            noiseScale = 1.

        if Tsys is None:
            print 'System temperature not provided. Assuming Tsys = 30 K'
            Tsys = 30 # Tsys for Arecibo in K

        xarr = self.spectralAxis
        for i in range(self.get_xsize()):
            for j in range(self.get_ysize()):
                if self.get_tile(i,j).has_validProfile():
                    self.get_tile(i,j).calculate_noise(xarr,velMin,velMax,
                            noiseScale=noiseScale)
                    self.get_tile(i,j).calculate_noiseScale(Tsys)

    def check_residuals(self,x,y,threshold):
        """ Checks tile at x,y, to find any significant peaks in the
        residuals that should be fitted. Threshold is in units of standard 
        deviations.
        """

        from numpy import where

        noiseProfile = self.get_tile(x,y).get_noiseProfile()
        if self.get_tile(x,y).get_residuals() is not None:
            residuals = self.get_tile(x,y).get_residuals()
        else:
            residuals = self.get_tile(x,y).get_profile()

        maxResid = residuals.max()

        if noiseProfile is not None:
            if maxResid >= noiseProfile[maxResid == residuals] * threshold:
                return (maxResid,
                    self.spectralAxis[where(residuals == maxResid)[0][0]])
            else:
                return None

    def choose_neighbor(self,x,y,direction,COcube=None):
        ''' Determines whether neighboring tile is suitable for fitting.
        '''

        returnNeighbor = False # If all conditions satisfied, set to True
        # Check if tile exists
        if self.get_neighbor(x,y,direction) is not None:
            # Check if tile is at the boundary
            if not self.at_boundary(x,y,direction):
                # Check if tile profile has been fit
                if self.get_neighbor(x,y,direction).profile is not None:
                    # Check if tile is in queue already
                    if not self.get_neighbor(x,y,direction).is_visited():
                        # Does the CO tile exist?
                        if COcube is not None:
                            (imageX,imageY) = COcube.grid2image(x,y)
                            #hiX,hiY = self.grid2image(x,y)
                            #COx,COy = COcube.image2grid(hiX,hiY)
                            try:
                                COncomps = COcube.get_neighbor(imageY,
                                                        imageX,direction,
                                                        coords='image').ncomps
                                #COncomps = COcube.get_neighbor(COx,
                                #                        COy,direction).ncomps
                                if COncomps is not None:
                                    returnNeighbor = True
                            except AttributeError:
                                returnNeighbor = True
                        elif COcube is None:
                            returnNeighbor = True
        if returnNeighbor:
            self.get_neighbor(x,y,direction).set_visited()
            return self.get_neighbor(x,y,direction)
        else:
            return None

    def fit_profiles(self,guesses=None,ncomp=None,growPos=None,coords=None,
                    alpha=None,filename=None,tileSaveFreq=None,
                    threshold=None,numberOfFits=None,verbose=True,
                    outputFilename=None,COcube=None,COwidthScale=1.,
                    veryverbose=False):
        """ Fits spectral axis for each pixel in the cube.

        Parameters
        ----------
        guesses : list
            Default = [0,0,0]
            1D List of length ncomp * 3
            Contains initial guesses for the gaussian fits.
            Format: [shift,amplitude,width,shift,amplitude,width]

        ncomp : int
            Default = 1
            Number of gaussian components to fit.

        growPos : tuple
            2D tuple consisting of (x,y) positions for grow method to begin fitting

        coords : string, optional
            Default = 'grid'
            Either 'grid' or 'image' coordinates.

        alpha : float, optional
            Default = 0.05
            1 - alpha = confidence for F-test
            Rather, p-values below alpha are considered significant

        filename : string, optional
            Name of file to save SpectralGrid instance to.

        tileSaveFreq : int, optional
            After tileSaveFreq number of tiles fit, the SpectralGrid instance will
            be saved to filename.

        threshold : float, optional
            Default = 1
            When fitting additional components, only allow residuals to be
            considered significant if amplitude > threshold * noise in channel.

        verbose : boolean, optional
            Default = True

        outputFilename : string, optional
            Default is no output file will be written.

        COcube : SpectralGrid
            SpectralGrid instance of a CO on the same spatial grid as the input
            HI cube. Fitted profiles in the CO cube will be used to guesses
            where HI absorption is present in the HI cube. Channels in the HI
            cube with HI absorption will be excluded from a profile fit.

        COwidthScale : float
            Scalar of the CO velocity width. 
        """

        import Queue as queue
        from scipy.stats import f
        from mystats import ftest,fvalue
        from pickle import Pickler
        from sys import exit
        from agpy.gaussfitter import multigaussfit as gfit
        import random

        if growPos is None:
            growPos = (0,0)
        if ncomp is None:
            ncomp = 1
        if guesses is None:
            guesses = [0.1,0.1,0.1]
            ncomp = 1
        if alpha is None:
            alpha = 0.05
        if tileSaveFreq is None:
            tileSaveFreq = 10
        if threshold is None:
            threshold = 1

        # Set number of fits
        if COcube is None:
            totalTiles = (self.get_xsize())*(self.get_ysize())-1
        elif COcube.fitProfileCount < (self.get_xsize())*(self.get_ysize()):
            totalTiles = COcube.fitProfileCount - 1
        else:
            totalTiles = (self.get_xsize())*(self.get_ysize())-1

        if numberOfFits is None or numberOfFits > totalTiles:
            if COcube is None:
                numberOfFits = totalTiles - self.invalidProfileCount
            else:
                numberOfFits = totalTiles - self.invalidProfileCount

        # Convert coordinates to grid if image / set default
        if coords is None or coords is 'grid':
            try:
                x = growPos[0]
                y = growPos[1]
            except IndexError:
                print '\n Tile at grid position ' + str(x) + ', ' + str(y) + \
                        ' does not exist.'
                exit()

        if coords is 'image':
            try:
                x,y = self.image2grid(growPos[0],growPos[1])
            except IndexError:
                print '\n Tile at image position ' + str(x) + ', ' + str(y) + \
                        ' does not exist.'
                exit()
        if verbose:
            print 'Initial positions: x = ' + str(x) + ', y = ' + str(y)

        # Initialize arrays and counts
        tileQueue = queue.Queue()
        tileParamQueue = queue.Queue()
        tileNCompQueue = queue.Queue()
        tileCount = 0
        directions = ['east','west','north','south',]
        tileQueue.put(self.get_tile(x,y))
        tileParamQueue.put(guesses)

        while tileCount <= numberOfFits:
            if verbose:
                print '\nBeginning fit #' + str(tileCount) + ' of ' + \
                        str(numberOfFits)

            # Get the next tile, but if empty the fitting is complete
            if not tileQueue.empty():
                tile = tileQueue.get()
                if tileCount != 0:
                    try:
                        guesses = self.writeToList(tileParamQueue.get())
                    except TypeError:
                        guesses = []
                else:
                    guesses = tileParamQueue.get()
                    #####################
                    #guesses = []
                    ###################

            else:
                print 'No further profiles to fit.'
                break


            # Perform first fit
            if tileCount == 0:
                nCompInit = len(guesses) / 3
            else:
                nCompInit = tileNCompQueue.get()

            ##############
            #nCompInit = 0
            ##############


            x = tile.get_position()[0]
            y = tile.get_position()[1]
            (imageX,imageY) = self.grid2image(x,y)

            # Get the CO tile
            if COcube is not None:
                COtile = COcube.get_tile(imageX,imageY,coords='image')
            else:
                COtile = None

            if verbose:
                print 'Fitting tile at Grid position ' + \
                        str(x) + ', ' + str(y)
                print 'Fitting tile at Image position ' + \
                        str(self.grid2image(x,y)[0]) + \
                            ', ' + str(self.grid2image(x,y)[1])

            #print 'nCompInit = ' + str(nCompInit)
            #print 'initial guesses = ' + str(guesses)

            # Fitting begins:
            # Steps to take:
            # Perform fit with previous pixel's parameters as initial guesses
            # Perform fit with same parameters - 
                # parameters for smallest gaussian
                # Calculate p-value between 1st and 2nd fits
            # Perform fit with previous pixels's parameters as intial guesses +
                # guess from highest residual
                # Calculate p-value betwen 1st and 3rd fits
            # If p-value is signifcant, choose fit with smaller p-value
            # Continue performing f-test with 
                # either fewer components or more components
                # until p-value is insignificant.

            # Define significant p-value
            pSignificant = alpha

            # Perform fit with initial # of components
            # ----------------------------------------
            # First create a temporary tile
            if verbose:
                print 'Initial # of components = ' + str(nCompInit)
            tempTile_init = SpectralTile(tile=tile)

            if nCompInit == 0:
                profile = tempTile_init.profile
                noiseProfile = tempTile_init.noiseProfile
                chiInit = sum((profile - noiseProfile)**2 / noiseProfile)
                initNComps = 0
                initParams = []
            else:
                tempTile_init.fit_profile(self.spectralAxis,
                    guesses=guesses,
                    ncomp=nCompInit,
                    COtile=COtile,
                    COwidthScale=COwidthScale)
                tempTile_init.set_visited()
                chiInit = tempTile_init.get_chi2()
                initNComps = tempTile_init.get_ncomps()
                initParams = tempTile_init.get_params()
                initCompAreaList = tempTile_init.get_compAreaList()

            # Change the first tile guesses to the fit parameters
            if tileCount == 0:
                guesses = self.writeToList(tempTile_init.params)

            # Perform fit with one fewer components
            # -------------------------------------
            if initNComps > 1:
                # First create a temporary tile
                tempTile_minus = SpectralTile(tile=tile)

                # Create initial guesses
                newParams = initParams[1:]
                guesses = self.writeToList(newParams)

                # Perform the fit
                tempTile_minus.fit_profile(self.spectralAxis,
                    guesses=guesses,
                    ncomp=initNComps - 1,
                    COtile=COtile,
                    COwidthScale=COwidthScale)

                chiMinus = tempTile_minus.get_chi2()
            else:
                profile = tempTile_init.profile
                noiseProfile = tempTile_init.noiseProfile
                chiMinus = sum((profile - noiseProfile)**2 / noiseProfile)

            # Perform fit with an additional component
            # ----------------------------------------     
            # First create a temporary tile

            tempTile_plus = SpectralTile(tile=tile)
            self.gridlist[x][y][0] = tempTile_init

            peakResidual = self.check_residuals(x,y,threshold)
            if peakResidual is not None:
                if initNComps == 0 or initNComps == 1:
                    self.get_tile(x,y).residuals = self.get_tile(x,y).profile
                # Guess next component parameters
                guesses = self.writeToList(initParams)
                guesses.append(peakResidual[0])
                guesses.append(peakResidual[1])
                width = self.guess_width(x,y,peakResidual[1])
                guesses.append(width) # guess of 5 km/s width

                # Perform fit
                tempTile_plus.fit_profile(self.spectralAxis,
                        guesses=guesses,
                        ncomp=initNComps+1,
                        COtile=COtile,
                        COwidthScale=COwidthScale)

                chiPlus = tempTile_plus.get_chi2()
            elif peakResidual is None:
                chiPlus = chiInit

            # Perform F-test on three fits
            # ----------------------------------------     
            dofInit = len(self.spectralAxis) - initNComps * 3
            dofMinus = dofInit + 3
            dofPlus = dofInit - 3
            pvaluePlus = ftest(chiInit,chiPlus,dofInit,dofPlus)
            pvalueMinus = ftest(chiMinus,chiInit,dofMinus,dofInit)

            veryverbose=False
            if veryverbose:
                print 'chiPlus: ' + str(chiPlus)
                print 'chiInit: ' + str(chiInit)
                print 'chiMinus: ' + str(chiMinus)
                print 'pvalue-: ' + str(pvalueMinus)
                print 'pvalue+: ' + str(pvaluePlus)

            significant = False
            if pvaluePlus < pSignificant or pvalueMinus > pSignificant:
                significant = True

            # Perform fits with more components
            # ---------------------------------    
            if significant and pvaluePlus < 1 - pvalueMinus:
                if initNComps == 0:
                    tempTile_init.ncomps = 0
                    tempTile_init.residuals = tempTile_init.profile
                    tempTile_init.params = []
                    tempTile_init.visited = True
                    tempTile_init.fitSuccess = True
                    tempTile_init.fit_chi2 = chiMinus

                while pvaluePlus < pSignificant:
                    if verbose:
                        print 'Additional components significant.'
                    # Create temporary tile
                    tempTile_plus = SpectralTile(tile=tempTile_init)

                    # Get initial parameters
                    chiInit = tempTile_init.get_chi2()
                    initNComps = tempTile_init.get_ncomps()
                    initParams = tempTile_init.get_params()

                    self.gridlist[x][y][0] = tempTile_init
                    peakResidual = self.check_residuals(x,y,1)
                    if peakResidual is not None:
                        # Write new guess from largest residual
                        guesses = self.writeToList(initParams)
                        guesses.append(peakResidual[0])
                        guesses.append(peakResidual[1])
                        width = self.guess_width(x,y,peakResidual[1])
                        guesses.append(width) # guess of 5 km/s width

                        # Fit the new tile                    
                        tempTile_plus.fit_profile(self.spectralAxis,
                                guesses=guesses,
                                ncomp=initNComps+1,
                                COtile=COtile,
                                COwidthScale=COwidthScale)

                        chiPlus = tempTile_plus.get_chi2()

                    if verbose:
                        print 'Number of components in fit = ' + \
                            str(tempTile_plus.get_ncomps())

                    # Calculate p-value            
                    dofInit = len(self.spectralAxis) - initNComps * 3
                    dofPlus = dofInit - 3
                    tempTile_init.fValuePlusList.append(fvalue(chiInit,chiPlus,
                                                                dofInit,dofPlus))
                    pvaluePlus = ftest(chiInit,chiPlus,dofInit,dofPlus)
                    if veryverbose:
                        print 'pvalue = ' + str(pvaluePlus)
                    # If p-value is still significant, use new tile and repeat
                    tempTile_init = tempTile_plus

            # Perform fits with fewer components
            # ----------------------------------
            if significant and 1 - pvalueMinus <= pvaluePlus:
                pvaluePlus = pvalueMinus #extract p-value from previous fit
                pvalueMinus = 0

                if initNComps <= 1:
                    if verbose and initNComps == 1:
                            print 'Fewer components significant'
                    tempTile_init.ncomps = 0
                    tempTile_init.residuals = tempTile_init.profile
                    tempTile_init.params = []
                    tempTile_init.visited = True
                    tempTile_init.fitSuccess = True
                    tempTile_init.fit_chi2 = chiMinus

                else:
                    while pvaluePlus > pvalueMinus:
                        if verbose:
                            print 'Fewer components significant'

                        # Get initial parameters
                        chiInit = tempTile_init.get_chi2()
                        initNComps = tempTile_init.get_ncomps()
                        initParams = tempTile_init.get_params()

                         # Create temporary tile
                        tempTile_minus = SpectralTile(tile=tempTile_init)

                        # Use params from previous tile for guesses
                        newParams = initParams[1:]
                        guesses = self.writeToList(newParams)

                        # Fit the new tile
                        tempTile_minus.fit_profile(self.spectralAxis,
                                guesses=guesses,
                                ncomp=initNComps - 1,
                                COtile=COtile,
                                COwidthScale=COwidthScale)

                        if verbose:
                            print 'Number of components in fit = ' + \
                                    str(tempTile_minus.get_ncomps())

                        # Calculate p-value
                        chiMinus = tempTile_minus.get_chi2()
                        dofInit = len(self.spectralAxis) - initNComps*3
                        dofMinus = dofInit + 3

                        # Append fvalue to list
                        tempTile_init.fValueMinusList.append(
                                fvalue(chiMinus,chiInit,dofMinus,dofInit))

                        pvalueMinus = ftest(chiMinus,chiInit,dofMinus,dofInit)

                        if veryverbose:
                            print 'pvalue = ' + str(pvalueMinus)
                        #print 'chiMinus : ' + str(chiMinus)
                        #print 'chiInit : ' + str(chiInit)
                        #print 'dofMinus : ' + str(dofMinus)
                        #print 'dofInit : ' + str(dofInit)

                        # If p-value is still significant, use new tile and
                        # repeat
                        tempTile_init = tempTile_minus
                        pvaluePlus = pvalueMinus

            # If no fits were significant, then adjust the parameters
            if initNComps == 0:
                tempTile_init.ncomps = 0
                tempTile_init.residuals = tempTile_init.profile
                tempTile_init.params = []
                tempTile_init.visited = True
                tempTile_init.fitSuccess = True
                tempTile_init.fit_chi2 = chiMinus


            # Write best-fit tile to the grid
            self.gridlist[x][y][0],tile = tempTile_init,tempTile_init

            # add fit to total fit count
            self.fitProfileCount += 1

            if verbose:
                print 'Tile at ' + str(x) + ', ' + str(y) + ' fitted with ' + \
                        str(tile.ncomps) + ' components'

            # Write grid to file
            if filename is not None and tileCount%tileSaveFreq == 0:
                if verbose:
                    print '\n Saving SpectralGrid instance as ' + str(filename)
                write_grid(self,filename)
            elif filename is not None and tileCount == numberOfFits:
                if verbose:
                    print '\n Saving SpectralGrid instance as ' + str(filename)
                write_grid(self,filename)

            # Now add new tiles to list
            # MUST randomize directions so that initial conditions do not
            #   affect the output!!!!!
            random.shuffle(directions)
            for i, direction in enumerate(directions):
                # Check if tile exists
                neighbor = self.choose_neighbor(x,y,direction,COcube=COcube)
                if neighbor is not None:
                    tileQueue.put(self.get_neighbor(x,y,direction,
                            verbose=verbose))
                    tileNCompQueue.put(tile.get_ncomps())
                    tileParamQueue.put(tile.params)

            tileCount = tileCount + 1

    def fit_profile(self,guesses,ncomp,x,y,coords='grid',COcube=None,
                    COwidthScale=1.):
        """ Fits gaussians to a pixel in the cube.
        """

        # Convert coordinates to grid if image / set default
        if coords is 'image':
            try:
                x,y = self.image2grid(x,y)
            except IndexError:
                print '\n Tile at image position ' + str(x) + ', ' + str(y) + \
                        ' does not exist.'
                exit()

        if COcube is not None:
            COtile = COcube.get_tile(x,y)
        else:
            COtile = None

        self.get_tile(x,y).fit_profile(self.spectralAxis,
                                           guesses=guesses,
                                           ncomp=ncomp,
                                           COtile=COtile,
                                           COwidthScale=COwidthScale)

    def get_component(self,x,y,compNumber):
        """ Returns a gaussian component from a fit.
        """

        tile = self.get_tile(x,y)
        return tile.make_component(compNumber,
                                   self.spectralAxis)
    def get_imagexsize(self):
        return len(self.cube[0][:])

    def get_imageysize(self):
        return len(self.cube[0][0])

    def get_neighbor(self,x,y,direction,coords='grid',verbose=True):
        """ Gets a neighboring tile East, West, North, or South of the
        reference tile. The North-West corner of the grid has coordinates 
        (x,y) = (0,0).
        """
        # Convert coordinates to grid if image / set default
        if coords is 'image':
            try:
                x,y = self.image2grid(x,y)
            except IndexError:
                if verbose:
                    print '\n get_neighbor: Tile at image position ' \
                            + str(x) + ', ' + str(y) + ' does not exist.'

        if direction is 'south' and y != 0:
            yNeighbor = y - 1
            return self.get_tile(x,yNeighbor)
        if direction is 'north' and y != self.get_ysize():
            yNeighbor = y + 1
            return self.get_tile(x,yNeighbor)
        if direction is 'west' and x != 0:
            xNeighbor = x - 1
            return self.get_tile(xNeighbor,y)
        if direction is 'east' and x != self.get_xsize():
            xNeighbor = x + 1
            return self.get_tile(xNeighbor,y)

    def get_spectralAxis(self):
        """ Returns the spectral axis of the cube.
        """

        return self.spectralAxis

    def get_tile(self,x, y, coords = None):
        """ Returns a tile at position x,y.

        Parameters
        ----------
        x,y : int
            Positions
        coords : string
            Default = 'grid'
            Either 'grid' or 'image' coordinates.
        """

        if coords is None or coords is 'grid':
            try:
                return self.gridlist[x][y][0]
            except IndexError:
                print '\n Tile at grid position ' + str(x) + ', ' + str(y) + \
                        ' does not exist.'

        if coords is 'pixel':
            try:
                x,y = self.image2grid(x,y)
                return self.gridlist[x][y][0]
            except IndexError:
                print '\n Tile at image position ' + str(x) + ', ' + str(y) + \
                        ' does not exist.'

        if coords is 'image':
            try:
                x,y = self.image2grid(x,y)
                return self.gridlist[x][y][0]
            except IndexError:
                print '\n Tile at image position ' + str(x) + ', ' + str(y) + \
                        ' does not exist.'


    def get_xsize(self):
        return len(self.gridlist[:])

    def get_ysize(self):
        return len(self.gridlist[0])

    def grid2image(self,gridXPos,gridYPos):
        """ Converts grid positions to image pixel positions.
        """

        return (gridXPos + self.region[0], gridYPos + self.region[1])

    def guess_width(self,x,y,maxpos):
        ''' Guesses the width of a residual peak.
        maxpos : int
            In velocity units.
        '''
        from numpy import where

        tile = self.get_tile(x,y)
        residuals = tile.get_residuals()
        # Find position in array of maximum velocity position
        maxgridpos = where(self.spectralAxis == maxpos)[0][0]
        # Find amplitude of maximum residual
        maxResid = residuals[maxgridpos]
        # Define width as maxpos minus position where residual is 0.5 maxpos
        try:
            minPoses = where(residuals[maxgridpos:-1] < 0.5 * maxResid)[0]
            widthPos = abs(min(minPoses) - maxgridpos)
            width = abs(self.spectralAxis[widthPos] - \
                    self.spectralAxis[maxgridpos])
        except ValueError:
            minPoses = where(residuals[0:maxgridpos] < 0.5 * maxResid)[0]
            widthPos = abs(min(minPoses) - maxgridpos)
            width = abs( - self.spectralAxis[widthPos] + \
                    self.spectralAxis[maxgridpos])
        return width

    def image2grid(self,imageXPos,imageYPos):
        """ Converts image positions to grid pixel positions.
        """

        return (imageXPos - self.region[0], imageYPos - self.region[1])

    def pix2image(self,pixXPos,pixYPos):
        """ Converts grid positions to image pixel positions.
        """

        return (self.raAxis[pixXPos],self.decAxis[pixYPos])


    def writeToList(self,inputList):
        """ Returns a list of lists as a list.
        """

        newList = []

        if inputList is not None:
           for i in xrange(0,len(inputList)):
            for j in xrange(0,len(inputList[0])):
                newList.append(inputList[i][j])

        return newList

class SpectralTile(object):
    """ Subclass of Tile.
    """

    guesses = None
    ncomps = None
    params = None
    residuals = None
    gridXPos = None
    gridYPos = None
    imageXPos = None
    imageYPos = None
    xOffset = None
    yOffset = None
    visited = False
    fitSuccess = False
    profile = None
    noise = None
    noiseProfile = None
    compAreaList = []
    fit_chi2 = None
    fValuePlusList = []
    fValueMinusList = []
    spectralAxis = None
    coVelocities = None

    def __init__(self,gridXPos = None,gridYPos = None,
                 imageXPos = None, imageYPos = None,
                 profile = None, tile = None,
                 box = None):
        """ Create a tile use a pyspeckit Cube instance.

        tile : SpectralTile
            Copies parameters from input tile to new tile.
        """
        import numpy as np

        if tile is None:
            self.guesses = None
            self.ncomps = None
            self.params = None
            self.gridXPos = gridXPos
            self.gridYPos = gridYPos
            self.imageXPos = imageXPos
            self.imageYPos = imageYPos
            self.xOffset = box[0]
            self.yOffset = box[1]
            self.profile = profile
            self.fitSuccess = False
            self.visited = False

        # Copy tile parameters to new tile
        if tile is not None:
            self.guesses = tile.guesses
            self.ncomps = tile.ncomps
            self.params = tile.params
            self.residuals = tile.residuals
            self.visited = tile.visited
            self.fitSuccess = tile.fitSuccess
            self.profile = tile.profile
            self.noise = tile.noise
            self.noiseProfile = tile.noiseProfile
            self.compAreaList = tile.compAreaList
            self.fit_chi2 = tile.fit_chi2
            self.gridXPos = tile.gridXPos
            self.gridYPos = tile.gridYPos
            self.imageXPos = tile.imageXPos
            self.imageYPos = tile.imageYPos
            self.xOffset = tile.yOffset
            self.yOffset = tile.xOffset
            self.fValuePlusList = tile.fValuePlusList
            self.fValueMinusList = tile.fValueMinusList
            self.spectralAxis = tile.spectralAxis
            self.coVelocities = tile.coVelocities



    def calculate_compArea(self,compNumber,spectralAxis):
        """ Caclutes area under gaussian.
        """

        from scipy.integrate import quad

        a = spectralAxis.min()
        b = spectralAxis.max()

        params = self.get_compParams(compNumber)

        return quad(self.gaussian,a,b,args=(params[0],params[1],params[2]))[0]

    def calculate_compAreaList(self,spectralAxis):
        """ Calculates areas of gaussian components in profile fit.
        """

        compAreaList = []

        for i in xrange(0,self.ncomps):
            area = self.calculate_compArea(i,spectralAxis)
            compAreaList.append(area)

        self.compAreaList = sorted(compAreaList)

        tempParams = []
        # Sort gaussian parameters by area
        for i in xrange(0,self.ncomps):
            newParams = self.params[compAreaList.index(sorted(compAreaList)[i])]
            tempParams.append(newParams)

        self.params = tempParams

    def calculate_noise(self,xarr,lowrange,highrange):
        """ Calculates rms noise of Tile profile given velocity ranges.
        """
        velMin = [lowrange[0],highrange[0]]
        velMax = [lowrange[1],highrange[1]]

        std = 0
        if self.profile is not None:
            for k in xrange(len(velMin)):
                noiseRegion = np.where((xarr >= velMin[k]) & \
                                (xarr <= velMax[k]))
                std = np.std(self.get_profile()[noiseRegion]) + std

            std = std / len(velMin)
            self.noise = std
        else:
            self.noise = np.NaN

    def calculate_noiseScale(self,Tsys):
        """ Creates an array for scaling the noise by (Tsys + Tb) / Tsys
        """

        if self.has_validProfile():
            profile = self.get_profile()
            n = np.zeros(len(profile))

            for i, Tb in enumerate(profile):
                n[i] = (Tsys + Tb) / Tsys * self.noise

            self.noiseProfile = n

    def fit_profile(self,spectralAxis,guesses=None,ncomp=None,COtile=None,
                    COwidthScale=1.):
        """ Perform a gaussian decomposition on the profile.
        """

        from agpy.gaussfitter import multigaussfit as gfit
        from numpy import where,zeros,linspace,array,unique
        from numpy import concatenate as concat
        import numpy as np

        self.guesses = guesses
        # define mask for spectral axis
        mask = np.ones(spectralAxis.shape,dtype=bool)

        # Determine extent of CO
        if COtile is not None:
            if COtile.ncomps is not None and COtile.ncomps > 0:
                #vels = zeros(len(COtile.params))
                #widths = zeros(len(COtile.params))
                ## Get velocities and widths of components
                #for i in range(len(COtile.params)):
                #    vels[i] = COtile.params[i][1]
                #    widths[i] = COtile.params[i][2]
                vels = []
                widths = []
                for i in range(len(COtile.params)):
                    if COtile.params[i][0] > 3*COtile.noise:
                        vels.append(COtile.params[i][1])
                        widths.append(COtile.params[i][2])
                vels = np.asarray(vels)
                widths = np.asarray(widths)
                # lower and higher velocity extent of components:
                minVels = vels - COwidthScale*widths
                maxVels = vels + COwidthScale*widths
                self.coVelocities = []
                velpos = array([])

                for i in range(len(minVels)):
                    # Define pixels excluded from CO region
                    velpos_temp = where((spectralAxis > minVels[i]) & \
                                        (spectralAxis < maxVels[i]))[0]

                    mask[velpos_temp] = 0

                    # Define velocities included in CO region
                    coVels = spectralAxis[(spectralAxis > minVels[i]) & \
                                          (spectralAxis < maxVels[i])]
                    self.coVelocities.append(coVels)

                # now determine unique pixels
                # velpos = unique(velpos).astype(int)
                # Remove profile pixels within CO region
                xarray = spectralAxis[mask]
                profile = self.profile[mask]
                noiseProfile = self.noiseProfile[mask]
            else:
                xarray = spectralAxis
                profile = self.profile
                noiseProfile = self.noiseProfile
                velpos=linspace(1,len(profile),len(profile))
        else:
            xarray = spectralAxis
            profile = self.profile
            noiseProfile = self.noiseProfile
            velpos=linspace(1,len(profile),len(profile))

        # Fit the profile
        if guesses != [0,0,0]:
            # force high amplitude guesses to a threshold
            for i in range(len(guesses)/3):
                if guesses[i*3] > 1.1*self.profile.max():
                    guesses[i*3] = self.profile.max()
            fit = gfit(xarray,
                    profile,
                    ngauss=ncomp,
                    err=noiseProfile,
                    params=self.guesses,
                    limitedmin=[True,True,True]*ncomp,
                    limitedmax=[True,True,False]*ncomp,
                    minpars=[0,spectralAxis.min(),0],
                    maxpars=[1.1*self.profile.max(),spectralAxis.max(),0])

            self.fitSuccess = True
        else:
            self.residuals = self.profile
            self.fitSuccess = False

        # Set class variables
        if self.fitSuccess:
            self.residuals = zeros(len(self.profile))
            for i in range(len(profile)):
                if mask[i]:
                    self.residuals[i] = profile[i] - fit[1][i]
            tempParams = fit[0]
            self.params = []
            self.ncomps = 0
            self.fit_chi2 = fit[3]
            self.spectralAxis = spectralAxis

            # Write gaussian component parameters to Tile, ignore comps with
            # amp=0 or width=0
            for i in xrange(0,ncomp):
                if tempParams[i*3] > self.noise * 0. \
                        and tempParams[i*3 + 2] != 0:
                    self.params.append(tempParams[i*3:i*3+3])
                    self.ncomps = self.ncomps + 1

            self.calculate_compAreaList(spectralAxis)


    def gaussian(self,x,amp,shift,width):
        """ Returns gaussian with amp, shift, width, at position x.
        """

        return amp * np.exp( - (x - shift)**2 / (2. * width**2))

    def get_chi2(self):
        return self.fit_chi2

    def get_compAreaList(self):
        return self.compAreaList

    def get_compParams(self,compNumber):
        """ Extracts parameters of fit from compNumber.
        Parameters are returned as a list in the format:
        [amplitude,shift,width]
        """

        if self.params is not None:
            if compNumber > self.ncomps:
                print 'Number of components in spectrum is ' + str(self.ncomps)
                print '\n' + 'Please choose a valid component.'
            else:
                return self.params[compNumber]

    def get_guesses(self):
        return self.guesses


    def get_ncomps(self):
        """ Returns number of components in a fit.
        """

        return self.ncomps

    def get_noise(self):

        return self.noise

    def get_noiseProfile(self):

        return self.noiseProfile

    def get_params(self):

        return self.params

    def get_position(self):
        """ Returns a tuple of the x and y position of the Tile.
        """

        return (self.gridXPos,self.gridYPos)

    def get_profile(self):

        return self.profile

    def get_residuals(self):
        """ Returns residuals.
        """

        return self.residuals

    def get_componentVelocities(self):
        """ Returns velocities of components.
        """

        vels = np.zeros(self.ncomps)
        for i in range(0,self.ncomps):
            vels[i] = self.params[i][1]
        return vels

    def has_validProfile(self):
        """ Determines whether a profile exists for the tile.
        """

        if self.profile is None:
            return False
        else:
            return True

    def is_fitSuccessful(self):
        """ Returns whether the tile profile fit was successful.
        """

        return self.fitSuccess

    def is_visited(self):
        """ Returns whether a gaussian has been fit to the tile profile.
        """

        return self.visited

    def make_component(self,compNumber,spectralAxis):
        """ Creates a component using fit parameters.
        Returns a 1D numpy array.
        """

        if self.params is None:
            print '\n Warning: Fit has not been performed, \n' + \
                ' or fit has failed. \n' + \
                ' No parameters are available.'
        else:
            amp = self.params[compNumber][0]
            shift = self.params[compNumber][1]
            width = self.params[compNumber][2]

            return self.make_gaussian(spectralAxis,amp,shift,width)

    def make_gaussian(self,xarr,amp,shift,width):
        """ Creates a gaussian shape over across an axis using gaussian
        parameters.
        """

        gaussArray = np.zeros(len(xarr))

        for i, x in enumerate(xarr):
            gaussArray[i] = self.gaussian(x,amp,shift,width)

        return gaussArray

    def set_noise(self,noise):
        """ Sets noise of profile.
        """

        self.noise = noise

    def set_visited(self):
        """ Sets the Tile as visited.
        """

        self.visited = True

    def subtract_baseline(self,xarr,lowrange,highrange):
        """ Subtracts a first-degree polynomial baseline from the profile.
        """

        from pylab import polyfit

        if self.profile is not None:
            tempProfile = self.profile

            fitRegion = []
            velMin = [lowrange[0],highrange[0]]
            velMax = [lowrange[1],highrange[1]]

            for k in xrange(len(velMin)):
                tempRegion = np.where((xarr >= velMin[k]) & \
                        (xarr <= velMax[k]))[0]
                for i in range(len(tempRegion)):
                    fitRegion.append(tempRegion[i])

            m,b = polyfit(xarr[fitRegion], self.profile[fitRegion], 1)

            self.profile = tempProfile - (m * xarr + b)

    def calculate_noise(self,xarr,lowrange,highrange,noiseScale=None):
        """ Calculates rms noise of Tile profile given velocity ranges.
        """

        if noiseScale is None:
            noiseScale = 1.
        std = 0
        velMin = [lowrange[0],highrange[0]]
        velMax = [lowrange[1],highrange[1]]


        if self.profile is not None:
            for k in xrange(len(velMin)):
                noiseRegion = np.where((xarr >= velMin[k]) & \
                                (xarr <= velMax[k]))
                std = np.std(self.get_profile()[noiseRegion]) + std

            std = std / len(velMin)
            self.noise = std * noiseScale
        else:
            self.noise = np.NaN


################################################################################
# Functions
################################################################################

def compare_grids(grid0, grid1, velocityrange=[], plotallpixels=False,
        savedir='./', filename=None, show=True):

    ''' Plots heat map of N(HI) residuals between SpectralGrid instances
    between grid0 and grid1. The two SpectralGrids must be on the same
    coordinate grid.
    '''

    from kapteyn import wcs as kwcs
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs

    if len(velocityrange) != 2:
        sp = grid0.spectralAxis
        velocityrange = [min(sp),max(sp)]

    fig = plt.figure(figsize=(8,8))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=grid0.header)),
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    # Grid0
    image = np.empty((grid0.get_imagexsize(),grid0.get_imageysize()))
    image[:,:] = np.NaN
    region = grid0.region
    tilecount=0
    for i in xrange(region[0],region[2]):
        for j in xrange(region[1],region[3]):
            tile = grid0.get_tile(i,j,coords='pixel')
            if tile.ncomps > 0:
                for k in range(tile.ncomps):
                    vel = tile.params[k][1]
                    if vel < velocityrange[1] and vel > velocityrange[0]:
                        image[i,j] = tile.compAreaList[k]
            else:
                image[i,j] = np.NaN
            tilecount = tilecount + 1

    NHIimage0 = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2

    # Grid1
    image = np.empty((grid1.get_imagexsize(),grid1.get_imageysize()))
    image[:,:] = np.NaN
    region = grid1.region
    tilecount=0
    for i in xrange(region[0],region[2]):
        for j in xrange(region[1],region[3]):
            tile = grid1.get_tile(i,j,coords='pixel')
            if tile.ncomps > 0:
                for k in range(tile.ncomps):
                    vel = tile.params[k][1]
                    if vel < velocityrange[1] and vel > velocityrange[0]:
                        image[i,j] = tile.compAreaList[k]
            else:
                image[i,j] = np.NaN
            tilecount = tilecount + 1

    NHIimage1 = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2

    # Define plot properties
    ax = imagegrid[0]
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")
    ax.set_xlabel('Right Ascension (J2000)', 
              size = 'small', 
              family='serif')
    ax.set_ylabel('Declination (J2000)', 
              size = 'small',
              family='serif')
    if not plotallpixels:
        ax.set_xlim(grid0.region[0],grid0.region[2])
        ax.set_ylim(grid0.region[1],grid0.region[3])

    #cmap = plt.cm.get_cmap('rainbow',NHIimage.max()-NHIimage.min()+1)
    #cmap = ax.cm.get_cmap('gray')
    #cmap.set_bad('w',1.)
    #ax.imshow(NHIimage, interpolation='nearest',cmap=cmap,origin='lower')    #cbar = ax.colorbar()
    #cbar.set_label(r'N(H\,\textsc{i}) $\times$ 10$^{18}$ cm$^{-2}$')


    im = ax.imshow(NHIimage0 - NHIimage1,
                   interpolation='nearest',
                   origin='lower')
    cb = ax.cax.colorbar(im)
    # Write label to colorbar
    cb.set_label_text(r'N(H\,\textsc{i}) $\times$ 10$^{20}$ cm$^{-2}$',
                   size='small',
                   family='serif')
    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def display_image(grid, plotallpixels=False, savedir='./', filename=None,
                  show=True):
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs

    fig = plt.figure(figsize=(8,8))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=grid.header)),
                 aspect=True,
                 #label_mode='L',
                 share_all=True)

    ax = imagegrid[0]

    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',
              size = 'small',
              family='serif')
    ax.set_ylabel('Declination (J2000)',
              size = 'small',
              family='serif')

    #cmap = plt.cm.get_cmap('rainbow',NHIimage.max()-NHIimage.min()+1)
    #cmap = ax.cm.get_cmap('gray')
    #cmap.set_bad('w',1.)
    #ax.imshow(NHIimage, interpolation='nearest',cmap=cmap,origin='lower')    #cbar = ax.colorbar()
    #ax.rc('text', usetex=True)
    #ax.rc('font', family='serif')
    #cbar.set_label(r'N(H\,\textsc{i}) $\times$ 10$^{18}$ cm$^{-2}$')

    image = grid.cube.sum(axis=0)

    im = ax.imshow(image, interpolation='nearest',origin='lower')
    cb = ax.cax.colorbar(im)
    # Write label to colorbar
    cb.set_label_text(grid.header['bunit'],
                   size='small',
                   family='serif')
    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def load_grid(filename):
    """ Loads saved SpectralGrid instance.
    """

    from pickle import load

    with open(filename + '/SpectralGrid','rb') as f:
        SpectralGrid = tile=load(f)

    with open(filename + '/gridlist','rb') as f:
        SpectralGrid.gridlist = tile=load(f)

    return SpectralGrid

def make_velocityAxis(h):
    """ Creates the velocity axis given a header.
    """
    from numpy import arange

    array = (arange(h['NAXIS3']) - h['CRPIX3'] + 1) * h['CDELT3'] + h['CRVAL3']

    return array / 1000.

def make_radecAxes(h):
    """ Creates the RA and Dec axis given a header.
    """

    from numpy import arange

    ra = (arange(h['NAXIS1']) - h['CRPIX1'] + 1) * h['CDELT1'] + h['CRVAL1']
    dec = (arange(h['NAXIS2']) - h['CRPIX2'] + 1) * h['CDELT2'] + h['CRVAL2']

    return ra,dec

def plot_NHI(grid, velocityrange=[], plotallpixels=False, savedir='./',
             filename=None,show=True, returnimage=False):
    from kapteyn import wcs as kwcs
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs

    if len(velocityrange) != 2:
        sp = grid.spectralAxis
        velocityrange = [min(sp),max(sp)]

    fig = plt.figure(figsize=(8,8))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=grid.header)),
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    image = np.empty((grid.get_imagexsize(),grid.get_imageysize()))
    image[:,:] = np.NaN
    region = grid.region
    tilecount=0
    for i in xrange(region[0],region[2]):
        for j in xrange(region[1],region[3]):
            tile = grid.get_tile(i,j,coords='pixel')
            if tile.ncomps > 0:
                for k in range(tile.ncomps):
                    lowVel = tile.params[k][1] - tile.params[k][2]
                    highVel = tile.params[k][1] + tile.params[k][2]
                    if highVel < velocityrange[1] and lowVel > velocityrange[0]:
                        image[i,j] = 0.
                        image[i,j] += tile.compAreaList[k]
            else:
                image[i,j] = np.NaN
            tilecount = tilecount + 1

    #print tilecount

    NHIimage = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2

    #NHIimage[NHIimage == NHIimage.max()] = 1

    ax = imagegrid[0]
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)', 
              size = 'small', 
              family='serif')
    ax.set_ylabel('Declination (J2000)', 
              size = 'small',
              family='serif')
    if not plotallpixels:
        ax.set_ylim(grid.region[0],grid.region[2])
        ax.set_xlim(grid.region[1],grid.region[3])


    #cmap = plt.cm.get_cmap('rainbow',NHIimage.max()-NHIimage.min()+1)
    #cmap = ax.cm.get_cmap('gray')
    #cmap.set_bad('w',1.)
    #ax.imshow(NHIimage, interpolation='nearest',cmap=cmap,origin='lower')    #cbar = ax.colorbar()
    #cbar.set_label(r'N(H\,\textsc{i}) $\times$ 10$^{18}$ cm$^{-2}$')

    plt.rc('text', usetex=True)
    im = ax.imshow(NHIimage, interpolation='nearest',origin='lower')
    cb = ax.cax.colorbar(im)
    # Write label to colorbar
    cb.set_label_text(r'N(H\,\textsc{i}) $\times$ 10$^{20}$ cm$^{-2}$',
                   size='small',
                   family='serif')

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        fig.show()
    if returnimage:
        return NHIimage

def plot_fit(grid,x,y,coords=None,coRegion=True,savedir='./',
             filename=None,show=True):
    """ Plot the components of a gaussian fit.
    """

    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from sys import exit

    if coords is 'image':
        x,y = grid.image2grid(x,y)

    tile = grid.get_tile(x,y)
    sp = tile.profile
    noise = tile.noiseProfile
    #print noise
    xarr = grid.get_spectralAxis()
    residualXarr = tile.spectralAxis
    #print residuals.shape
    #print residualXarr.shape
    #print len(xarr)
    totGaussFit = np.zeros(len(xarr))
    ncomps = tile.get_ncomps()

    if ncomps is None:
        print '\n Profile has not been successfully fit for tile'
        exit()

    fig = plt.figure(figsize=[20,10])
    ax = fig.add_subplot(1,1,1)

    for i in xrange(0,tile.get_ncomps()):
        tileParams = grid.get_tile(x,y).get_compParams(i)
        print tileParams
        comp = grid.get_component(x,y,i)
        ax.plot(xarr,comp,'r--',markersize=1)
        totGaussFit = totGaussFit + comp

    #residuals = sp - totGaussFit
    residuals = sp - totGaussFit


    # Begin plotting
    if coRegion:
        if tile.coVelocities is not None:
            coVels = tile.coVelocities
            for i, coVel in enumerate(coVels):
                if i==0:
                    ax.fill_between(coVel,
                            -1e100,1e100,
                            facecolor='orange', alpha=0.5,
                            label='HI absorption region')
                else:
                    ax.fill_between(coVel,
                            -1e100,1e100,
                            facecolor='orange', alpha=0.5)
    ax.plot(xarr,residuals, 'g-',label='Residuals')
    ax.plot(residualXarr,residuals,'g-')
    ax.plot(xarr,totGaussFit,'y-.',label='Composite Model Spectrum')
    ax.plot(xarr,sp, 'b',label='Spectrum')
    ax.plot(0,0,'r--',label='Model Component')
    plt.legend(loc='upper right')
    ax.annotate('Pos: ' + str(tile.imageXPos) + ', ' \
            + str(tile.imageYPos) + '\n' + \
            '# components = ' + str(ncomps),
            xy=(0.02,0.90),xycoords='axes fraction')
    ax.set_ylim(2*residuals.min(), 1.05*sp.max())


    #from matplotlib.ticker import MultipleLocator
    #ml = MultipleLocator(5)
    #ax.xaxis.set_minor_locator(ml)
    #ax.xaxis.get_ticklines(minor=True)

    plt.setp(ax.get_xticklabels(), visible=False)
    plt.ylabel('T$_b$ (K)')

    # Plot residuals:
    divider = make_axes_locatable(ax)
    residAx = divider.append_axes("bottom", size=1.2, pad=0.1, sharex=ax)
    residAx.plot(xarr,residuals, 'g')
    residAx.plot(xarr,noise, 'k')
    residAx.plot(xarr,-noise,'k')
    if noise.max() > residuals.max(): maxy = noise.max()
    else: maxy = residuals.max()
    if noise.min() < residuals.min(): miny = noise.min()
    else: miny = residuals.min()
    residAx.set_ylim(1.1*miny,1.1*maxy)
    plt.setp(residAx.get_xticklabels(), visible=False)
    plt.ylabel('T$_b$ (K)')

    # Plot noise array:
    Tsys = 30. # Assume for Arecibo
    Tb = sp
    noise = (Tsys + Tb) / Tsys
    #divider2 = make_axes_locatable(residAx)
    noiseAx = divider.append_axes("bottom", size=1.2, pad=0.1, sharex=ax)
    noiseAx.plot(xarr,noise, 'k')
    noiseAx.set_ylim(noise.min() - 1,noise.max() + 1)
    plt.setp(noiseAx.get_xticklabels(), visible=True)
    plt.ylabel('$\eta_{eff}$')
    '''
    # Plot noise array - residuals:
    noiseScale = residuals * noise
    noiseScaleAx = divider.append_axes("bottom", size=1.2, pad=0.1, sharex=ax)
    noiseScaleAx.plot(xarr,noiseScale, 'k',label='Scaled Residuals')
    noiseScaleAx.plot(xarr,residuals, 'g',label='Raw Residuals')
    plt.legend(loc='upper right')
    noiseScaleAx.set_ylim(noiseScale.min() - 1,noiseScale.max() + 1)
    plt.ylabel('T$_b$ (K)')
    '''



    plt.xlim(xarr.min(),xarr.max())
    #plt.ylim(residuals.min() - 5, sp.max() + 5)

    plt.xlabel('Velocity (km/s)')

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def plot_fits(grid,xlist,ylist,coords=None,legend=False,coRegion=True,
              savedir='./',filename=None,show=True):
    """ Plot the components of a gaussian fit.
    """

    from matplotlib.ticker import MultipleLocator
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from sys import exit
    from mpl_toolkits.axes_grid1 import Grid

    fig = plt.figure(figsize=[20,10])
    plotGrid = Grid(fig, (1,1,1),
                 nrows_ncols=(1,len(xlist)),
                 ngrids=len(xlist),
                 direction='row',
                 axes_pad=0,
                 label_mode='L',
                 share_all=True)

    for i in range(len(xlist)):
        x = xlist[i]
        y = ylist[i]
        if coords is 'image':
            x,y = grid.image2grid(x,y)

        tile = grid.get_tile(x,y)
        sp = tile.profile
        noise = tile.noiseProfile
        xarr = grid.get_spectralAxis()
        totGaussFit = np.zeros(len(xarr))
        ncomps = tile.get_ncomps()

        if ncomps is None:
            print '\n Profile has not been successfully fit for tile'
            exit()

        ax = plotGrid[i]

        for i in xrange(0,tile.get_ncomps()):
            tileParams = grid.get_tile(x,y).get_compParams(i)
            print tileParams
            comp = grid.get_component(x,y,i)
            ax.plot(xarr,comp,'r',markersize=1)
            totGaussFit = totGaussFit + comp

        residuals = sp - totGaussFit

        if coRegion:
            if tile.coVelocities is not None:
                coVels = tile.coVelocities
                ax.fill_between(coVels,
                                -1e100,1e100,
                                facecolor='orange', alpha=0.5,
                                label='HI absorption region')
        ax.plot(xarr,residuals, 'g-',label='Residuals')
        ax.plot(xarr,totGaussFit,'y-.',label='Composite Model Spectrum')
        ax.plot(xarr,sp, 'b',label='Spectrum')
        ax.plot(0,0,'r',linestyle='solid',label='Model Component')

        if legend:
            plt.legend(loc='upper right')
        ax.annotate('Pos: ' + str(tile.imageXPos) + ', ' \
            + str(tile.imageYPos) + '\n' + \
            '# components = ' + str(ncomps),
            xy=(0.02,0.90),xycoords='axes fraction')
        ax.set_ylim(2*residuals.min(), 1.05*sp.max())


        plt.setp(ax.get_xticklabels(), visible=False)
        plt.xlim(xarr.min(),xarr.max())
    plt.ylabel('T$_b$ (K)')
    plt.xlabel('Velocity (km/s)')

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def plot_ncompImage(grid, plotallpixels=False, savedir='./', filename=None,
                    show=True):
    from kapteyn import wcs as kwcs
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs

    fig = plt.figure(figsize=(8,8))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=grid.header)),
                 aspect=False,
                 label_mode='L',
                 share_all=False)


    ncompImage = np.empty((grid.get_imagexsize(),grid.get_imageysize()))
    ncompImage[:,:] = np.NaN
    region = grid.region
    tilecount = 0

    for i in xrange(region[0],region[2]):
        for j in xrange(region[1],region[3]):
            tile = grid.get_tile(i,j,coords='pixel')
            if tile.get_ncomps() is not None:
                ncompImage[i,j] = tile.get_ncomps()
                tilecount = tilecount + 1
            else:
                ncompImage[i,j] = np.NaN

    image = np.ma.array(ncompImage,mask=np.isnan(ncompImage))

    '''
    cmap = plt.cm.get_cmap('rainbow',image.max()-image.min()+1)
    cmap.set_bad('w',1.)
    plt.imshow(image, interpolation='nearest',cmap=cmap,origin='lower')
    plt.colorbar()
    plt.show()
    '''

    ax = imagegrid[0]
    ax.set_display_coord_system("fk5")
    #ax.set_ticklabel_type("hms", "dms")


    ax.set_xlabel('Right Ascension (J2000)', 
              size = 'small', 
              family='serif')
    ax.set_ylabel('Declination (J2000)', 
              size = 'small',
              family='serif')
    if not plotallpixels:
        ax.set_ylim(grid.region[0],grid.region[2])
        ax.set_xlim(grid.region[1],grid.region[3])


    im = ax.imshow(image, interpolation='nearest',origin='lower')
    cb = ax.cax.colorbar(im)
    # Write label to colorbar
    cb.set_label_text(r'# of components',
                   size='small',
                   family='serif')

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def plot_componentPV(SpectralGrid,xslice=None,yslice=None,width=None,
                     savedir='./',filename=None,show=True):
    import matplotlib.pyplot as plt
    import numpy as np

    if width is None:
        width = 1

    colors = ['r','g','b','y','k','c']
    grid = SpectralGrid.gridlist

    if xslice is not None:
        print 'nothing'
    if yslice is not None:
        for j in range(0,width):
            for i in range(0,SpectralGrid.get_xsize()):
                tile = grid[i][yslice + j][0]
                if tile.get_ncomps() is not None:
                    vels = np.zeros(tile.get_ncomps())
                    areas = np.zeros(tile.get_ncomps())
                    for k in range(0,len(vels)):
                        vel = tile.params[k][1]
                        area = tile.compAreaList[k]
                        if vel > -100 and vel < 100:
                            vels[k] = vel
                            areas[k] = area

                    #velocities = tile.get_componentVelocities()
                    plt.scatter([i]*len(vels),vels,c=colors[j%len(colors)],
                            alpha=0.5,s=areas)

        plt.xlabel('Pixel')
        plt.ylabel('Velocity (km $s^{-1}$)')

        if filename is not None:
            plt.savefig(savedir + filename)
        if show:
            plt.show()

def plot_residualsImage(grid, plotallpixels=False, savedir='./', filename=None,
                        show=True):
    from kapteyn import wcs as kwcs
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs


    fig = plt.figure(figsize=(8,8))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=grid.header)),
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    residImage = np.zeros((grid.get_xsize(),grid.get_ysize()))

    for i in xrange(0,residImage.shape[0]):
        for j in xrange(0,residImage.shape[1]):
            try:
                residImage[i,j] = grid.get_tile(i,j).get_residuals().std()
            except AttributeError:
                residImage[i,j] = 0

    ax = imagegrid[0]
    ax.set_display_coord_system("fk5")
    #ax.set_ticklabel_type("hms", "dms")
    ax.set_xlabel('Right Ascension (J2000)', 
              size = 'small', 
              family='serif')
    ax.set_ylabel('Declination (J2000)', 
              size = 'small',
              family='serif')
    if not plotallpixels:
        ax.set_xlim(grid.region[0],grid.region[2])
        ax.set_ylim(grid.region[1],grid.region[3])

    im = ax.imshow(image, interpolation='nearest',origin='lower')
    cb = ax.cax.colorbar(im)
    # Write label to colorbar
    cb.set_label_text(r'# of components',
                   size='small',
                   family='serif')


    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def plot_spectrum(grid,x,y,coords=None,savedir='./',filename=None,show=True):
    """ Plot the components of a gaussian fit.
    """

    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from sys import exit

    if coords is 'image':
        x,y = grid.image2grid(x,y)

    tile = grid.get_tile(x,y)
    sp = tile.profile
    noise = tile.noiseProfile
    xarr = grid.get_spectralAxis()

    fig = plt.figure(figsize=[20,10])
    ax = fig.add_subplot(1,1,1)

    ax.plot(xarr,sp, 'b',label='Spectrum')
    ax.plot(0,0,'r--',label='Model Component')
    plt.legend(loc='upper right')
    ax.annotate('Pos: ' + str(tile.imageXPos) + ', ' + str(tile.imageYPos),
            xy=(0.02,0.90),xycoords='axes fraction')
    ax.set_ylim(0.95*sp.min(), 1.05*sp.max())


    #from matplotlib.ticker import MultipleLocator
    #ml = MultipleLocator(5)
    #ax.xaxis.set_minor_locator(ml)
    #ax.xaxis.get_ticklines(minor=True)

    plt.setp(ax.get_xticklabels(), visible=False)
    plt.ylabel('T$_b$ (K)')

    # Plot noise array:
    Tsys = 30. # Assume for Arecibo
    Tb = sp
    noise = (Tsys + Tb) / Tsys
    #divider2 = make_axes_locatable(residAx)
    divider = make_axes_locatable(ax)
    noiseAx = divider.append_axes("bottom", size=1.2, pad=0.1, sharex=ax)
    noiseAx.plot(xarr,noise, 'k')
    noiseAx.set_ylim(noise.min() - 1,noise.max() + 1)
    plt.setp(noiseAx.get_xticklabels(), visible=True)
    plt.ylabel('$\eta_{eff}$')
    '''
    # Plot noise array - residuals:
    noiseScale = residuals * noise
    noiseScaleAx = divider.append_axes("bottom", size=1.2, pad=0.1, sharex=ax)
    noiseScaleAx.plot(xarr,noiseScale, 'k',label='Scaled Residuals')
    noiseScaleAx.plot(xarr,residuals, 'g',label='Raw Residuals')
    plt.legend(loc='upper right')
    noiseScaleAx.set_ylim(noiseScale.min() - 1,noiseScale.max() + 1)
    plt.ylabel('T$_b$ (K)')
    '''



    plt.xlim(xarr.min(),xarr.max())
    #plt.ylim(residuals.min() - 5, sp.max() + 5)

    plt.xlabel('Velocity (km/s)')

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def plot_velocityHistograms(grid,savedir='./',filename=None,show=True):
    import matplotlib.pyplot as plt
    import numpy as np

    velocityList = []

    for i in xrange(0,grid.get_xsize()):
        for j in xrange(0,grid.get_ysize()):
            params = grid.get_tile(i,j).params
            ncomps = grid.get_tile(i,j).get_ncomps()
            if ncomps is not None:
                for k in xrange(0,ncomps):
                    if params[k][1] > -100 and params[k][1] < 100:
                        velocityList.append(params[k][1])

    velocityArray = np.asarray(velocityList)
    plt.xlim([velocityArray.min()-2,velocityArray.max()+2])
    plt.hist(velocityList,bins=100)
    plt.xlabel('Velocity (km/s)')
    plt.ylabel('Frequency')

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()

def write_grid(SpectralGrid,filename,verbose=True):
    """ Writes gridlist to file.
    """

    from pickle import dump
    from os import system,path

    save = True

    if path.isdir(filename):
        save = True
        if verbose:
            print ' Folder exists, SpectralGrid will be saved under as ' + \
                   filename
    elif path.isfile(filename):
        save = False
        if verbose:
            print ' File exists, SpectralGrid will not be saved.'
    else:
        system('mkdir ' + filename)

    if save:
        f = open(filename + '/gridlist','w')
        dump(SpectralGrid.gridlist,f,protocol=2)
        f.close()

        f = open(filename + '/SpectralGrid','w')
        dump(SpectralGrid,f,protocol=2)
        f.close()

def write_fits_NHI(SpectralGrid,filename,velocityrange=[]):
    ''' Calculates the column density of the SpectralGrid model profiles and
    writes as a fits file.

    Parameters
    ----------

    Returns
    -------
    '''

    # External modules
    import pyfits as pf

    # Use all velocities if not specified
    if len(velocityrange) != 2:
        sp = SpectralGrid.spectralAxis
        velocityrange = [min(sp),max(sp)]

    # Choose components
    image = np.empty((SpectralGrid.get_imagexsize(),
        SpectralGrid.get_imageysize()))
    image[:,:] = np.NaN
    region = SpectralGrid.region
    tilecount=0
    for i in xrange(region[0],region[2]):
        for j in xrange(region[1],region[3]):
            tile = SpectralGrid.get_tile(i,j,coords='pixel')
            if tile.ncomps > 0:
                for k in range(tile.ncomps):
                    lowVel = tile.params[k][1] - tile.params[k][2]
                    highVel = tile.params[k][1] + tile.params[k][2]
                    if highVel < velocityrange[1] and lowVel > velocityrange[0]:
                        image[i,j] = 0.
                        image[i,j] += tile.compAreaList[k]
            else:
                image[i,j] = np.NaN
            tilecount = tilecount + 1

    NHIimage = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2

    header = SpectralGrid.header
    header['NAXIS3'] = 1

    hdu = pf.PrimaryHDU(NHIimage,header=header)
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(filename)



