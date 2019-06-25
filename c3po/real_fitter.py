import os
import time
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy.io import ascii
import sed_config
import temps_config
import bound_options
from fitting_options import *


# Main loop definition.
def run_fits(starName):
    # Convenience variables
    wa, fx, er, kw, va = 'wavelength', 'flux', 'error', 'keywords', 'value'

    # Compositions of the grains.
    grainComps = ['AstroSil', 'DirtyIceAstroSil']

    # Instrument names
    mj, mh, mk      = '2MASSJ', '2MASSH', '2MASSK'
    w1, w2, w3, w4  = 'WISE1', 'WISE2', 'WISE3', 'WISE4'
    mp24, mp70      = 'MIPS24', 'MIPS70',
    mp24ul, mp70ul  = 'MIPS24UL', 'MIPS70UL'
    h70, h100, h160 = 'HerschelPACS70', 'HerschelPACS100', 'HerschelPACS160'
    h70ul, h100ul   = 'HerschelPACS70UL', 'HerschelPACS100UL'
    h160ul          = 'HerschelPACS160UL'
    sl2, sl1        = 'SpitzerIRS-SL2', 'SpitzerIRS-SL1'
    ll2, ll1        = 'SpitzerIRS-LL2', 'SpitzerIRS-LL1'

    # Grain compositions. This directly impacts the files that are loaded.
    innerGrain = 'AstroSil'
    outerGrain = 'DirtyIceAstroSil'

    # Grain densities. Important for calculating the blowout size of the grains.
    densities = {
        # For the density of the grains:
        # Mixture of H20:NH3 is 3:1 and 10% vlfr of aC to make dirty ice.
        # then dirty ice is 50% vlfr of AstroSil to make DIAS
        # Assuming:
        # 1.0 g/cm^3 for H20
        # 0.817 g/cm^3 for NH3
        # 3.0 g/cm^3 for AstroSil
        'DirtyIce': 1.07,
        'AstroSil': 2.7,
        'DirtyIceAstroSil': 1.885, # 1:1 DirtyIce:AstroSil
        'WaterIce': 1
        # 'AstroSil': 3.0,           # Correlates to DIAS -> 2.034 g/cm^3
        # 'DirtyIceAstroSil': 2.034, # Correlates to AS -> 3.0 g/cm^3
        }

    # Read in the star file
    starData = ascii.read(sed_config.STAR_FILES+'{}_stitched.txt'.format(starName))

    # Read in the stellar properties
    starT    = starData.meta[kw]['TEMP'][va]    # Stellar Temperature
    starL    = starData.meta[kw]['starL'][va]   # Stellar Luminosity
    starM    = starData.meta[kw]['starM'][va]   # Stellar Mass
    specType = starData.meta[kw]['SpType'][va]  # Spectral Type
    starD    = starData.meta[kw]['DIST_pc'][va] # Distance to Earth

    # Fillers for in case the data is missing.
    starL = 1 if np.isnan(starL) else starL
    starM = 1 if np.isnan(starM) else starM
    starD = 1 if np.isnan(starD) else starD

    # Concerning the nextgen and kurucz models: the wavelengths are in units of
    # angstrom, which equals 10^-10 meters. Flux is in erg/cm^2/s/a.
    # Create temp array that matches the temperatures of the nextgen stellar
    # models. From 2600-4100, uses a step of 100. From 4200-10000, uses a step
    # of 200. Grabs the nearest model to the star temp in the star file.
    ngTemps    = np.arange(2600, 4100, 100)
    ngTemps    = np.append(ngTemps, np.arange(4200, 10200, 200))
    TEMP       = sed_config.find_nearest(ngTemps, starT)
    ngfilename = 'xp00_'+str(TEMP)+'g40.txt'

    # Temperature label for the graph of the SED fit
    starLabel = TEMP

    # Load the nextgen stellar model
    ngWave, ngFlux   = np.loadtxt(sed_config.NEXTGEN+ngfilename, unpack=True)
    ngWave, ngwIndex = np.unique(ngWave, return_index=True)
    ngFlux           = ngFlux[ngwIndex]

    # Convert the nextgen model to janskies
    c_cgs  = 2.99792458e10                          #cm/s
    ngFnu  = ngFlux*(ngWave**2)/c_cgs/1.0e8*1.0e23  #Jy
    ngWave = ngWave*1.0e-4                          #cm -> um

    if TEMP > 10000:
        # Load kurucz model
        kzTemps = np.arange(3500, 10000, 250)
        kzTemps = np.append(kzTemps, np.arange(10000, 13000, 500))
        kzTemps = np.append(kzTemps, np.arange(13000, 35000, 1000))
        kzTemps = np.append(kzTemps, np.arange(35000, 52500, 2500))
        kzTemp = find_nearest(kzTemps, starT)
        starLabel = kzTemp

        kzfilename = 'kp00_'+str(kzTemp)+'g40.txt'
        kzWave, kzFlux = np.loadtxt(sed_config.KURUCZ+kzfilename, unpack=True)
        kzWave, kzwIndex = np.unique(kzWave, return_index=True)
        kzFlux = kzFlux[kzwIndex]

        # Convert the kurucz model to janskies
        c_cgs = 2.99792458e10                          #cm/s
        kzFnu = kzFlux*(kzWave**2)/c_cgs/1.0e8*1.0e23  #Jy
        kzWave = kzWave*1.0e-4                         #cm -> um
        index = np.where(kzWave < 2.154)
        kzFnu = kzFnu[index]
        kzWave = kzWave[index]

        index = np.where(ngWave >= 2.154)
        ngFnu = ngFnu[index]
        ngWave = ngWave[index]

        s_norm = np.nanmean(kzFnu[-5:])/np.nanmean(ngFnu[:5])
        ngFnu *= s_norm

        ngWave = np.array(list(kzWave)+list(ngWave))
        ngFnu = np.array(list(kzFnu)+list(ngFnu))

    # Log all instruments and sort data by instrument. sData is a dictionary,
    # the keys of which are the instrument names. Within that are the fluxes,
    # the wavelengths, and errors.
    insts, sData = sed_config.sort_by_instrument(starData)

    # Boolean flags used to simplify accounting for which modules are present
    # in the data
    mjf, mhf, mkf    = mj in insts, mh in insts, mk in insts
    w1f, w2f         = w1 in insts, w2 in insts
    w3f, w4f         = w3 in insts, w4 in insts
    sl2f, sl1f       = sl2 in insts, sl1 in insts
    ll2f, ll1f       = ll2 in insts, ll1 in insts
    mp24f, mp24ulf   = mp24 in insts, mp24ul in insts
    mp70f, mp70ulf   = mp70 in insts, mp70ul in insts
    h70f, h100f      = h70 in insts, h100 in insts
    h160f, h70ulf    = h160 in insts, h70ul in insts
    h100ulf, h160ulf = h100ul in insts, h160ul in insts

    # Iterator lists for the flags
    # NOTE: must be in matching order
    spitzFlags = [sl2f, sl1f, ll2f, ll1f]
    spitzNames = [sl2,  sl1,  ll2,  ll1 ]
    instFlags  = [mjf,mhf,mkf,w1f,w2f,w3f,w4f,mp24f,mp70f,h70f,h100f,h160f]
    instNames  = [mj, mh, mk, w1, w2, w3, w4, mp24, mp70, h70, h100, h160 ]
    ulFlags    = [mp24ulf, mp70ulf, h70ulf, h100ulf, h160ulf]
    ulNames    = [mp24ul,  mp70ul,  h70ul,  h100ul,  h160ul ]

    # Saturation limits for each instrument
    satLims = {mj:10.057, mh:10.24, mk:10.566, w1:0.18, w2:0.36, w3:0.88,
        w4:12.0, h70: 220., h100:510., h160:1125., mp24:np.inf, mp70:np.inf}
    # Create a list of instrument names for non-spitzer/non-upper limit data
    totalInsts = list()
    for i, f in enumerate(instFlags):
        if f:
            totalInsts.append(instNames[i])
    # Create a list of instrument names for upper limits data
    ulInsts = list()
    for i, f in enumerate(ulFlags):
        if f:
            ulInsts.append(ulNames[i])
    # Create a list of instruments names for spitzer data
    spitzInsts = list()
    for i, f in enumerate(spitzFlags):
        if f:
            spitzInsts.append(spitzNames[i])
    # Colors for plotting accessed by instrument name
    plotColors = {
        mj: 'r', mh: 'r', mk: 'r',
        w1: 'b', w2: 'b', w3: 'b', w4: 'b',
        mp24:   'g', mp70:   'g',
        mp24ul: 'g', mp70ul: 'g',
        h70:   'purple', h100:   'purple', h160:   'purple',
        h70ul: 'purple', h100ul: 'purple', h160ul: 'purple',
        }
    # Begin the arrays that will be used in the fitting optimization function.
    fitWaves = np.array([])
    fitFlux = np.array([])
    fitError = np.array([])
    if totalInsts:
        for inst in totalInsts:
            if saturation_limits:
                ind = np.where(sData[inst][fx] < satLims[inst])
                fitWaves = np.append(fitWaves, sData[inst][wa][ind])
                fitFlux  = np.append(fitFlux, sData[inst][fx][ind])
                fitError = np.append(fitError, sData[inst][er][ind])
            else:
                fitWaves = np.append(fitWaves, sData[inst][wa])
                fitFlux  = np.append(fitFlux, sData[inst][fx])
                fitError = np.append(fitError, sData[inst][er])
    if ulInsts:
        for inst in ulInsts:
            fitWaves = np.append(fitWaves, sData[inst][wa])
            fitFlux  = np.append(fitFlux, sData[inst][fx])
            fitError = np.append(fitError, sData[inst][er])

    fitters = np.where((fitWaves/fitFlux) < 3)
    fitWaves = fitWaves[fitters]
    fitFlux = fitFlux[fitters]
    fitError = fitError[fitters]

    ind = np.isfinite(fitFlux)
    fitWaves = fitWaves[ind]
    fitFlux = fitFlux[ind]
    fitError = fitError[ind]

    ind = np.isfinite(fitError)
    fitWaves = fitWaves[ind]
    fitFlux = fitFlux[ind]
    fitError = fitError[ind]

    # Create arrays for the stitched data
    if spitzInsts:
        spitzWaves = np.array([])
        spitzFlux  = np.array([])
        spitzError = np.array([])
        for ins in spitzInsts:
            spitzWaves = np.append(spitzWaves, sData[ins][wa])
            spitzFlux  = np.append(spitzFlux, sData[ins][fx])
            spitzError = np.append(spitzError, sData[ins][er])

    # Do convolving if there's MIPS24 data and IRS
    if mp24f and len(spitzInsts) > 0:
        # Convolve IRS data to the MIPS24 data
        MIPS24W = sData[mp24][wa]
        MIPS24F = sData[mp24][fx]
        mipsw, mipsr = np.loadtxt(sed_config.FILTERS_DIR + 'mips24_frf.txt',
            unpack=True)
        IRS24      = sed_config.convolve(mipsw, mipsr, spitzWaves, spitzFlux)
        spitzFlux *= (MIPS24F/IRS24)
        fitWaves   = np.append(fitWaves, spitzWaves)
        fitFlux    = np.append(fitFlux, spitzFlux)
        fitError   = np.append(fitError, spitzError)
    # Else, just add the spitzer data
    elif not mp24f and len(spitzInsts) > 0:
        fitWaves  = np.append(fitWaves, spitzWaves)
        fitFlux   = np.append(fitFlux, spitzFlux)
        fitError  = np.append(fitError, spitzError)

    # Organize all data by increasing wavelength.
    ind = np.argsort(fitWaves)
    fitWaves = fitWaves[ind]
    fitFlux  = fitFlux[ind]
    fitError = fitError[ind]

    # Normalize stellar model from either SL2 or 2MASSK data.
    if sl2f:
        ind          = np.where(np.logical_and(fitWaves>5, fitWaves<5.5))
        dataFluxNorm = np.nanmean(fitFlux[ind])
        ngFluxNorm   = np.nanmean(np.interp(fitWaves[ind], ngWave, ngFnu))
        n_3   = (dataFluxNorm/ngFluxNorm)
        ngFnu = n_3 * ngFnu
    elif mkf:
        dataFluxNorm = sData[mk][fx]
        ngFluxNorm   = np.nanmean(np.interp(sData[mk][wa], ngWave, ngFnu))
        n_3   = (dataFluxNorm/ngFluxNorm)
        ngFnu = n_3 * ngFnu

    # Interpolate the stellar model to the fitting data wavelengths
    ngFnu_fit = np.e**np.interp(np.log(fitWaves), np.log(ngWave), np.log(ngFnu))

    # Create the star objects
    # Calculate blowoutsize given grain density
    blowoutSize1 = sed_config.blowout_size(densities[innerGrain], starL, starM)*bosScalar1
    blowoutSize2 = sed_config.blowout_size(densities[outerGrain], starL, starM)*bosScalar2

    # Index of grains greater than the blowout size
    graindex1 = sed_config.find_nearest_ind(sed_config.GRAINSIZES, blowoutSize1)
    graindex2 = sed_config.find_nearest_ind(sed_config.GRAINSIZES, blowoutSize2)

    # Load emissivities per grain comp (use these for hi res plots)
    TOTAL_EMISSIVITIES = dict()
    for grain in grainComps:
        TOTAL_EMISSIVITIES[grain] = np.load(sed_config.ARR_DIR+grain+'Emissivities.npy')

    # NOTE: this part is possibly inefficient? can fix later
    # hi res arrays needed to calculate the grain temps
    star2 = sed_config.Star(starD, starL, starT, None, blowoutSize1,
        blowoutSize2, TOTAL_EMISSIVITIES, sed_config.GRAINSIZES, graindex1,
        graindex2)

    # Load the grain temperatures per grain composition
    grainTemps = dict()
    # These files should be made already, but if not, then the script will
    # create them and then save them for subsequent fits.
    for grain in grainComps:
        try:
            grainTemps[grain] = np.load(
                sed_config.GRAIN_TEMPS_DIR+"%s_%s.npy"%(starName, grain))
        except:
            grainTemps[grain] = temps_config.calcTemps(star2, grain)
            np.save(sed_config.GRAIN_TEMPS_DIR+"%s_%s.npy"%(starName, grain),
                grainTemps[grain], allow_pickle=False)

    # hi res arrays, now with grain temps
    star2 = sed_config.Star(starD, starL, starT, grainTemps, blowoutSize1,
        blowoutSize2, TOTAL_EMISSIVITIES, sed_config.GRAINSIZES, graindex1,
        graindex2)

    # Grab the minimum radial location that's below the temp
    # of sublimation for volatiles. (Icy belt radius)
    # This could potentially break in the future if the grain temps happen to be
    # broken. In that case, we would NOT want a for/else loop because it would
    # be ideal to catch a case of incorrect temps instead
    for r in range(1000):
        # If mean temp < 120, then that becomes the min radius
        if np.nanmean(grainTemps[outerGrain][r]) < 120:
            minRad = r
            break

    minRad = np.logspace(-1, 3, 1000)[minRad]

    print( '----------------------------------------' )
    print( '      AS blowout size: %.2f'  % blowoutSize1 )
    print( '     IMP blowout size: %.2f' % blowoutSize2 )
    print( '       r0 lower limit: %.2f'  % (0.5*starL) )
    print( ' Minimum radius for an icy belt: %.2f' % (minRad) )

    # Interp emissivities to fitWaves (use these for fitting)
    emis = {}
    emis[innerGrain] = np.empty((sed_config.GRAINSIZES.size, fitWaves.size))
    for g in range(sed_config.GRAINSIZES.size):
        emis[innerGrain][g] = np.interp(fitWaves, sed_config.WAVELENGTHS,
            TOTAL_EMISSIVITIES[innerGrain][g])
    emis[outerGrain] = np.empty((sed_config.GRAINSIZES.size, fitWaves.size))
    for g in range(sed_config.GRAINSIZES.size):
        emis[outerGrain][g] = np.interp(fitWaves, sed_config.WAVELENGTHS,
            TOTAL_EMISSIVITIES[outerGrain][g])

    # Warm belt guess is half the distance between the cold belt min & the star
    bbr1 = minRad*0.5

    # Cold belt guess is either from spatially resolved radii or 1.5*minRad
    if useSpatialRadii:
        for i in range(len(spatialRadii)):
            if starName == spatialRadii[i][0] or starName == spatialRadii[i][1]:
                bbr2 = np.float(spatialRadii[i][2])
                bbr2_unc = np.float(spatialRadii[i][3])
                try:
                    if 3*bbr2_unc > bbr2 or np.isnan(bbr2):
                        raise ValueError("No valid radius from Herschel")
                    sTrigger = 1
                except:
                    bbr2 = minRad*1.5
                    sTrigger = 0
                break
        else:
            bbr2 = minRad*1.5
            sTrigger = 0
    else:
        bbr2 = minRad*1.5
        sTrigger = 0

    if bbr2 > 1000:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(f"The initial guess for the cold belt is {bbr2} which is greater")
        print("than the highest possible value of 1000. Continuing with next")
        print("star. In future, we may expand radial locations?")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return

    print( f"initial guess for RW: {bbr1: 0.2f}" )
    print( f"initial guess for RC: {bbr2: 0.2f}" )

    # Low res for fitting star
    star1 = sed_config.Star(starD, starL, starT, grainTemps, blowoutSize1,
        blowoutSize2, emis, sed_config.GRAINSIZES, graindex1, graindex2)
    # Normalize warm and cold dust
    if oneWander:
        bb1 = star1.calcFluxBoSWarm(fitWaves, bbr1, blowoutSize1)
        bb2 = star1.calcFluxCold(fitWaves, bbr2)
    elif twoWander:
        bb1 = star1.calcFluxBoSWarm(fitWaves, bbr1, blowoutSize1)
        bb2 = star1.calcFluxBoSCold(fitWaves, bbr2, blowoutSize2)
    elif noWander:
        bb1 = star1.calcFluxWarm(fitWaves, bbr1)
        bb2 = star1.calcFluxCold(fitWaves, bbr2)
    elif twoWarmBelts:
        bb1 = star1.calcFluxWarm(fitWaves, bbr1)
        bb2 = star1.calcFluxCold(fitWaves, bbr2)
    # warm
    n_1 = fitFlux[np.where(fitFlux>0)].min() / bb1.max() # new test
    # cold
    n_2 = fitFlux[-3:].max() / bb2.max()
    # Reset norm factor for ngfNu
    n_3 = 1

    if showNormedBelts:
        # Show the normalized belts with observed flux and stellar model
        fig = plt.figure(figsize=(8, 6))
        plt.plot(fitWaves, n_1*bb1, label='Warm Belt')
        plt.plot(fitWaves, n_2*bb2, label='Cold Belt')
        plt.scatter(fitWaves, fitFlux, label='Observed Data')
        plt.plot(ngWave, ngFnu, label='Stellar Model')
        plt.xlim(1, 300)
        plt.ylim(n_1*bb1.max()*0.1, n_2*bb2.max()*100 )
        plt.title(starName+': for seeing the normalized belts', fontsize=22)
        plt.legend()
        plt.xlabel(r'$\lambda$ ($\mu m$)')
        plt.ylabel(r'$F_{\nu}$ ($Jy$)')
        plt.semilogx()
        plt.semilogy()
        plt.show()
        return


    ############################################################################
    ####      Fitting functions for optimization routine
    ############################################################################
    # option: noWander
    def warmCold(waves, r0warm, r0cold, n1, n2, n3):
        ''' Standard realistic grain fit: one warm belt, one cold belt '''
        return n1*star1.calcFluxWarm(waves, r0warm) \
            + n2*star1.calcFluxCold(waves, r0cold) \
            + n3*ngFnu_fit
    # option: oneWander
    def oneWarmWander(waves, r0warm, r0cold, bos, n1, n2, n3):
        ''' one warm belt with a wandering blowout grain size, one cold belt '''
        return n1*star1.calcFluxBoSWarm(waves, r0warm, bos) \
            + n2*star1.calcFluxCold(waves, r0cold) \
            + n3*ngFnu_fit
    # option: twoWander
    def warmColdWander(waves, r0warm, r0cold, bos1, bos2, n1, n2, n3):
        ''' one warm belt with a wandering blowout grain size,
        one cold belt with a wandering blowout grain size'''
        return n1*star1.calcFluxBoSWarm(waves, r0warm, bos1) \
            + n2*star1.calcFluxBoSCold(waves, r0cold, bos2) \
            + n3*ngFnu_fit
    # option: twoWarmBelts
    def twoWarmCold(waves, r0warm, r0cold, bos1, n1, n2, n3, n4):
        return n1*star1.calcFluxWarm(waves, r0warm) \
            + n2*star1.calcFluxCold(waves, r0cold) \
            + n3*ngFnu_fit \
            + n4*star1.calcFluxWarmMinGrains(waves, r0warm, bos1)

    ############################################################################
    ####             Run the optimization routine
    ############################################################################
    print("----------------------------------------")
    print( "      Begin Optimizing Parameters" )
    print("----------------------------------------")
    rw = bbr1
    rc = bbr2
    # Set initial guess for the blowout size. If the bounded region includes the
    # calculated blowout size, use that, else go in between the bounds.
    if  bosHighInner < 1 or 1./bosLowInner > 1:
        p0bos1 = blowoutSize1*np.average([bosHighInner, 1./bosLowInner])
    else:
        p0bos1 = blowoutSize1
    if  bosHighOuter < 1 or 1./bosLowOuter > 1:
        p0bos2 = blowoutSize2*np.average([bosHighOuter, 1./bosLowOuter])
    else:
        p0bos2 = blowoutSize2
    if oneWander:
        # parameters: r0warm, r0cold, bos, n1, n2, n3
        if useSpatialRadii and sTrigger:
            lBounds = [0.3, bbr2-bbr2_unc, 1e-3,
                    n_1/beltBound, n_2/beltBound, n_3*0.8]
            uBounds = [minRad, bbr2+bbr2_unc, bosHighInner*blowoutSize1,
                    n_1*beltBound, n_2*beltBound, n_3*1.2]
        else:
            lBounds = [0.3, minRad, 1e-3, n_1/beltBound, n_2/beltBound, n_3*0.8]
            uBounds = [minRad, maxRad, bosHighInner*blowoutSize1,
                n_1*beltBound, n_2*beltBound, n_3*1.2]
        p0=(rw, rc, p0bos1, n_1, n_2, n_3)
    elif twoWander:
        # parameters: r0warm, r0cold, bos1, bos2, n1, n2, n3
        if useSpatialRadii and sTrigger:
            lBounds = [0.3, bbr2 - bbr2_unc, blowoutSize1/bosLowInner, blowoutSize2/bosLowOuter,
                    n_1/beltBound, n_2/beltBound, n_3*0.8]
            uBounds = [minRad, bbr2 + bbr2_unc, blowoutSize1*bosHighInner, blowoutSize2*bosHighOuter,
                    n_1*beltBound, n_2*beltBound, n_3*1.2]
        else:
            lBounds = [0.3, minRad, blowoutSize1/bosLowInner, blowoutSize2/bosLowOuter,
                    n_1/beltBound, n_2/beltBound, n_3*0.8]
            uBounds = [minRad, maxRad, blowoutSize1*bosHighInner, blowoutSize2*bosHighOuter,
                    n_1*beltBound, n_2*beltBound, n_3*1.2]
        p0=(rw, rc, p0bos1, p0bos2, n_1, n_2, n_3)
    elif noWander:
        # parameters: r0warm, r0cold, n1, n2, n3
        if useSpatialRadii and sTrigger:
            lBounds = [0.3, bbr2-bbr2_unc, n_1/beltBound, n_2/beltBound, n_3*0.8]
            uBounds = [minRad, bbr2+bbr2_unc, n_1*beltBound, n_2*beltBound, n_3*1.2]
        else:
            lBounds = [0.3, minRad, n_1/beltBound, n_2/beltBound, n_3*0.8]
            uBounds = [minRad, maxRad, n_1*beltBound, n_2*beltBound, n_3*1.2]
        p0=(rw, rc, n_1, n_2, n_3)
    elif twoWarmBelts:
        # params: r0warm, r0cold, bos1, n1, n2, n3, n4
        if useSpatialRadii and sTrigger:
            lBounds = [0.3, bbr2-bbr2_unc, 0.001, n_1/beltBound, n_2/beltBound, n_3*0.8, n_1/beltBound]
            uBounds = [minRad, bbr2+bbr2_unc, blowoutSize1, n_1*beltBound, n_2*beltBound, n_3*1.2, n_1*beltBound]
        else:
            lBounds = [0.3, minRad, 0.001, n_1/beltBound, n_2/beltBound, n_3*0.8, n_1/beltBound]
            uBounds = [minRad, 500, blowoutSize1, n_1*beltBound, n_2*beltBound, n_3*1.2, n_1*beltBound]
        p0=(rw, rc, blowoutSize1/2, n_1, n_2, n_3, n_1)
    bounds =[lBounds, uBounds]

    before = time.perf_counter() # Timer for the routine

    # Call the optimization routine here.
    try:
        if oneWander:
            popt, pcov = curve_fit(
                oneWarmWander, fitWaves, fitFlux, sigma=fitError,
                p0=p0,
                # absolute_sigma=True,
                bounds=bounds
                )
        elif twoWander:
            popt, pcov = curve_fit(
                warmColdWander, fitWaves, fitFlux, sigma=fitError,
                p0=p0,
                # absolute_sigma=True,
                bounds=bounds
                )
        elif noWander:
            popt, pcov = curve_fit(
                warmCold, fitWaves, fitFlux, sigma=fitError,
                p0=p0,
                # absolute_sigma=True,
                bounds=bounds
                )
        elif twoWarmBelts:
            popt, pcov = curve_fit(
                    twoWarmCold, fitWaves, fitFlux, sigma=fitError,
                    p0=p0,
                    bounds=bounds
                    )
    except RuntimeError:
        print("RuntimeError: The least squares minimization failed.")
        print("Continuing with next star")
        return
    except ValueError:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("ValueError: Either xdata or ydata contains NaN")
        print("Check flux sigma for nan values.")
        print("Use those points as UL data instead.")
        print("Continuing with next star")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return
    except OptimizeWarning:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("The covariance of the parameters cannot be estimated")
        print("Continuing with next star")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return


    # Unpack the parameters and print them to the terminal for the user
    if oneWander:
        RW, RC, bos, n1, n2, n3 = popt
        with open(PARAMS_DIR + f"{starName}.txt", 'w+') as f:
            f.write(f"Warm radius: {RW}\n"
                    f"Cold radius: {RC}\n"
                    f"Warm Blowoutsize: {blowoutSize1}\n"
                    f"Cold Blowoutsize: {blowoutSize2}\n"
                    f"Warm Blowoutsize Scalar Multiple: {bosScalar1}\n"
                    f"Cold Blowoutsize Scalar Multiple: {bosScalar2}\n"
                    f"n1: {n1}\n"
                    f"n2: {n2}\n"
                    f"n3: {n3}\n")
    elif twoWander:
        RW, RC, bos1, bos2, n1, n2, n3 = popt
        with open(PARAMS_DIR + f"{starName}.txt", 'w+') as f:
            f.write(f"Warm radius: {RW}\n"
                    f"Cold radius: {RC}\n"
                    f"Warm Blowoutsize: {blowoutSize1}\n"
                    f"Cold Blowoutsize: {blowoutSize2}\n"
                    f"Warm Blowoutsize Scalar Multiple: {bosScalar1}\n"
                    f"Cold Blowoutsize Scalar Multiple: {bosScalar2}\n"
                    f"n1: {n1}\n"
                    f"n2: {n2}\n"
                    f"n3: {n3}\n")
    elif noWander:
        RW, RC, n1, n2, n3 = popt
        with open(PARAMS_DIR + f"{starName}.txt", 'w+') as f:
            f.write(f"Warm radius: {RW}\n"
                    f"Cold radius: {RC}\n"
                    f"Warm Blowoutsize: {blowoutSize1}\n"
                    f"Cold Blowoutsize: {blowoutSize2}\n"
                    f"Warm Blowoutsize Scalar Multiple: {bosScalar1}\n"
                    f"Cold Blowoutsize Scalar Multiple: {bosScalar2}\n"
                    f"n1: {n1}\n"
                    f"n2: {n2}\n"
                    f"n3: {n3}\n")
    elif twoWarmBelts:
        RW, RC, bos1, n1, n2, n3, n4 = popt
        with open(PARAMS_DIR + f"{starName}.txt", 'w+') as f:
            f.write(f"Warm radius: {RW}\n"
                    f"Cold radius: {RC}\n"
                    f"Warm Blowoutsize: {blowoutSize1}\n"
                    f"Cold Blowoutsize: {blowoutSize2}\n"
                    f"Warm Blowoutsize Scalar Multiple: {bosScalar1}\n"
                    f"Cold Blowoutsize Scalar Multiple: {bosScalar2}\n"
                    f"n1: {n1}\n"
                    f"n2: {n2}\n"
                    f"n3: {n3}\n"
                    f"n4: {n4}\n")
    print( '  Time to optimize parameters: %.2fs' % (time.perf_counter() - before) )
    print( '        Warm radius: %.2f' % RW )
    print( '        Cold radius: %.2f' % RC )
    print( '          Warm norm: %.16f' % n1 )
    print( '          Cold norm: %.16f' % n2 )
    print( '       Stellar norm: %.2f' % n3 )
    if oneWander:
        print( '      BlowoutSize: %.2f' % bos )
        print( f"Warm Blowoutsize Scalar Multiple: {bosScalar1}\n"
               f"Cold Blowoutsize Scalar Multiple: {bosScalar2}\n")
    elif twoWander:
        print( f'  Warm BlowoutSize: {bos1:.2f}')
        print( f'  Cold BlowoutSize: {bos2:.2f}')
        print(f"Warm Blowoutsize Scalar Multiple: {bosScalar1}\n"
              f"Cold Blowoutsize Scalar Multiple: {bosScalar2}\n")
    elif noWander:
        pass
    elif twoWarmBelts:
        print( f'  Minimum BlowoutSize: {bos1:.2f}')
    print( '----------------------------------------' )

    # Calculate chi square value
    if oneWander:
        resid = (fitFlux-oneWarmWander(fitWaves, RW, RC, bos, n1, n2, n3))/fitError
        degsFreedom = fitWaves.size - 6
    elif twoWander:
        resid = (fitFlux-warmColdWander(fitWaves, RW, RC, bos1, bos2, n1, n2, n3))/fitError
        degsFreedom = fitWaves.size - 7
    elif noWander:
        resid = (fitFlux-warmCold(fitWaves, RW, RC, n1, n2, n3))/fitError
        degsFreedom = fitWaves.size - 5
    elif twoWarmBelts:
        resid = (fitFlux-twoWarmCold(fitWaves, RW, RC, bos1, n1, n2, n3, n4))/fitError
        degsFreedom = fitWaves.size - 7
    chiSqr = np.dot(resid, resid)/degsFreedom # Reduced chi square

    SMALL_SIZE = 8
    MEDIUM_SIZE = 12
    BIGGER_SIZE = 18
    if showIRSVariance:
        # Create a figure with two plots, one for the model + data, one for the
        # IRS variance
        fig = plt.figure(figsize=(8, 10))
        ax = fig.add_subplot(211)
    else:
        # Otherwise just create a figure with a single plot for the model + data
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
    plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)   # fontsize of the figure title

    # Calculate the high resolution totalFlux array given the optimized
    # parameters
    if oneWander:
        y1 = n1 * star2.calcFluxBoSWarm(sed_config.WAVELENGTHS, RW, bos)
        y2 = n2 * star2.calcFluxCold(sed_config.WAVELENGTHS, RC)
        totalFlux = y1 + y2 + np.e**np.interp(np.log(sed_config.WAVELENGTHS), np.log(ngWave),
            np.log(n3*ngFnu))
    elif twoWander:
        y1 = n1 * star2.calcFluxBoSWarm(sed_config.WAVELENGTHS, RW, bos1)
        y2 = n2 * star2.calcFluxBoSCold(sed_config.WAVELENGTHS, RC, bos2)
        totalFlux = y1 + y2 + np.e**np.interp(np.log(sed_config.WAVELENGTHS), np.log(ngWave),
            np.log(n3*ngFnu))
    elif noWander:
        y1 = n1 * star2.calcFluxWarm(sed_config.WAVELENGTHS, RW)
        y2 = n2 * star2.calcFluxCold(sed_config.WAVELENGTHS, RC)
        totalFlux = y1 + y2 + np.e**np.interp(np.log(sed_config.WAVELENGTHS), np.log(ngWave),
            np.log(n3*ngFnu))
    elif twoWarmBelts:
        y1 = n1 * star2.calcFluxWarm(sed_config.WAVELENGTHS, RW)
        y2 = n2 * star2.calcFluxCold(sed_config.WAVELENGTHS, RC)
        y3 = n4 * star2.calcFluxBoSWarm(sed_config.WAVELENGTHS, RW, bos1)
        totalFlux = y1 + y2 + y3 + np.e**np.interp(np.log(sed_config.WAVELENGTHS), np.log(ngWave),
            np.log(n3*ngFnu))

    # Plot realistic grain fluxes
    if oneWander:
        plt.plot(sed_config.WAVELENGTHS, y1, ls='--', color='blue',
            label=r'R$_0$ Warm: %.2f AU'%RW)
        plt.plot(sed_config.WAVELENGTHS, y2, ls='--', color='r',
            label='R$_0$ Cold: %.2f AU'%RC)
    elif twoWander:
        plt.plot(sed_config.WAVELENGTHS, y1, ls='--', color='blue',
            label=r'R$_0$ Warm: %.2f AU'%RW)
        plt.plot(sed_config.WAVELENGTHS, y2, ls='--', color='r',
            label='R$_0$ Cold: %.2f AU'%RC)
    elif noWander:
        plt.plot(sed_config.WAVELENGTHS, y1, ls='--', color='blue',
            label=r'R$_0$ Warm: %.2f AU'%RW)
        plt.plot(sed_config.WAVELENGTHS, y2, ls='--', color='r',
            label='R$_0$ Cold: %.2f AU'%RC)
    elif twoWarmBelts:
        plt.plot(sed_config.WAVELENGTHS, y1, ls='--', color='blue',
            label=r'R$_0$ Warm: %.2f AU'%RW)
        plt.plot(sed_config.WAVELENGTHS, y2, ls='--', color='r',
            label=r'R$_0$ Cold: %.2f AU'%RC)
        plt.plot(sed_config.WAVELENGTHS, y3, ls='-.', color='blue',
            label=r'R$_0$ SmallGrains: %.2f AU'%RW)

    # Plot stellar model, total flux, and IRS data
    plt.plot(ngWave, n3*ngFnu, color = 'gray',
        label='Next Gen T: %i K'%starLabel)
    plt.plot(sed_config.WAVELENGTHS, totalFlux, color='lime',
        label='Total Flux')
    if spitzInsts:
        plt.plot(spitzWaves, spitzFlux, color='black',
            label='Spitzer IRS')

    # Plot the rest of the real data and error bars
    # We don't want distinct labels for upper limits, so we use this list
    # of labels to check if those instruments have already been plotted.
    labels = []
    if totalInsts:
        for inst in totalInsts:
            label = inst.rstrip('0123456789JHK')
            if not label in labels:
                labels.append(label)
            else:
                label = None
            plt.scatter(sData[inst][wa], sData[inst][fx],
                s=20, marker='D', color=plotColors[inst], label=label,
                zorder=10)
            plt.errorbar(sData[inst][wa], sData[inst][fx],
                yerr=sData[inst][er], color=plotColors[inst])
    if ulInsts:
        for inst in ulInsts:
            label = inst.rstrip('0123456789JHKUL')
            if not label in labels:
                labels.append(label)
            else:
                label = None
            plt.scatter(sData[inst][wa], 3*sData[inst][er]+sData[inst][fx],
                s=20, marker='D', color=plotColors[inst], label=label,
                zorder=10)
            plt.errorbar(sData[inst][wa], 3*sData[inst][er]+sData[inst][fx],
                yerr=sData[inst][er], uplims=True, color=plotColors[inst])

    # Plot formatting
    plt.title(starName, fontsize=22)
    plt.xlabel(r'$\lambda$ ($\mu m$)')
    plt.ylabel(r'$F_{\nu}$ ($Jy$)')
    plt.semilogx()
    plt.semilogy()

    # Calculate the y limits based on the data
    if y1.max() < fitFlux[np.where(fitFlux>0)].min() or \
        y2.max() < fitFlux[np.where(fitFlux>0)].min():
        lowerylimit = min([y1.max(), y2.max()])*0.5
    else:
        lowerylimit = fitFlux[np.where(fitFlux>0)].min()*0.5

    if y1.max() > fitFlux.max() or y2.max() > fitFlux.max():
        upperylimit = max([y1.max(), y2.max()])*1.5
    else:
        upperylimit = fitFlux.max()*1.5
    plt.ylim(lowerylimit, upperylimit)

    # The x limits are standard (we only have data for certain wavelengths)
    plt.xlim(.5,300)
    plt.text(0.98, 0.98, r"Reduced $\chi^2$: %0.1f"%(chiSqr),
        ha = 'right', va = 'top', transform = ax.transAxes)

    # Show the resolved radial location (if available and selected)
    if showResolved and sTrigger:
        plt.text(0.98, 0.92, f"{r'r$_{herschel}$'}: ({bbr2} +/- {bbr2_unc}) AU",
            ha = 'right', va = 'top', transform = ax.transAxes)

    # Show the minimum grain size
    if showMinGrain and sTrigger:
        if oneWander:
            plt.text(0.98, 0.86, f"a{r'$_{min}$'}: {bos: 0.4f} {chr(956)}m",
                ha = 'right', va = 'top', transform = ax.transAxes)
        elif twoWander:
            plt.text(0.98, 0.86, f"inner a{r'$_{min}$'}: {bos1: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.80, f"outer a{r'$_{min}$'}: {bos2: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.74, f"inner bos: {blowoutSize1: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.68, f"outer bos: {blowoutSize2: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.62, f"inner {r'f$_{mb}$'}: {bos1/blowoutSize1: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.56, f"outer {r'f$_{mb}$'}: {bos2/blowoutSize2: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
        elif noWander:
            plt.text(0.98, 0.86, f"a{r'$_{min}$'}: {blowoutSize1: 0.4f} {chr(956)}m",
                ha = 'right', va = 'top', transform = ax.transAxes)
        elif twoWarmBelts:
            plt.text(0.98, 0.86, f"small a{r'$_{min}$'}: {bos1: 0.4f} {chr(956)}m",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.80, f"a{r'$_{min}$'}: {blowoutSize1: 0.4f} {chr(956)}m",
                ha = 'right', va = 'top', transform = ax.transAxes)
    elif showMinGrain and not sTrigger:
        if oneWander:
            plt.text(0.98, 0.92, f"a{r'$_{min}$'}: {bos: 0.4f} {chr(956)}m",
                ha = 'right', va = 'top', transform = ax.transAxes)
        elif twoWander:
            plt.text(0.98, 0.92, f"inner a{r'$_{min}$'}: {bos1: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.86, f"outer a{r'$_{min}$'}: {bos2: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.80, f"inner bos: {blowoutSize1: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.74, f"outer bos: {blowoutSize2: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.68, f"inner {r'f$_{mb}$'}: {bos1/blowoutSize1: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.62, f"outer {r'f$_{mb}$'}: {bos2/blowoutSize2: 0.4f}",
                ha = 'right', va = 'top', transform = ax.transAxes)
        elif noWander:
            plt.text(0.98, 0.92, f"a{r'$_{min}$'}: {blowoutSize1: 0.4f} {chr(956)}m",
                ha = 'right', va = 'top', transform = ax.transAxes)
        elif twoWarmBelts:
            plt.text(0.98, 0.92, f"small a{r'$_{min}$'}: {bos1: 0.4f} {chr(956)}m",
                ha = 'right', va = 'top', transform = ax.transAxes)
            plt.text(0.98, 0.86, f"a{r'$_{min}$'}: {blowoutSize1: 0.4f} {chr(956)}m",
                ha = 'right', va = 'top', transform = ax.transAxes)

    plt.legend(loc='lower left')
    if showIRSVariance:
        ax = fig.add_subplot(212)
        indices = np.where(np.logical_and(fitWaves>= spitzWaves[0], fitWaves<= spitzWaves[-1]))
        variance = fitFlux - warmColdWander(fitWaves, RW, RC, bos1, bos2, n1, n2, n3)
        plt.plot(fitWaves[indices], variance[indices])
        plt.axhline()

    if saveFigure:
        plt.savefig(IMG_DIR+starName+'.png', bbox_inches='tight')
    if showFigure:
        plt.show()
    plt.close('all')

################################################################################
####                Run the code on the list of star names
################################################################################

if __name__ == '__main__':

    # starNames = starNames[7:]
    for stn, starName in enumerate(starNames):
        # Choose the bounds on the blowout sizes given the conditional
        # statements in the bound_options.py file.
        bosLowInner, bosHighInner, bosLowOuter, bosHighOuter = (
            bound_options.boundify(starName)
            )

        # Let's the user know which star is running and how many are being
        # processed at a time.
        print("++++++++++++++++++++++++++++++++++++++++")
        print(f" Running fit #{stn+1} of {len(starNames)} for star: {starName}")
        print("++++++++++++++++++++++++++++++++++++++++")
        # Call the fitting module on the star name.
        run_fits(starName)
