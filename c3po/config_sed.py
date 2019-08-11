# Functions for SED fitting
import os
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import constants
from scipy.interpolate import UnivariateSpline as uni_spline
from astropy.io import ascii

# This tells NumPy to stop issuing warnings everytime an invalid value is
# encountered in a division (zero division error). NumPy still puts NaN
# where this occurs but it happens silently! At last.
np.seterr(all='ignore')

# Make the path operating system independent
def to_path(path_string):
    return os.sep.join( path_string.split('/') )

# File Paths
ARR_DIR         = to_path('Data/Arrays/')
STAR_FILES      = to_path('Data/StarFiles/')
KURUCZ          = to_path('Data/StellarModels/kurucz/')
NEXTGEN         = to_path('Data/StellarModels/nextgen/')
RES_DIR         = to_path('Results/')
PLOTS_DIR       = to_path('Results/Plots/')
PARAMS_DIR      = to_path('Results/Params/')
FRATIOS_DIR     = to_path('Results/FluxRatios/')
FARRAYS_DIR     = to_path('Results/FluxArrays/')
INTERPS_DIR     = to_path('Data/Arrays/InterpGrainTemps/')
FILTERS_DIR     = to_path('Data/FilterResponse/')
GRAIN_TEMPS_DIR = to_path('Data/GrainTemps/')
FLUX_CUBES_DIR  = to_path('Data/FluxCubes/')

# Frequently used arrays
GRAINSIZES = np.loadtxt(ARR_DIR + 'GrainSizes.dat')
WAVELENGTHS = np.loadtxt(ARR_DIR + 'Wavelengths.dat')
STAR_TEMPS = np.linspace(2000, 15000, 14)
WAVES = np.logspace(-3, 3, 1000)
TEMPS_RADII = np.logspace(-3, 3, 1000)

# Grain temperatures per grain composition
# Deprecated soon?
grainComps = ['AstroSil', 'DirtyIceAstroSil']
# GRAIN_TEMPS_TOTAL = dict()
EMISSIVITIES_TOTAL = dict()
for grain in grainComps:
    EMISSIVITIES_TOTAL[grain] = np.load(ARR_DIR+grain+'Emissivities.npy')

################################################################################
#                                                                     Functions
################################################################################
# Find nearest value in an array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

# Find the array of the nearest element in an array
def find_nearest_ind(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# Wrapper for the IPAC reader function
def read_star_file(starFile):
    '''
    Returns a dictionary with all of the comments included as keys
    '''
    x = ascii.read(starFile)
    data = {}
    for column in x.itercols():
        data[column.name] = np.array(column)
    for key in x.meta['keywords'].keys():
        data[key] = x.meta['keywords'][key]['value']
    return data

def stitch_spitz(allSpitz):
    # Organize by ascending wavelength
    # idx = np.argsort(allSpitz['wavelength'])
    # wav = allSpitz['wavelength'][idx]
    wav, idx = np.unique(allSpitz['wavelength'], return_index=True)
    flx = allSpitz['flux'][idx]
    err = allSpitz['error'][idx]
    ist = allSpitz['instrument'][idx]
    # Trim to less than 35 microns wavelengths
    flx = flx[wav<35]
    err = err[wav<35]
    ist = ist[wav<35]
    wav = wav[wav<35]
    # Grab the unique Spitzer modules
    uist = np.unique(ist)
    # Options for the stitching. cutsize refers to how much of the module is
    # used to construct the UnivariateSpline. e.g. 2 == one half of the data.
    # pt2 and pt1 are the points in the arrays where the midpoints begin
    cutsize =  2
    pt2     = -4
    pt1     =  3
    # Stitch from right to left
    if 'SpitzerIRS-LL2' in uist and 'SpitzerIRS_LL1' in uist:
        lidx = np.where(ist=='SpitzerIRS-LL2')
        ridx = np.where(ist=='SpitzerIRS-LL1')
        sl = int(wav[lidx].size/3) # last part of left module
        sr = int(wav[ridx].size/3) # first part of right module
        # sl = int(wav[lidx].size/cutsize) # last part of left module
        # sr = int(wav[ridx].size/cutsize) # first part of right module
        leftfit = uni_spline(wav[lidx][sl:], flx[lidx][sl:])
        rightfit = uni_spline(wav[ridx][:sr], flx[ridx][:sr])
        midpoints = np.linspace(wav[lidx][pt2], wav[ridx][pt1], 100)
        lefty = leftfit(midpoints)
        righty = rightfit(midpoints)
        stitchfactor = righty.sum()/lefty.sum()
        # print("ll2 to ll1: {}".format(stitchfactor))
        flx[lidx] *= stitchfactor
    if 'SpitzerIRS-SL1' in uist and 'SpitzerIRS-LL2' in uist:
        lidx = np.where(ist=='SpitzerIRS-SL1')
        ridx = np.where(ist=='SpitzerIRS-LL2')
        sl = int(wav[lidx].size/cutsize) # last part of left module
        sr = int(wav[ridx].size/cutsize) # first part of right module
        leftfit = uni_spline(wav[lidx][sl:], flx[lidx][sl:])
        rightfit = uni_spline(wav[ridx][:sr], flx[ridx][:sr])
        midpoints = np.linspace(wav[lidx][pt2], wav[ridx][pt1], 100)
        lefty = leftfit(midpoints)
        righty = rightfit(midpoints)
        stitchfactor = righty.sum()/lefty.sum()
        # print("sl1 to ll2: {}".format(stitchfactor))
        flx[lidx] *= stitchfactor
    if 'SpitzerIRS-SL2' in uist and 'SpitzerIRS-SL1' in uist:
        lidx = np.where(ist=='SpitzerIRS-SL2')
        ridx = np.where(ist=='SpitzerIRS-SL1')
        sl = int(wav[lidx].size/cutsize) # last part of left module
        sr = int(wav[ridx].size/cutsize) # first part of right module
        leftfit = uni_spline(wav[lidx][sl:], flx[lidx][sl:])
        rightfit = uni_spline(wav[ridx][:sr], flx[ridx][:sr])
        midpoints = np.linspace(wav[lidx][pt2], wav[ridx][pt1], 100)
        lefty = leftfit(midpoints)
        righty = rightfit(midpoints)
        stitchfactor = righty.sum()/lefty.sum()
        # print("sl2 to sl1: {}".format(stitchfactor))
        flx[lidx] *= stitchfactor
    # Finally return the stitched data
    allSpitzStitched = {
        'wavelength': wav, 'flux': flx, 'error': err, 'instrument': ist
    }
    return allSpitzStitched

def sort_spitz_data(starData):
    # Sort Data: Pull SpitzerIRS
    spindexes = np.concatenate((
        np.nonzero('SpitzerIRS-SL2' == starData['instrument'])[0],
        np.nonzero('SpitzerIRS-SL1' == starData['instrument'])[0],
        np.nonzero('SpitzerIRS-LL2' == starData['instrument'])[0],
        np.nonzero('SpitzerIRS-LL1' == starData['instrument'])[0]
        ))
    # if len(spindexes) == 0:
    #     return [], starData
    spitzWaves = starData['wavelength'][spindexes]
    spitzFlux  = starData['flux'][spindexes]
    spitzError = starData['error'][spindexes]
    spitzInsts = starData['instrument'][spindexes]
    allSpitz   = {'wavelength': spitzWaves[spitzWaves<35],
        'flux': spitzFlux[spitzWaves<35], 'error': spitzError[spitzWaves<35],
        'instrument': spitzInsts[spitzWaves<35]}
    # Sort Data: Pull all other Photometry Points
    phindexes = np.nonzero(np.logical_and(
        np.logical_and('SpitzerIRS-SL2' != starData['instrument'],
            'SpitzerIRS-SL1' != starData['instrument']),
        np.logical_and('SpitzerIRS-LL2' != starData['instrument'],
            'SpitzerIRS-LL1' != starData['instrument'])
        ))
    totalWaves = starData['wavelength'][phindexes]
    totalFlux  = starData['flux'][phindexes]
    totalError = starData['error'][phindexes]
    totalInsts = starData['instrument'][phindexes]
    nonSpitz   = {'wavelength': totalWaves, 'flux': totalFlux,
        'error': totalError, 'instrument': totalInsts}
    return allSpitz, nonSpitz

def sort_wise_data(starData):
    # Separate the AllSky data
    skindexes = np.concatenate((
        np.nonzero('AllSkyW1' == starData['instrument'])[0],
        np.nonzero('AllSkyW2' == starData['instrument'])[0],
        np.nonzero('AllSkyW3' == starData['instrument'])[0],
        np.nonzero('AllSkyW4' == starData['instrument'])[0]
        ))
    skyWaves = starData['wavelength'][skindexes]
    skyFlux  = starData['flux'][skindexes]
    skyError = starData['error'][skindexes]
    skyInsts = starData['instrument'][skindexes]
    allSky = {'wavelength': skyWaves, 'flux': skyFlux,
        'error': skyError, 'instrument': skyInsts}
    # Separate the AllWise data
    windexes = np.concatenate((
        np.nonzero('AllWiseW1' == starData['instrument'])[0],
        np.nonzero('AllWiseW2' == starData['instrument'])[0],
        np.nonzero('AllWiseW3' == starData['instrument'])[0],
        np.nonzero('AllWiseW4' == starData['instrument'])[0]
        ))
    allWiseWaves = starData['wavelength'][windexes]
    allWiseFlux  = starData['flux'][windexes]
    allWiseError = starData['error'][windexes]
    allWiseInsts = starData['instrument'][windexes]
    allWise = {'instrument': allWiseWaves, 'flux': allWiseFlux,
        'error': allWiseError, 'instrument': allWiseInsts}
    remaindexes = np.nonzero(
        np.logical_and(
            np.logical_and(
                np.logical_and('AllSkyW1' != starData['instrument'],
                    'AllSkyW2' != starData['instrument']),
                np.logical_and('AllSkyW3' != starData['instrument'],
                    'AllSkyW4' != starData['instrument'])
                ),
            np.logical_and(
                np.logical_and('AllWiseW1' != starData['instrument'],
                    'AllWiseW2' != starData['instrument']),
                np.logical_and('AllWiseW3' != starData['instrument'],
                    'AllWiseW4' != starData['instrument'])
                )
            )
        )
    nonWiseWaves = starData['wavelength'][remaindexes]
    nonWiseFlux  = starData['flux'][remaindexes]
    nonWiseError = starData['error'][remaindexes]
    nonWiseInsts = starData['instrument'][remaindexes]
    nonWise = {'wavelength': nonWiseWaves, 'flux': nonWiseFlux,
        'error': nonWiseError, 'instrument': nonWiseInsts}
    return allSky, allWise, nonWise

# Convolution function for calibrating data
def convolve(filterwaves, filterresponse, datawaves, dataflux):
    top = integrate.simps(
      np.interp(filterwaves, datawaves, dataflux)*filterresponse,
        filterwaves)
    bottom = integrate.simps(filterresponse, filterwaves)
    return top/bottom

def load_stellar_model(starT):
    '''Load the appropriate NextGen stellar model that is closest to the given
    temperature of the star. If the temperature exceeds 10,000 degrees Kelvin,
    then stitch together a NextGen model with a Kurucz model.

    Returns
    value       :  units
    --------------------
    wavelengths : micron
           flux : jansky
      starLabel : Kelvin'''
    # Concerning the nextgen and kurucz models: the wavelengths are in units of
    # angstrom, which equals 10^-10 meters. Flux is in erg/cm^2/s/a.
    # Create temp array that matches the temperatures of the nextgen stellar
    # models. From 2600-4100, uses a step of 100. From 4200-10000, uses a step
    # of 200. Grabs the nearest model to the star temp in the star file.
    ngTemps    = np.arange(2600, 4100, 100)
    ngTemps    = np.append(ngTemps, np.arange(4200, 10200, 200))
    TEMP       = find_nearest(ngTemps, starT)
    ngfilename = 'xp00_'+str(TEMP)+'g40.txt'
    # Temperature label for the graph of the SED fit
    starLabel = TEMP
    # Load the nextgen stellar model
    ngWave, ngFlux   = np.loadtxt(NEXTGEN+ngfilename, unpack=True)
    ngWave, ngwIndex = np.unique(ngWave, return_index=True)
    ngFlux           = ngFlux[ngwIndex]
    # Convert the nextgen model to janskies
    c_cgs  = 2.99792458e10                          #cm/s
    ngFnu  = ngFlux*(ngWave**2) / c_cgs / 1.0e8 * 1.0e23  #Jy
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
        kzWave, kzFlux = np.loadtxt(KURUCZ+kzfilename, unpack=True)
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
        ngWave = np.append(kzWave, ngWave)
        ngFnu = np.append(kzFnu, ngFnu)
    return ngWave, ngFnu, starLabel

def calc_luminosity(ngWave, ngFnu, starD):
    '''Requires ngWave in microns, ngFnu in Jy, starD in parsecs.

    Returns in units of solar luminosities
    '''
    absolute_lum = ngFnu/1e23 * 4. * np.pi * (starD*3.086e18)**2
    absolute_lum = -integrate.simps(absolute_lum, (3e10/(ngWave/1e4)))
    absolute_lum /= 3.828e33
    return absolute_lum

# Blackbody radiation function, in terms of frequency
def b_nu(wavelengths, temperature):
    '''Blackbody radiation per frequency. Expects input wavelengths to be in
    units of microns and temperature in put to be in degrees Kelvin.

    Returns in units of erg / (sr * m^2 * Hz)'''
    H = constants.h
    C = constants.c
    K = constants.k
    wavelengths_m = wavelengths/1.0e6
    nu = C/wavelengths_m
    top = 2 * H * (nu**3)
    exponent = (H*nu)/(K*temperature)
    bottom = C**2 * (np.exp(exponent) - 1)
    return top/bottom

# Blackbody radiation function, in terms of wavelength
def b_lam(wavelengths, temperature):
    '''Blackbody radiation per wavelength. Expects input wavelengths to be in
    units of microns and temperature in put to be in degrees Kelvin.

    Returns in units of W / (sr * m^3)'''
    H = constants.h
    C = constants.c
    K = constants.k
    wavelengths_m = wavelengths/1.0e6
    top = 2 * H * C**2
    exponent = (H*C)/(wavelengths_m*K*temperature)
    bottom = wavelengths_m**5 * (np.exp(exponent) - 1)
    return top/bottom

# Calculate the minimum blowout size given stellar properties
def blowout_size(grainDensity, starL=1., starM=1., qRad=0.9):
    '''Calculate blowoutsize given the properties of a star.

    Returns in units of microns.'''
    # Calculate blowout grain size
    density = grainDensity / 1000 / (0.01**3) # Density converted to SI
    nume = 3 * starL * 3.826e26 * qRad
    deno = 8 * np.pi*starM* 1.989e30 *constants.G*constants.c*density
    return nume/deno * 1e6 # microns

def lum_ratio(starFlux, dustFlux, waves):
    '''Calculate the luminosity ratio of the given dust belt to the star.

    Returns in no units, because it's a ratio.
    '''
    l_dust = integrate.simps(dustFlux, (3e14/waves))
    l_star = integrate.simps(starFlux, (3e14/waves))
    return l_dust/l_star

def flux_ratios_herschel(y1=None, y2=None, ngModel=None, y3=None):
    '''Calculate the flux ratios for Herschel frequency reponse functions.
    '''
    # Herschel frequency response functions
    h70w, h70r = np.loadtxt(FILTERS_DIR+'pacs-blue-70.dat', unpack=True)
    h100w, h100r = np.loadtxt(FILTERS_DIR+'pacs-green-100.dat', unpack=True)
    h160w, h160r = np.loadtxt(FILTERS_DIR+'pacs-red-160.dat', unpack=True)

    # Calc flux for the star
    fh70star  = convolve(h70w, h70r, WAVELENGTHS, ngModel)
    fh100star = convolve(h100w, h100r, WAVELENGTHS, ngModel)
    fh160star = convolve(h160w, h160r, WAVELENGTHS, ngModel)

    # Calculate flux ratios for the warm belt
    if np.any(y1):
        fh70warm  = convolve(h70w, h70r, WAVELENGTHS, y1)/fh70star
        fh100warm = convolve(h100w, h100r, WAVELENGTHS, y1)/fh100star
        fh160warm = convolve(h160w, h160r, WAVELENGTHS, y1)/fh160star

    # Calculate flux ratios for the cold belt
    fh70cold  = convolve(h70w, h70r, WAVELENGTHS, y2)/fh70star
    fh100cold = convolve(h100w, h100r, WAVELENGTHS, y2)/fh100star
    fh160cold = convolve(h160w, h160r, WAVELENGTHS, y2)/fh160star

    if np.any(y3):
        # y3 is only for fits with two co-located warm belts
        fh70warmsmall   = convolve(h70w, h70r, WAVELENGTHS, y3)
        fh70warmsmall  /= fh70star
        fh100warmsmall  = convolve(h100w, h100r, WAVELENGTHS, y3)
        fh100warmsmall /= fh100star
        fh160warmsmall  = convolve(h160w, h160r, WAVELENGTHS, y3)
        fh160warmsmall /= fh160star
    else:
        fh70warmsmall, fh100warmsmall, fh160warmsmall = np.nan, np.nan, np.nan
    if not np.any(y1):
        fh70warm, fh100warm, fh160warm = [np.nan]*3
    return (fh70warm, fh70cold,  fh100warm, fh100cold, fh160warm, fh160cold,
        fh70warmsmall, fh100warmsmall, fh160warmsmall)

def flux_ratios_mips(y1=None, y2=None, ngModel=None, y3=None):
    '''Calculate the flux ratios for Spitzer MIPS frequency response functions.
    '''
    # MIPS frequency response functions
    m24w, m24r = np.loadtxt(FILTERS_DIR+'mips24_frf.txt', unpack=True)
    m70w, m70r = np.loadtxt(FILTERS_DIR+'mips70_frf.txt', unpack=True)
    m160w, m160r = np.loadtxt(FILTERS_DIR+'mips160_frf.txt', unpack=True)
    # Calc flux for the star
    fm24star  = convolve(m24w, m24r, WAVELENGTHS, ngModel)
    fm70star  = convolve(m70w, m70r, WAVELENGTHS, ngModel)
    fm160star = convolve(m160w, m160r, WAVELENGTHS, ngModel)

    # Calculate flux ratios for the warm belt
    if np.any(y1):
        fm24warm  = convolve(m24w, m24r, WAVELENGTHS, y1)/fm24star
        fm70warm  = convolve(m70w, m70r, WAVELENGTHS, y1)/fm70star
        fm160warm = convolve(m160w, m160r, WAVELENGTHS, y1)/fm160star

    # Calculate flux ratios for the cold belt
    fm24cold  = convolve(m24w, m24r, WAVELENGTHS, y2)/fm24star
    fm70cold  = convolve(m70w, m70r, WAVELENGTHS, y2)/fm70star
    fm160cold = convolve(m160w, m160r, WAVELENGTHS, y2)/fm160star

    if np.any(y3):
        # y3 is only for fits with two co-located warm belts
        fm24warmsmall = convolve(m24w, m24r, WAVELENGTHS, y3)/fm24star
        fm70warmsmall = convolve(m70w, m70r, WAVELENGTHS, y3)/fm70star
        fm160warmsmall = convolve(m160w, m160r, WAVELENGTHS, y3)/fm160star
    else:
        fm24warmsmall, fm70warmsmall, fm160warmsmall = np.nan, np.nan, np.nan
    if not np.any(y1):
        fm24warm, fm70warm, fm160warm = [np.nan]*3
    return (fm24warm, fm24cold, fm70warm, fm70cold, fm160warm, fm160cold,
        fm24warmsmall, fm70warmsmall, fm160warmsmall)

# Star object handles all of the flux calculations
class Star:
    '''We use this Star object to hold frequently accessed data for the various
    functions used in SED fitting. The use of this object is to enable
    computational efficiency while also speeding up calculation completion time.'''
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
        'WaterIce': 1,
        # 'AstroSil': 3.0,           # Correlates to DIAS -> 2.034 g/cm^3
        # 'DirtyIceAstroSil': 2.034, # Correlates to AS -> 3.0 g/cm^3
        }
    RHO_AS   = densities['AstroSil']
    RHO_DIAS = densities['DirtyIceAstroSil']
    GRAINS   = GRAINSIZES/1e6
    WAVES    = WAVELENGTHS/1e6
    RADII    = TEMPS_RADII*1.4959787066e11
    EMIS     = EMISSIVITIES_TOTAL
    A_MAX    = 0.001 # 1000 microns == 1 mm == 0.001 m
    SIGMA    = 0.1
    Q        = -3.5

    def __init__(self, name, starD, starL, starT, starM, gTemps, emis=None, fitWaves=None):
        # Input data
        self.starD      = starD
        self.starL      = starL
        self.starT      = starT
        self.starM      = starM
        self.emis       = emis
        self.fitWaves   = fitWaves
        self.gTemps     = gTemps
        self.name       = name

        # Calculate blowout sizes/indexes
        self.bos1       = blowout_size(  self.RHO_AS, starL, starM)
        self.bos2       = blowout_size(self.RHO_DIAS, starL, starM)
        self.graindex1  = find_nearest_ind(self.GRAINS, self.bos1/1e6)
        self.graindex2  = find_nearest_ind(self.GRAINS, self.bos2/1e6)

        print('****************************************')
        print('       Initializing flux cubes')
        # Compute 'flux cubes' for calculations
        # Low resolutions (at only the fitting wavelengths)
        self.warm_high_res_cube = self.init_flux_cube('AstroSil', 'high')
        self.cold_high_res_cube = self.init_flux_cube('DirtyIceAstroSil','high')
        # Low resolutions (at the full wavelengths array)
        if emis is not None and fitWaves is not None:
            self.warm_low_res_cube  = self.init_flux_cube('AstroSil', 'low')
            self.cold_low_res_cube  = self.init_flux_cube('DirtyIceAstroSil', 'low')
        print('           All cubes loaded')
        print('****************************************')


    def init_flux_cube(self, grainComp, resolution):
        '''Calculate the flux cube a priori. Uses all grainsizes, all radii, and
        all wavelengths (depending on fitWaves vs standard wavelengths array).
        Later, this cube will be accessed using the calculated indexes from the
        function call.

        Shape returned:
            (# of grains, # of radii, # of wavelengths): (121, 1000, 901)'''
        fname = '{}_{}_{}.npy'.format(self.name, grainComp, resolution)
        if resolution == 'low':
            waves1 = self.fitWaves[np.newaxis, :]
            emis1  = self.emis[grainComp][:, np.newaxis, :]
        elif resolution == 'high':
            waves1 = self.WAVES[np.newaxis, :]*1e6 # want in Microns
            emis1  = self.EMIS[grainComp][:, np.newaxis, :]
        grains1    = self.GRAINS[:, np.newaxis, np.newaxis]
        temps1     = self.gTemps[grainComp][np.newaxis,:,:].T
        radii1     = self.RADII[np.newaxis, :, np.newaxis]
        flux_cube  = emis1 * (1./grains1/np.sqrt(grains1)) * radii1
        flux_cube *= b_nu(waves1, temps1)
        return flux_cube

    def calcDustMass(self,  n_dust, r_0, a_min, a_max, density):
        '''
        Calculate the mass of the small grains in the dust belt.

        Returns in units of kg
        '''
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.RADII<1.4*r0, self.RADII>0.6*r0))
        radii1 = self.RADII[rindex]
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        asdf = np.exp(exponent) * 2.0 * np.pi * radii1
        money = integrate.simps(asdf, radii1) * (4./3.) * (density*1e3)
        dust_mass = np.sqrt(a_min*a_max) * money * n_dust
        return dust_mass

    def low_res_warm_small(self, waves, _r0, _amin1):
        '''
        This function is included so that we can do stacking analysis of the
        small grains from the belt. Shouldn't be used in fitting at all.
        '''
        r0 = _r0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS & a_min to SI units
        amin1 = _amin1 / 1e6
        bos1 = self.bos1 / 1e6
        # Find index of small grain
        graindex = find_nearest_ind(self.GRAINS, amin1)
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( amin1**(3+self.Q) - (bos1)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.warm_low_res_cube[graindex:self.graindex2, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[graindex:self.graindex2], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def high_res_warm_small(self, waves, _r0, _amin1):
        '''
        This function is included so that we can do stacking analysis of the
        small grains from the belt. Shouldn't be used in fitting at all.
        '''
        r0 = _r0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS & a_min to SI units
        amin1 = _amin1 / 1e6
        bos1 = self.bos1 / 1e6
        # Find index of small grain
        graindex = find_nearest_ind(self.GRAINS, amin1)
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( amin1**(3+self.Q) - (bos1)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.warm_high_res_cube[graindex:self.graindex2, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[graindex:self.graindex2], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def low_res_warm_amin(self, waves, _r0, _amin1):
        r0 = _r0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS to SI units
        amin1 = _amin1 / 1e6
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( amin1**(3+self.Q) - (self.A_MAX)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.warm_low_res_cube[self.graindex2:, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[self.graindex2:], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def high_res_warm_amin(self, waves, _r0, _amin1):
        r0 = _r0*1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS to SI units
        amin1 = _amin1 / 1e6
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( amin1**(3+self.Q) - (self.A_MAX)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.cold_high_res_cube[self.graindex2:, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[self.graindex2:], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def low_res_cold_amin(self, waves, _r0, _amin2):
        r0 = _r0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS to SI units
        amin2 = _amin2 / 1e6
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( amin2**(3+self.Q) - (self.A_MAX)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.cold_low_res_cube[self.graindex2:, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[self.graindex2:], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def high_res_cold_amin(self, waves, _r0, _amin2):
        r0 = _r0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS to SI units
        amin2 = _amin2 / 1e6
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( amin2**(3+self.Q) - (self.A_MAX)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.cold_high_res_cube[self.graindex2:, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[self.graindex2:], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def low_res_warm(self, waves, r_0):
        r0 = r_0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS to SI units
        bos = self.bos2/1e6
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( bos**(3+self.Q) - (self.A_MAX)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.warm_low_res_cube[self.graindex2:, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[self.graindex2:], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def high_res_warm(self, waves, r_0):
        r0 = r_0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS to SI units
        bos = self.bos2/1e6
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( bos**(3+self.Q) - (self.A_MAX)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.warm_high_res_cube[self.graindex2:, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[self.graindex2:], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def low_res_cold(self, waves, r_0):
        r0 = r_0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS to SI units
        bos = self.bos2/1e6
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( bos**(3+self.Q) - (self.A_MAX)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.cold_low_res_cube[self.graindex2:, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[self.graindex2:], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def high_res_cold(self, waves, r_0):
        r0 = r_0 *1.4959787066e11
        # Find the index of radii within the thickness of the belt
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        # Convert BOS to SI units
        bos = self.bos2/1e6
        # Calculate the distribution of grains per radii
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = np.exp(exponent) / ( bos**(3+self.Q) - (1e-3)**(3+self.Q) )
        ca1 = ca[np.newaxis, :, np.newaxis]
        # Multiply the integrand using the flux cube computed a priori
        f = ca1 * self.cold_high_res_cube[self.graindex2:, rindex, :]
        # Integrate with respect to grains. Collapses shape -> (#radii, #waves)
        f = integrate.simps(f, self.GRAINS[self.graindex2:], axis=0)
        # Integrate with respect to radii. Collapses shape -> (#waves)
        f  = integrate.simps(f, self.RADII[rindex], axis=0)
        # Return the result multiplied by the derived constant
        return f*3.299e-7/(self.starD**2)

    def deprecated_calcFluxCold(self, waves, r_0, T_0 = 1):
        '''
        This is kept for use as a comparison to the new algorithm.
        '''
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.RADII<1.4*r0,
            self.RADII>0.6*r0))[0]
        radii1 = self.RADII[rindex]
        grainTemps = self.gTemps['DirtyIceAstroSil'][rindex]
        grains = self.GRAINS[self.graindex2:]
        bos = self.bos2/1e6
        exponent = -0.5 * ((radii1 - r0) / (self.SIGMA*r0))**2
        ca = T_0*np.exp(exponent)\
            / (np.power(bos,3+self.Q)-np.power(.001,3+self.Q))

        # Nested loop in Python
        fw = np.empty(waves.size)
        fr = np.empty(radii1.size)
        for w in range(waves.size):
            for r in range(radii1.size):
                flux = b_nu(waves[w], grainTemps[r, self.graindex2:])
                flux *= self.EMIS['DirtyIceAstroSil'][self.graindex2:,w] * \
                    (grains**-1.5) * radii1[r] * ca[r]
                fr[r] = integrate.simps(flux, grains)
            fw[w] = integrate.simps(fr, radii1)
        return fw*3.299e-7/(self.starD**2)
