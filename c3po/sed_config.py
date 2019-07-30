# Functions for SED fitting
import os
import numpy as np
import pandas as pd
from scipy import integrate
from scipy import constants
from scipy.interpolate import UnivariateSpline as uni_spline
from astropy.io import ascii


# File Paths
ARR_DIR     = os.sep.join('Data/Arrays/'.split('/'))
# STAR_FILES  = os.sep.join('Data/IPAC Files/'.split('/'))
STAR_FILES  = os.sep.join('Data/StarFiles/'.split('/'))
KURUCZ      = os.sep.join('Data/StellarModels/kurucz/'.split('/'))
NEXTGEN     = os.sep.join('Data/StellarModels/nextgen/'.split('/'))
RES_DIR     = os.sep.join('Results/'.split('/'))
PLOTS_DIR   = os.sep.join('Results/Plots/'.split('/'))
PARAMS_DIR  = os.sep.join('Results/Params/'.split('/'))
FRATIOS_DIR = os.sep.join('Results/FluxRatios/'.split('/'))
FARRAYS_DIR = os.sep.join('Results/FluxArrays/'.split('/'))
INTERPS_DIR = os.sep.join('Data/Arrays/InterpGrainTemps/'.split('/'))
FILTERS_DIR = os.sep.join('Data/FilterResponse/'.split('/'))
GRAIN_TEMPS_DIR = os.sep.join('Data/GrainTemps/'.split('/'))

# Frequently used arrays
GRAINSIZES = np.loadtxt(ARR_DIR + 'GrainSizes.dat')
WAVELENGTHS = np.loadtxt(ARR_DIR + 'Wavelengths.dat')
STAR_TEMPS = np.linspace(2000, 15000, 14)
DISK_RADII = np.logspace(-1, 3, 121)
WAVES = np.logspace(-3, 3, 1000)
TEMPS_RADII = np.logspace(-1, 3, 1000)

# Grain temperatures per grain composition
# Deprecated soon?
grainComps = ['AstroSil', 'DirtyIceAstroSil']
GRAIN_TEMPS_TOTAL = dict()
EMISSIVITIES_TOTAL = dict()
for grain in grainComps:
    GRAIN_TEMPS_TOTAL[grain] = np.load(ARR_DIR+grain+'GrainTemps.npy')
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
    '''
    Requires ngWave in microns, ngFnu in Jy, starD in parsecs.
    Returns in units of solar luminosities
    '''
    absolute_lum = ngFnu/1e23 * 4. * np.pi * (starD*3.086e18)**2
    absolute_lum = -integrate.simps(absolute_lum, (3e10/(ngWave/1e4)))
    absolute_lum /= 3.828e33
    return absolute_lum

# Normalize a black body
def norm_blackbodies(sData, fitWaves, nWarm, nCold, t_1, t_2):
    wa, fx = 'wavelength', 'flux'
    er = 'error'
    w4 = 'WISE4'
    mp24, mp70 = 'MIPS24', 'MIPS70'
    h70, h100 = 'HerschelPACS70', 'HerschelPACS100'
    h160 = 'HerschelPACS160'
    ll1  = 'SpitzerIRS-LL1'
    # Which values are present in the data?
    mp24f, ll1f, w4f = nWarm
    mp70f, h70f, h100f, h160f = nCold
    wav = np.logspace(-1, 2.3, 1000) # 0.1 (um) to 199.5 (um)
    bb1 = b_lam(wav, t_1)
    bb2 = b_lam(wav, t_2)
    # Normalize warm dust
    if mp24f:
        n_1 = np.nanmean(sData[mp24][fx]) / bb1.max()
    elif ll1f:
        index2 = np.where(np.logical_and(sData[ll1][wa]>20, sData[ll1][wa]<25))
        n_1 = np.nanmean(sData[ll1][fx][index2]) / bb1.max()
    elif w4f:
        n_1 = np.nanmean(sData[w4][fx]) / bb1.max()
    else:
        n_1 = 1
    # Normalize cold dust
    if mp70f:
        n_2 = np.nanmean(sData[mp70][fx]) / bb2.max()
    elif h70f:
        n_2 = np.nanmean(sData[h70][fx]) / bb2.max()
    elif h100f:
        n_2 = np.nanmean(sData[h100][fx]) / bb2.max()
    elif h160f:
        n_2 = np.nanmean(sData[h160][fx]) / bb2.max()
    else:
        n_2 = n_1
    return n_1, n_2

# Blackbody radiation function, in terms of frequency
def b_nu(wavelengths, temperature):
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
    # Calculate blowout grain size
    density = grainDensity / 1000 / (0.01**3) # Density converted to SI
    nume = 3 * starL * 3.826e26 * qRad
    deno = 8 * np.pi*starM* 1.989e30 *constants.G*constants.c*density
    return nume/deno * 1e6 # microns

# Separates data according to instrument. Input: dict of IPAC Tables
def sort_by_instrument(data):
    wa = 'wavelength'
    fx = 'flux'
    er = 'error'
    unique_insts = np.unique(data['instrument'])
    separatedStarData = dict()
    for inst in unique_insts:
        index = np.where(data['instrument']==inst)
        separatedStarData[inst] = dict()
        separatedStarData[inst][wa] = np.array(data[wa][index])
        separatedStarData[inst][fx] = np.array(data[fx][index])
        separatedStarData[inst][er] = np.array(data[er][index])
    return unique_insts, separatedStarData

def lum_ratio(starFlux, dustFlux, waves):
    l_dust = integrate.simps(dustFlux, (3e14/waves))
    l_star = integrate.simps(starFlux, (3e14/waves))
    return l_dust/l_star

def flux_ratios_herschel(y1=None, y2=None, ngModel=None, y3=None):
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
    def __init__(self, starD, starL, starT, gTemps, blowoutSize1, blowoutSize2,
        emis, grains, graindex1, graindex2):
        self.starD = starD
        self.starL = starL
        self.starT = starT
        self.grainTemps = gTemps
        self.blowoutSize1 = blowoutSize1
        self.blowoutSize2 = blowoutSize2
        self.emis = emis
        self.grains = grains
        self.graindex1 = graindex1
        self.graindex2 = graindex2
        self.radii = np.logspace(-1, 3, 1000)*1.4959787066e11
        self.sigma = 0.1
        self.q = -3.5     # collisional cascade
        # For the flux calculations, use q = -3.5, so the constant at the end
        # is affected by this number but is hardcoded in
        self.graindist = self.makeDist(0.6, 1.4)

    def makeDist(self, low, high):
        # this is always the same regardless of radial location
        x = np.linspace(low, high, 1000)
        return np.exp( -50 * (1 - x)**2 )

    def calcDustMass(self,  n_dust, r_0, a_min, a_max, density):
        '''
        Calculate the mass of the small grains in the dust belt.
        '''
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        exponent = -0.5 * ((radii1 - r0) / (self.sigma*r0))**2
        # print(np.exp(exponent))
        asdf = np.exp(exponent) * 2.0 * np.pi * radii1
        money = integrate.simps(asdf, radii1) * (4./3.) * (density*1e3)
        dust_mass = np.sqrt(a_min*a_max) * money * n_dust
        return dust_mass

    def calcFluxWarmMinGrains(self, waves, r_0, amin1, T_0=1):
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['AstroSil'][rindex]
        graindex = find_nearest_ind(GRAINSIZES, amin1)
        grains = GRAINSIZES[graindex:self.graindex1]/1e6
        amin1 /= 1e6
        exponent = -0.5 * ((radii1 - r0) / (self.sigma*r0))**2
        ca = T_0*np.exp(exponent) \
            / (np.power(amin1,3+self.q)-np.power(self.blowoutSize1,3+self.q))
        ca1 = np.reshape(ca, (1, ca.size, 1))
        grains1 = np.reshape(grains, (grains.size, 1, 1))
        waves1 = np.broadcast_to(waves, (1, waves.size))
        temps1 = np.broadcast_to(
            grainTemps[:,graindex:self.graindex1],
                (1, grainTemps[:,graindex:self.graindex1].shape[0],
                grainTemps[:,graindex:self.graindex1].shape[1])).T
        emis1 = np.reshape(self.emis['AstroSil'][graindex:self.graindex1],
            (self.emis['AstroSil'][graindex:self.graindex1].shape[0], 1,
            self.emis['AstroSil'][graindex:self.graindex1].shape[1]))
        radii1 = np.reshape(radii1, (1, radii1.size, 1))
        flux = emis1 * ca1 * grains1**(self.q+2) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(integrate.simps(flux, grains, axis=0),
            self.radii[rindex], axis=0)
        # return f*1.649407760419599e-07/(self.starD**2)
        return f*3.299e-7/(self.starD**2)

    def calcFluxWarmAmin(self, waves, r_0, amin1, T_0=1):
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['AstroSil'][rindex]
        graindex = find_nearest_ind(GRAINSIZES, amin1)
        grains = GRAINSIZES[graindex:]/1.0e6
        amin1 /= 1e6
        exponent = -0.5 * ((radii1 - r0) / (self.sigma*r0))**2
        ca = T_0*np.exp(exponent)\
            / (np.power(amin1,3+self.q)-np.power(.001,3+self.q))
        ca1 = np.reshape(ca, (1, ca.size, 1))
        grains1 = np.reshape(grains, (grains.size, 1, 1))
        waves1 = np.broadcast_to(waves, (1, waves.size))
        temps1 = np.broadcast_to(
            grainTemps[:,graindex:],
                (1, grainTemps[:,graindex:].shape[0],
                grainTemps[:,graindex:].shape[1])).T
        emis1 = np.reshape(self.emis['AstroSil'][graindex:],
            (self.emis['AstroSil'][graindex:].shape[0], 1,
            self.emis['AstroSil'][graindex:].shape[1]))
        radii1 = np.reshape(radii1, (1, radii1.size, 1))
        flux = emis1 * ca1 * grains1**(self.q+2) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(
            integrate.simps(flux, grains, axis=0), self.radii[rindex], axis=0)
        # Convert to Jy -> 1e26
        # Convert star distance to meters -> 3.086e16
        # The below constant is: 1e26 / (3.086e16**2) * np.pi
        # computation is faster if we just skip to this.
        return f*3.299e-7/(self.starD**2)

    def calcFluxColdAmin(self, waves, r_0, amin2, T_0=1):
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['DirtyIceAstroSil'][rindex]
        graindex = find_nearest_ind(GRAINSIZES, amin2)
        grains = GRAINSIZES[graindex:]/1.0e6
        amin2 /= 1e6
        exponent = -0.5 * ((radii1 - r0) / (self.sigma*r0))**2
        ca = T_0*np.exp(exponent)\
            / (np.power(amin2,3+self.q)-np.power(.001,3+self.q))
        ca1 = np.reshape(ca, (1, ca.size, 1))
        grains1 = np.reshape(grains, (grains.size, 1, 1))
        waves1 = np.broadcast_to(waves, (1, waves.size))
        temps1 = np.broadcast_to(
            grainTemps[:,graindex:],
                (1, grainTemps[:,graindex:].shape[0],
                grainTemps[:,graindex:].shape[1])).T
        emis1 = np.reshape(self.emis['DirtyIceAstroSil'][graindex:],
            (self.emis['DirtyIceAstroSil'][graindex:].shape[0], 1,
            self.emis['DirtyIceAstroSil'][graindex:].shape[1]))
        radii1 = np.reshape(radii1, (1, radii1.size, 1))
        flux = emis1 * ca1 * grains1**(self.q+2) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(integrate.simps(flux, grains, axis=0),
            self.radii[rindex], axis=0)
        # Convert to Jy -> 1e26
        # Convert star distance to meters -> 3.086e16
        # The below constant is: 1e26 / (3.086e16**2) * np.pi / 2
        # computation is faster if we just skip to this.
        return f*3.299e-7/(self.starD**2)

    def calcFluxWarm(self, waves, r_0, T_0=1):
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['AstroSil'][rindex]
        grains = GRAINSIZES[self.graindex1:]/1.0e6
        bos = self.blowoutSize1/1e6
        exponent = -0.5 * ((radii1 - r0) / (self.sigma*r0))**2
        ca = T_0*np.exp(exponent)\
            / (np.power(bos,3+self.q)-np.power(.001,3+self.q))
        ca1 = np.reshape(ca, (1, ca.size, 1))
        grains1 = np.reshape(grains, (grains.size, 1, 1))
        waves1 = np.broadcast_to(waves, (1, waves.size))
        temps1 = np.broadcast_to(
            grainTemps[:,self.graindex1:],
                (1, grainTemps[:,self.graindex1:].shape[0],
                grainTemps[:,self.graindex1:].shape[1])).T
        emis1 = np.reshape(self.emis['AstroSil'][self.graindex1:],
            (self.emis['AstroSil'][self.graindex1:].shape[0], 1,
            self.emis['AstroSil'][self.graindex1:].shape[1]))
        radii1 = np.reshape(radii1, (1, radii1.size, 1))
        flux = emis1 * ca1 * grains1**(self.q+2) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(integrate.simps(flux, grains, axis=0),
            self.radii[rindex], axis=0)
        return f*3.299e-7/(self.starD**2)

    def calcFluxCold(self, waves, r_0, T_0=1):
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['DirtyIceAstroSil'][rindex]
        grains = self.grains[self.graindex2:]/1.0e6
        bos = self.blowoutSize2/1e6
        exponent = -0.5 * ((radii1 - r0) / (self.sigma*r0))**2
        ca = T_0*np.exp(exponent)\
            / (np.power(bos,3+self.q)-np.power(.001,3+self.q))

        # NumPy Broadcasting enables calculations in C instead of Python
        ca1 = np.reshape(ca, (1, ca.size, 1))
        grains1 = np.reshape(grains, (grains.size, 1, 1))
        waves1 = np.broadcast_to(waves, (1, waves.size))
        temps1 = np.broadcast_to(
            grainTemps[:,self.graindex2:],
                (1, grainTemps[:,self.graindex2:].shape[0],
                grainTemps[:,self.graindex2:].shape[1])).T
        emis1 = np.reshape(self.emis['DirtyIceAstroSil'][self.graindex2:],
            (self.emis['DirtyIceAstroSil'][self.graindex2:].shape[0], 1,
            self.emis['DirtyIceAstroSil'][self.graindex2:].shape[1]))
        radii1 = np.reshape(radii1, (1, radii1.size, 1))
        flux = emis1 * ca1 * grains1**(-1.5) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(integrate.simps(flux, grains, axis=0),
            self.radii[rindex], axis=0)
        return f*3.299e-7/(self.starD**2)

    def deprecated_calcFluxCold(self, waves, r_0, T_0 = 1):
        '''
        This is kept for use as a comparison to the new algorithm.
        '''
        r0 = r_0 *1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['DirtyIceAstroSil'][rindex]
        grains = self.grains[self.graindex2:]/1.0e6
        bos = self.blowoutSize2/1e6
        exponent = -0.5 * ((radii1 - r0) / (self.sigma*r0))**2
        ca = T_0*np.exp(exponent)\
            / (np.power(bos,3+self.q)-np.power(.001,3+self.q))

        # Nested loop in Python
        fw = np.empty(waves.size)
        fr = np.empty(radii1.size)
        for w in range(waves.size):
            for r in range(radii1.size):
                flux = b_nu(waves[w], grainTemps[r, self.graindex2:])
                flux *= self.emis['DirtyIceAstroSil'][self.graindex2:,w] * \
                    (grains**-1.5) * radii1[r] * ca[r]
                fr[r] = integrate.simps(flux, grains)
            fw[w] = integrate.simps(fr, radii1)
        return fw*3.299e-7/(self.starD**2)
