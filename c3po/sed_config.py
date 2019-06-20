# Functions for SED fitting
import os
import numpy as np
from scipy import integrate
from scipy import constants
from astropy.io import ascii

'''
One thing I noticed is that the temperatures calculator does not work if it is
done using b_nu instead of b_lam. There's likely something that I'm missing,
having to do with the fact that emissivity is calculated with respect to the
wavelength instead of frequency? Not sure exactly but I'm inclined to also now
question the accuracy of the black body model using b_nu instead of b_lam
during fitting.
'''

# File Paths
ARR_DIR = os.sep.join('Data/Arrays/'.split('/'))
STAR_FILES  = os.sep.join('Data/StarFiles/'.split('/'))
KURUCZ      = os.sep.join('Data/StellarModels/kurucz/'.split('/'))
NEXTGEN     = os.sep.join('Data/StellarModels/nextgen/'.split('/'))
RES_DIR     = os.sep.join('Results/'.split('/'))
INTERPS_DIR = os.sep.join('Data/Arrays/InterpGrainTemps/'.split('/'))
FILTERS_DIR = os.sep.join('Data/FilterResponse/'.split('/'))
GRAIN_TEMPS_DIR = os.sep.join('Data/GrainTemps/'.split('/'))

# Frequently used arrays
GRAINSIZES = np.loadtxt(ARR_DIR + 'GrainSizes.dat')
WAVELENGTHS = np.loadtxt(ARR_DIR + 'Wavelengths.dat')
STAR_TEMPS = np.linspace(2000, 15000, 14)
DISK_RADII = np.logspace(-1, 3, 121)
WAVES = np.logspace(-3, 3, 1000)

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
def ascii_read(starFile):
    '''
    Returns a dictionary with all of the comments included as keys
    '''
    x = ascii.read(starFile)
    data = {}
    for column in x.itercols():
        data[column.name] = column
    for key in x.meta['keywords'].keys():
        data[key] = x.meta['keywords'][key]['value']
    return data
# Convolution function for calibrating data
def convolve(filterwaves, filterresponse, datawaves, dataflux):
    top = integrate.simps(
      np.interp(filterwaves, datawaves, dataflux)*filterresponse,
        filterwaves)
    bottom = integrate.simps(filterresponse, filterwaves)
    return top/bottom

# Load the radial locations as seen by Herschel
def load_spatial_radii():
    spatialRadii = []
    rp_fname = os.sep.join('Data/spatialRadii.csv'.split('/'))
    with open(rp_fname, 'r') as f:
        text = f.readlines()
    for te in text[1:]:
        line = te.split(',')
        radial_name = line[0].strip()
        radial_hipname = line[1].strip()
        radial_loc = line[2].strip()
        radial_err = line[3].strip()
        radial_loc = np.nan if radial_loc == '-' else float(radial_loc)
        radial_err = np.nan if radial_err == '-' else float(radial_err)
        spatialRadii.append( [radial_name, radial_hipname, radial_loc, radial_err] )
    return spatialRadii

# Normalize the black body
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
    grainDensity = grainDensity / 1000 / (0.01**3) # Density converted to SI
    nume = 3 * starL * 3.826e26 * qRad
    deno = 8 * np.pi*starM* 1.989e30 *constants.G*constants.c*grainDensity
    return nume/deno * 1e6 # microns

# Separates data according to instrument. Input: dict of IPAC Tables
def sort_by_instrument(data):
    wa = 'wavelength'
    fx = 'flux'
    er = 'error'
    unique_insts = list()
    for i in range(data['instrument'].size):
        if not data['instrument'][i] in unique_insts:
            unique_insts.append(data['instrument'][i])
    separatedStarData = dict()
    for inst in unique_insts:
        index = np.where(data['instrument']==inst)
        separatedStarData[inst] = dict()
        separatedStarData[inst][wa] = np.array(data[wa][index])
        separatedStarData[inst][fx] = np.array(data[fx][index])
        separatedStarData[inst][er] = np.array(data[er][index])
    return unique_insts, separatedStarData

# Create 1000 radii arrays for given star temp. Used in realistic fitting.
# Takes dictionary of grain temps and list of grain comps.
def interpTemps(starTemp, oldGrainTemps, grainComps):
    STAR_TEMPS = np.linspace(2000, 15000, 14)
    DISK_RADII = np.logspace(-1, 3, 121)
    radii = np.logspace(-1, 3, 1000)
    GRAINSIZES = np.loadtxt(ARR_DIR+'GrainSizes.dat')
    for grainComp in grainComps:
        abr = ''
        for letter in grainComp:
            if letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                abr += letter
        starIndices = np.where(np.logical_and(
            STAR_TEMPS<starTemp+3000,
            STAR_TEMPS>starTemp-3000
            ))
        newStarTempGrainTemps = np.empty((
            DISK_RADII.size,
            GRAINSIZES.size
            ))
        for r in range(DISK_RADII.size):
            for g in range(GRAINSIZES.size):
                newStarTempGrainTemps[r][g] = np.interp(
                    starTemp,
                    STAR_TEMPS[starIndices],
                    oldGrainTemps[grainComp][starIndices][:,r][:,g]
                    )
        newGrainTemps = np.empty((radii.size,GRAINSIZES.size))
        for r  in range(radii.size):
            for g in range(GRAINSIZES.size):
                newGrainTemps[r][g] = np.interp(
                    radii[r],
                    DISK_RADII,
                    newStarTempGrainTemps[:,g]
                    )
        np.save(ARR_DIR+'InterpGrainTemps/'+'%.0fK_%s.npy'%
            (starTemp,grainComp), newGrainTemps, allow_pickle=False)

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

    def calcFluxWarmMinGrains(self, waves, r0, bos, T_0=1):
        sigma = 0.10
        r0 *= 1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['AstroSil'][rindex]
        graindex = find_nearest_ind(GRAINSIZES, bos)
        grains = GRAINSIZES[graindex:self.graindex1]/1.0e6
        bos /= 1e6
        q = -3.5
        exponent = -0.5 * ((radii1 - r0) / (sigma*r0))**2
        ca = T_0*np.exp(exponent)*np.abs(3+q) \
            / (np.pi*(np.power(bos,3+q)-np.power(self.blowoutSize1,3+q)))
        ca *= 1e6
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
        flux = emis1 * ca1 * grains1**(-1.5) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(integrate.simps(flux, grains, axis=0),
            self.radii[rindex], axis=0)
        return f*1.649407760419599e-07/(self.starD**2)

    def calcFluxBoSWarm(self, waves, r0, bos, T_0=1):
        sigma = 0.10
        r0 *= 1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['AstroSil'][rindex]
        graindex = find_nearest_ind(GRAINSIZES, bos)
        grains = GRAINSIZES[graindex:]/1.0e6
        bos /= 1e6
        q = -3.5
        exponent = -0.5 * ((radii1 - r0) / (sigma*r0))**2
        ca = T_0*np.exp(exponent)*np.abs(3+q) \
            / (np.pi*(np.power(bos,3+q)-np.power(.001,3+q)))
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
        flux = emis1 * ca1 * grains1**(-1.5) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(
            integrate.simps(flux, grains, axis=0), self.radii[rindex], axis=0)
        # Convert to Jy -> 1e26
        # Convert star distance to meters -> 3.086e16
        # The below constant is: 1e26 / (3.086e16**2) * np.pi / 2
        # computation is faster if we just skip to this.
        return f*1.649407760419599e-07/(self.starD**2)

    def calcFluxBoSCold(self, waves, r0, bos, T_0=1):
        sigma = 0.10
        r0 *= 1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['DirtyIceAstroSil'][rindex]
        graindex = find_nearest_ind(GRAINSIZES, bos)
        grains = GRAINSIZES[graindex:]/1.0e6
        bos /= 1e6
        q = -3.5
        exponent = -0.5 * ((radii1 - r0) / (sigma*r0))**2
        ca = T_0*np.exp(exponent)*np.abs(3+q) \
            / (np.pi*(np.power(bos,3+q)-np.power(.001,3+q)))
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
        flux = emis1 * ca1 * grains1**(-1.5) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(integrate.simps(flux, grains, axis=0),
            self.radii[rindex], axis=0)
        # Convert to Jy -> 1e26
        # Convert star distance to meters -> 3.086e16
        # The below constant is: 1e26 / (3.086e16**2) * np.pi / 2
        # computation is faster if we just skip to this.
        return f*1.649407760419599e-07/(self.starD**2)

    def calcFluxWarm(self, waves, r0, T_0=1):
        sigma = 0.10
        r0 *= 1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['AstroSil'][rindex]
        grains = GRAINSIZES[self.graindex1:]/1.0e6
        bos = self.blowoutSize1/1e6
        q = -3.5
        exponent = -0.5 * ((radii1 - r0) / (sigma*r0))**2
        ca = T_0*np.exp(exponent)*np.abs(3+q) \
            / (np.pi*(np.power(bos,3+q)-np.power(.001,3+q)))
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
        flux = emis1 * ca1 * grains1**(-1.5) * radii1 * b_nu(waves1, temps1)
        f = integrate.simps(integrate.simps(flux, grains, axis=0),
            self.radii[rindex], axis=0)
        return f*1.649407760419599e-07/(self.starD**2)

    def calcFluxCold(self, waves, r0, T_0=1):
        sigma = 0.10
        r0 *= 1.4959787066e11
        rindex = np.where(np.logical_and(self.radii<1.4*r0,
            self.radii>0.6*r0))[0]
        radii1 = self.radii[rindex]
        grainTemps = self.grainTemps['DirtyIceAstroSil'][rindex]
        grains = self.grains[self.graindex2:]/1.0e6
        bos = self.blowoutSize2/1e6
        q = -3.5
        exponent = -0.5 * ((radii1 - r0) / (sigma*r0))**2
        ca = T_0*np.exp(exponent)*np.abs(3+q) \
            / (np.pi*(np.power(bos,3+q)-np.power(.001,3+q)))
        ca *= 1e6
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
        return f*1.649407760419599e-07/(self.starD**2)
