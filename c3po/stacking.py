import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def weighted_average(arr, weights, axis):
    wa = arr*weights
    return np.nansum(wa, axis=axis) / np.nansum(weights, axis=axis)

def to_path(path_string):
    return os.sep.join( path_string.split('/') )

direc = to_path('Results/FluxArrays/Both Wander/')

starNames = []
# only do stars from our sample of 49 solar and A type stars
with open(os.sep.join("Data/starNames.txt".split('/')), 'r') as f:
    names = f.readlines()
for name in names:
    name = name.strip()
    starNames.append(name)

# Do all of the stars that have spitzer data?
# starNames = os.listdir('Data/Spitzer/')

spitzWaves = pd.read_csv('Data/Arrays/MasterSpitzerWaves2.txt',
    index_col=False).loc[:, 'wavelength'].to_numpy()
# spitzWaves = spitzWaves.loc[:, 'wavelength'].to_numpy()
spitzFluxes     = np.full( (len(starNames), spitzWaves.size), np.nan)
spitzSigmas     = np.full( (len(starNames), spitzWaves.size), np.nan)
warmFluxes      = np.full( (len(starNames), spitzWaves.size), np.nan)
coldFluxes      = np.full( (len(starNames), spitzWaves.size), np.nan)
ngFluxes        = np.full( (len(starNames), spitzWaves.size), np.nan)
medWarmFluxes   = np.full( (len(starNames), spitzWaves.size), np.nan)
smallWarmFluxes = np.full( (len(starNames), spitzWaves.size), np.nan)

for n, name in enumerate(starNames):
    df = pd.read_csv(direc+'{}.csv'.format(name))

    for i in range(df.loc[:, 'wavelength'].size):
        if np.any(df.loc[i, 'wavelength'] == spitzWaves):
            ind = np.nonzero(spitzWaves == df.loc[i, 'wavelength'])
            spitzFluxes[n][ind]        = df.loc[i,             'flux']
            spitzSigmas[n][ind]        = df.loc[i,            'error']
            warmFluxes[n][ind]         = df.loc[i,       'warm_ideal']
            coldFluxes[n][ind]         = df.loc[i,       'cold_ideal']
            ngFluxes[n][ind]           = df.loc[i,       'star_ideal']
            medWarmFluxes[n][ind]      = df.loc[i,   'med_warm_ideal']
            smallWarmFluxes[n][ind]    = df.loc[i, 'small_warm_ideal']

# Mask the arrays according to the nans in the spitzer data
spitzFluxes     = np.ma.array(     spitzFluxes, mask=np.isnan(spitzFluxes))
spitzSigmas     = np.ma.array(     spitzSigmas, mask=spitzFluxes.mask)
warmFluxes      = np.ma.array(      warmFluxes, mask=spitzFluxes.mask)
coldFluxes      = np.ma.array(      coldFluxes, mask=spitzFluxes.mask)
ngFluxes        = np.ma.array(        ngFluxes, mask=spitzFluxes.mask)
# These next two arrays will have NaNs independent of spitzFluxes
medWarmFluxes   = np.ma.array(   medWarmFluxes, mask=np.isnan(medWarmFluxes)+spitzFluxes.mask)
smallWarmFluxes = np.ma.array( smallWarmFluxes, mask=np.isnan(smallWarmFluxes)+spitzFluxes.mask)

totalFluxes = warmFluxes + coldFluxes + ngFluxes

# Now, calculate the weighted averages for each method.
'''
1) spitz - star
2) spitz - star - cold
3) spitz - star - cold - warm (medium grains)
'''

weights = spitzSigmas**-2
stacked_error = np.sqrt(1./np.nansum(weights, axis=0))
spitz_stacked = np.ma.average(spitzFluxes, weights=weights, axis=0)

useAbs = 0

if useAbs:
    spitzStar = np.ma.average(
        np.abs(spitzFluxes-ngFluxes),
        weights=weights, axis=0)

    spitzStarCold = np.ma.average(
        np.abs(spitzFluxes-ngFluxes-coldFluxes),
        weights=weights, axis=0)

    spitzStarColdWarm = np.ma.average(
        np.abs(spitzFluxes-ngFluxes-coldFluxes-medWarmFluxes),
        weights=weights, axis=0)

else:
    spitzStar = np.ma.average(
        spitzFluxes-ngFluxes,
        weights=weights, axis=0)

    spitzStarCold = np.ma.average(
        spitzFluxes-ngFluxes-coldFluxes,
        weights=weights, axis=0)

    spitzStarColdWarm = np.ma.average(
        spitzFluxes-ngFluxes-coldFluxes-medWarmFluxes,
        weights=weights, axis=0)

# st_small = np.ma.average(smallWarmFluxes, weights=weights, axis=0)
# plt.plot(spitzWaves, st_small)
# plt.semilogx()
# plt.semilogy()
# plt.show()
# quit()

# Count the modules per wavelength
N_measures = np.sum(np.isfinite(spitzFluxes), axis=0)


fig, axes = plt.subplots(2, 2)#, sharey='row', sharex='col')
fig.suptitle('Spitzer Stacking for 49 sources')

axes[0, 0].plot(spitzWaves, spitz_stacked, label='Spitzer')
axes[0, 0].semilogx()
axes[0, 0].semilogy()
axes[0, 0].set_ylabel(r'$F_{\nu} (Jy)$')
axes[0, 0].get_xaxis().set_visible(False)

axes[0, 1].plot(spitzWaves, spitzStar, label='Spitzer2')
axes[0, 1].semilogx()
axes[0, 1].semilogy()
axes[0, 1].get_xaxis().set_visible(False)

axes[1, 0].plot(spitzWaves, spitzStarCold, label='Substellarcold')
axes[1, 0].set_ylabel(r'$F_{\nu} (Jy)$')
axes[1, 0].set_xlabel(r'$\lambda (\mu m)$')
axes[1, 0].semilogx()
axes[1, 0].semilogy()

axes[1, 1].plot(spitzWaves, spitzStarColdWarm, label='Substellarcoldwarm')
axes[1, 1].set_xlabel(r'$\lambda (\mu m)$')
axes[1, 1].semilogx()
axes[1, 1].semilogy()


# plt.show()
plt.savefig('Stacking Example')
