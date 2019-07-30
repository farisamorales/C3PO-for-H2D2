import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def weighted_average(arr, weights, axis):
    wa = arr*weights
    return np.nansum(wa, axis=axis) / np.nansum(weights, axis=axis)

def path(path_string):
    return os.sep.join( path_string.split('/') )

direc = path('Results/FluxArrays/Both Wander/')

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
spitzFluxes = np.zeros((len(starNames), spitzWaves.size)) + np.nan
spitzSigmas = np.zeros((len(starNames), spitzWaves.size)) + np.nan
warmFluxes  = np.zeros((len(starNames), spitzWaves.size)) + np.nan
coldFluxes  = np.zeros((len(starNames), spitzWaves.size)) + np.nan
ngFluxes    = np.zeros((len(starNames), spitzWaves.size)) + np.nan

for n, name in enumerate(starNames):
    df = pd.read_csv(direc+'{}.csv'.format(name))

    for i in range(df.loc[:, 'wavelength'].size):
        if np.any(df.loc[i, 'wavelength'] == spitzWaves):
            ind = np.nonzero(spitzWaves == df.loc[i, 'wavelength'])
            spitzFluxes[n][ind] = df.loc[i, 'flux']
            spitzSigmas[n][ind] = df.loc[i, 'error']
            warmFluxes[n][ind]  = df.loc[i, 'warm_ideal']
            coldFluxes[n][ind]  = df.loc[i, 'cold_ideal']
            ngFluxes[n][ind]    = df.loc[i, 'star_ideal']


# Mask the arrays according to the nans in the spitzer data
warmFluxes = np.ma.array(warmFluxes, mask=np.isnan(spitzFluxes))
coldFluxes = np.ma.array(coldFluxes, mask=np.isnan(spitzFluxes))
ngFluxes   = np.ma.array(ngFluxes, mask=np.isnan(spitzFluxes))
# Spitzer has to be last
spitzFluxes = np.ma.array(spitzFluxes, mask=np.isnan(spitzFluxes))
spitzSigmas = np.ma.array(spitzSigmas, mask=np.isnan(spitzSigmas))
spitzFluxes = warmFluxes + coldFluxes + ngFluxes

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
        # (spitzFluxes-ngFluxes),
        weights=weights, axis=0)

    spitzStarCold = np.ma.average(
        np.abs(spitzFluxes-ngFluxes-coldFluxes),
        # (spitzFluxes-ngFluxes-coldFluxes),
        weights=weights, axis=0)

    spitzStarColdWarm = np.ma.average(
        np.abs(spitzFluxes-ngFluxes-coldFluxes-warmFluxes),
        # (spitzFluxes-ngFluxes-coldFluxes-warmFluxes),
        weights=weights, axis=0)

else:
    spitzStar = np.ma.average(
        # np.abs(spitzFluxes-ngFluxes),
        (spitzFluxes-ngFluxes),
        weights=weights, axis=0)

    spitzStarCold = np.ma.average(
        # np.abs(spitzFluxes-ngFluxes-coldFluxes),
        (spitzFluxes-ngFluxes-coldFluxes),
        weights=weights, axis=0)

    spitzStarColdWarm = np.ma.average(
        # np.abs(spitzFluxes-ngFluxes-coldFluxes-warmFluxes),
        (spitzFluxes-ngFluxes-coldFluxes-warmFluxes),
        weights=weights, axis=0)

# Count the modules per wavelength
N_measures = np.sum(np.isfinite(spitzFluxes), axis=0)


fig = plt.figure()
plt.title('Spitzer Stacking: 49 stars')

fig.add_subplot(221)
plt.plot(spitzWaves, spitz_stacked, label='Spitzer')
plt.semilogx()
plt.semilogy()
plt.legend()

fig.add_subplot(222)
plt.plot(spitzWaves, spitzStar, label='Substellar')
plt.semilogx()
plt.semilogy()
plt.legend()

fig.add_subplot(223)
plt.plot(spitzWaves, spitzStarCold, label='Substellarcold')
plt.semilogx()
plt.semilogy()
plt.legend()

fig.add_subplot(224)
plt.plot(spitzWaves, spitzStarColdWarm, label='Substellarcoldwarm')
plt.semilogx()
plt.semilogy()
plt.legend()


plt.show()
# plt.savefig('Using Nonabsolute Values')
