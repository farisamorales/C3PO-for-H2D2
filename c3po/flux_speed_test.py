
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import ascii
import sed_config
import temps_config

starName = 'HD 105'
fname = "C:\\Users\\justi\\C3PO-for-H2D2\\c3po\\Data\\StarFiles\\"
fname += starName +"_stitched.txt"
grainComps = ['AstroSil', 'DirtyIceAstroSil']
innerGrain = 'AstroSil'
outerGrain = 'DirtyIceAstroSil'
densities = { 'DirtyIce': 1.07, 'AstroSil': 2.7, 'DirtyIceAstroSil': 1.885,
    'WaterIce': 1}
data = sed_config.ascii_read(fname)
starD = data["DIST_pc"]
starL = data["starL"]
starT = data["TEMP"]
starM = data["starM"]
blowoutSize1 = sed_config.blowout_size(densities[innerGrain], starL, starM)
blowoutSize2 = sed_config.blowout_size(densities[outerGrain], starL, starM)
graindex1 = sed_config.find_nearest_ind(sed_config.GRAINSIZES, blowoutSize1)
graindex2 = sed_config.find_nearest_ind(sed_config.GRAINSIZES, blowoutSize2)
TOTAL_EMISSIVITIES = dict()
for grain in grainComps:
    TOTAL_EMISSIVITIES[grain] = np.load(
        sed_config.ARR_DIR+grain+'Emissivities.npy')
star2 = sed_config.Star(starD, starL, starT, None, blowoutSize1,
    blowoutSize2, TOTAL_EMISSIVITIES, sed_config.GRAINSIZES, graindex1,
    graindex2)
grainTemps = dict()
for grain in grainComps:
    try:
        grainTemps[grain] = np.load(
            sed_config.GRAIN_TEMPS_DIR+"%s_%s.npy"%(starName, grain))
    except:
        grainTemps[grain] = temps_config.calcTemps(star2, grain)
        np.save(sed_config.GRAIN_TEMPS_DIR+"%s_%s.npy"%(starName, grain),
            grainTemps[grain], allow_pickle=False)
star2 = sed_config.Star(starD, starL, starT, grainTemps, blowoutSize1,
    blowoutSize2, TOTAL_EMISSIVITIES, sed_config.GRAINSIZES, graindex1,
    graindex2)

before1 = time.perf_counter()
y1 = star2.deprecated_calcFluxCold(sed_config.WAVELENGTHS, 10)
after1 = time.perf_counter()
before2 = time.perf_counter()
y2 = star2.calcFluxCold(sed_config.WAVELENGTHS, 10)
after2 = time.perf_counter()

print("--------------------------------------------------")
print(f"Deprecated function call time, using Python loops:\n{after1-before1} seconds\n")
print(f"New function call time, using NumPy broadcasting:\n{after2-before2} seconds\n")
print(f"Speed increase per call: {(after1-before1)/(after2-before2)}")
print("--------------------------------------------------")
plt.plot(sed_config.WAVELENGTHS, y1, label="Old Function", lw=5)
plt.plot(sed_config.WAVELENGTHS, y2, ls='dashdot', lw=5, label="New Function")
plt.xlabel(r'$\lambda$ ($\mu m$)')
plt.ylabel(r'$F_{\nu}$ ($Jy$)')
plt.title("A Comparison of Flux Calculating Algorithms.")
plt.semilogx()
plt.semilogy()
plt.xlim(0.1, 1000)
plt.legend(loc='lower right')
plt.show()
