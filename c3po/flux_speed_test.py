
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import ascii
import config_sed as conf_s
import config_temps as conf_t

starName = 'HD 105'
fname = "Data/StarFiles/HD 105_stitched.txt"
grainComps = ['AstroSil', 'DirtyIceAstroSil']
innerGrain = 'AstroSil'
outerGrain = 'DirtyIceAstroSil'
densities = { 'DirtyIce': 1.07, 'AstroSil': 2.7, 'DirtyIceAstroSil': 1.885,
    'WaterIce': 1}
data = conf_s.read_star_file(fname)
starD = data["DIST_pc"]
starL = data["starL"]
starT = data["TEMP"]
starM = data["starM"]
blowoutSize1 = conf_s.blowout_size(densities[innerGrain], starL, starM)
blowoutSize2 = conf_s.blowout_size(densities[outerGrain], starL, starM)
graindex1 = conf_s.find_nearest_ind(conf_s.GRAINSIZES, blowoutSize1)
graindex2 = conf_s.find_nearest_ind(conf_s.GRAINSIZES, blowoutSize2)

grainTemps = dict()
print('loading grains')
for grain in grainComps:
    try:
        grainTemps[grain] = np.load(
            conf_s.GRAIN_TEMPS_DIR+"%s_%s.npy"%(starName, grain))
    except:
        grainTemps[grain] = conf_t.calc_temps(starT, starL, grain)
        np.save(conf_s.GRAIN_TEMPS_DIR+"%s_%s.npy"%(starName, grain),
            grainTemps[grain], allow_pickle=False)

star2 = conf_s.Star(starName, starD, starL, starT, starM, grainTemps)

NLoops = 10

before2 = time.perf_counter()
for _ in range(NLoops):
    y2 = star2.high_res_cold(conf_s.WAVELENGTHS, 10)
after2 = time.perf_counter()
before1 = time.perf_counter()
for _ in range(NLoops):
    y1 = star2.deprecated_calcFluxCold(conf_s.WAVELENGTHS, 10)
after1 = time.perf_counter()
print("--------------------------------------------------")
print(f"New function call time, using NumPy broadcasting:\n{(after2-before2)/NLoops} seconds\n")
print(f"Deprecated function call time, using Python loops:\n{(after1-before1)/NLoops} seconds\n")
print(f"Speed increase per call: {(after1-before1)/(after2-before2)}")
print("--------------------------------------------------")

plt.plot(conf_s.WAVELENGTHS, y1, label="Old Function", lw=5)
plt.plot(conf_s.WAVELENGTHS, y2, ls='dashdot', lw=5, label="New Function")
plt.xlabel(r'$\lambda$ ($\mu m$)')
plt.ylabel(r'$F_{\nu}$ ($Jy$)')
plt.title("A Comparison of Flux Calculating Algorithms")
plt.semilogx()
plt.semilogy()
plt.xlim(5, 1000)
plt.ylim(1e-5, 1e2)
plt.legend(loc='lower right')
# plt.savefig('Fast and Furious')
plt.show()
