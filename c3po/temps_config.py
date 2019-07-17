import numpy as np
from scipy import integrate
from scipy import constants
import sed_config as sconf

def energin(star_t, star_l, grainComp):
    '''
    Calculate the energy entering the grain coming from the star.
        Returns shape:
            (N_radii, N_grains) ~ (1000, 121)
    '''
    radii = sconf.TEMPS_RADII * 1.4959787066e11 # AU -> m
    star_r = np.sqrt(star_l*3.828e26/4/np.pi/constants.sigma/(star_t**4))
    f = sconf.b_lam(sconf.WAVELENGTHS, star_t) * sconf.EMISSIVITIES_TOTAL[grainComp]
    e_in = integrate.simps(f, sconf.WAVELENGTHS/1e6)
    return (star_r/radii[:, np.newaxis])**2 * e_in

def energout(temps, grainComp):
    '''
    Calculate the energy of dust grains at each radial location.
    Returns an array of shape:
        (N_grains, N_temps) ~ (121, arbitrarily 100)
    '''
    emis = sconf.EMISSIVITIES_TOTAL[grainComp][:, np.newaxis, :]
    f = sconf.b_lam(sconf.WAVELENGTHS, temps[:, np.newaxis]) * emis
    e_out = integrate.simps(f, sconf.WAVELENGTHS/1e6) * 4
    return e_out

def calc_temps(star_t, star_l, grainComp):
    '''
    Interpolate to the star temperature given the properties of the star and the
    grain composition.
    Returns an array of shape:
        (N_radii, N_grains) ~ (1000, 121)
    '''
    temps = np.logspace(-1, 4, 100) # 0.1 - 10k
    energy_in = energin(star_t, star_l, grainComp)
    energy_out = energout(temps, grainComp)
    grainTemps = np.empty((sconf.TEMPS_RADII.size, sconf.GRAINSIZES.size))
    for r in range(sconf.TEMPS_RADII.size):
        for a in range(sconf.GRAINSIZES.size):
            grainTemps[r, a] = np.interp(np.log(energy_in[r, a]),
                np.log(energy_out[a]),np.log(temps))
    return np.exp(grainTemps)

if __name__ == '__main__':
    '''
    Test the functionality of the grain temps calculator. Shows plots of the
    grain temperatures vs grainsize and grain temperatures vs radial location.
    '''
    starT = 5800
    starL = 1
    grainComp = 'DirtyIceAstroSil'
    gTemps = calc_temps(starT, starL, grainComp)

    import matplotlib.pyplot as plt
    for i in range(10):
        plt.loglog(sconf.GRAINSIZES, gTemps[i*100], label=sconf.TEMPS_RADII[i*100])
    plt.title('GrainTemps vs GrainSize')
    plt.ylabel('GrainTemp (K)')
    plt.xlabel('GrainSize (um)')
    plt.legend(title='Radial Location (AU)')
    plt.show()

    for i in range(10):
        plt.loglog(sconf.TEMPS_RADII, gTemps[:, i*10], label=sconf.GRAINSIZES[i*10])
    plt.title('GrainTemps vs Radii')
    plt.ylabel('GrainTemps (K)')
    plt.xlabel('Radii (AU)')
    plt.legend(title='GrainSize (um)')
    plt.show()
