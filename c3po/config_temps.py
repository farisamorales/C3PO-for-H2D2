import numpy as np
from scipy import integrate
from scipy import constants
import config_sed as conf_s

def energin(star_t, star_l, grainComp):
    '''
    Calculate the energy entering the grain coming from the star.
        Returns shape:
            (N_radii, N_grains) ~ (1000, 121)
    '''
    radii = conf_s.TEMPS_RADII * 1.4959787066e11 # AU -> m
    star_r = np.sqrt(star_l*3.828e26/4/np.pi/constants.sigma/(star_t**4))
    f = conf_s.b_lam(conf_s.WAVELENGTHS, star_t) * conf_s.EMISSIVITIES_TOTAL[grainComp]
    e_in = integrate.simps(f, conf_s.WAVELENGTHS/1e6)
    return (star_r/radii[:, np.newaxis])**2 * e_in

def energout(temps, grainComp):
    '''
    Calculate the energy of dust grains at each radial location.
    Returns an array of shape:
        (N_grains, N_temps) ~ (121, arbitrarily 100)
    '''
    emis = conf_s.EMISSIVITIES_TOTAL[grainComp][:, np.newaxis, :]
    f = conf_s.b_lam(conf_s.WAVELENGTHS, temps[:, np.newaxis]) * emis
    e_out = integrate.simps(f, conf_s.WAVELENGTHS/1e6) * 4
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
    grainTemps = np.empty((conf_s.TEMPS_RADII.size, conf_s.GRAINSIZES.size))
    for r in range(conf_s.TEMPS_RADII.size):
        for a in range(conf_s.GRAINSIZES.size):
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
        plt.loglog(conf_s.GRAINSIZES, gTemps[i*100], label=conf_s.TEMPS_RADII[i*100])
    plt.title('GrainTemps vs GrainSize')
    plt.ylabel('GrainTemp (K)')
    plt.xlabel('GrainSize (um)')
    plt.legend(title='Radial Location (AU)')
    plt.show()

    for i in range(10):
        plt.loglog(conf_s.TEMPS_RADII, gTemps[:, i*10], label=conf_s.GRAINSIZES[i*10])
    plt.title('GrainTemps vs Radii')
    plt.ylabel('GrainTemps (K)')
    plt.xlabel('Radii (AU)')
    plt.legend(title='GrainSize (um)')
    plt.show()
