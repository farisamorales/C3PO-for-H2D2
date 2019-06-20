import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import constants
import sed_config


def energin(waves, emis, star_t, star_l):
    '''
    Calculate the energy entering the grain coming from the star.
        Returns shape:
            (N_radial_locations, N_grainsizes)
    '''
    # Convert to appropriate units
    waves_m = waves/1e6              # um -> m
    radii = np.logspace(-1, 3, 1000) * 1.4959787066e11 # AU -> m

    # Calculate stellar radius from luminosity/temp
    star_r = np.sqrt(star_l*3.828e26/4/np.pi/constants.sigma/(star_t**4))

    # Calculate flux wrt wavelength (one per grainsize)
    f = sed_config.b_lam(waves, star_t) * emis
    flux = integrate.simps(f, waves_m)

    # Calculate energy into each grain per radial location per frequency
    grains = np.broadcast_to(sed_config.GRAINSIZES/1e6,
        (radii.size, sed_config.GRAINSIZES.size))
    radii = np.broadcast_to(radii, (1, radii.size)).T
    return (star_r/radii)**2 * flux

def energout(waves, emis, temps):
    '''
    Calculate the energy of dust grains at each radial location.
    Returns an array of shape:
        (N_grainsizes, N_temps)
    '''
    # Convert units
    waves_m = waves/1e6

    # Broadcast/reshape arrays so that looping is done in C
    temps1 = np.broadcast_to(temps, (1, temps.size)).T
    emis1 = np.reshape(emis, (emis.shape[0], 1, emis.shape[1]))

    # Calculate flux at every temp/wavelength
    f = sed_config.b_lam(waves, temps1) * emis1
    e_out = integrate.simps(f, waves_m) * 4
    return np.broadcast_to(e_out, (1000, e_out.shape[0],
        e_out.shape[1]))

def interp_to_star(energy_in, energy_out, temps):
    '''
    Returns shape:
        (N_radial_locations, N_grainsizes)
    '''
    grainTemps = np.empty((energy_in.shape[0], energy_out.shape[1]))
    for r in range(1000):
        for a in range(sed_config.GRAINSIZES.size):
            grainTemps[r, a] = np.interp(
                    np.log(energy_in[r, a]),
                    np.log(energy_out[r, a]),
                    np.log(temps)
                    )
    return np.exp(grainTemps)


def calcTemps(star_obj, grainComp):
    '''
    This is function that will be used in sed_config to calculate the grain
    temps for each star IF the temps haven't been calculated already.
    '''
    star_t = star_obj.starT
    star_l = star_obj.starL
    emis = star_obj.emis[grainComp]
    temps = np.logspace(-1, 4, 100) # 0.1 - 10k
    # Energy from star
    energy_in = energin(sed_config.WAVELENGTHS, # um
        emis, # scalar multiple
        star_t, # Kelvin
        star_l) # solar luminosities
    # Energy that would come from each grain at these temps
    energy_out = energout(sed_config.WAVELENGTHS,
        sed_config.EMISSIVITIES_TOTAL['AstroSil'],
        temps
        )
    gTemps = interp_to_star(energy_in, energy_out, temps)
    return gTemps


temps = np.logspace(-1, 4, 100)
starT = 9750















if __name__ == '__main__':
    # Energy from star
    energy_in = energin(sed_config.WAVELENGTHS, # um
        sed_config.EMISSIVITIES_TOTAL['AstroSil'], # scalar multiple
        starT, # Kelvin
        5) # solar luminosities
    # print(energy_in.shape)

    # Energy that would come from each grain at these temps
    energy_out = energout(sed_config.WAVELENGTHS,
        sed_config.EMISSIVITIES_TOTAL['AstroSil'],
        temps
        )
    # print( energy_out.shape )

    gTemps = interp_to_star(energy_in, energy_out, temps)

    radii = np.logspace(-1, 3, 1000)

    for i in range(10):
        plt.loglog(sed_config.GRAINSIZES, gTemps[i*100], label=radii[i*100])
    plt.legend()
    plt.show()



    for i in range(10):
        plt.loglog(radii, gTemps[:, i*10], label=sed_config.GRAINSIZES[i*10])
    plt.legend()
    plt.show()
