'''
Determine the general options for the fitting routine. Can later be ignored in
the github repository, so that the client's code doesn't change when they
pull any updates to their machine.
'''

import os
import sed_config as sconf

# Grab star names
starNames = []
with open(os.sep.join("Data/starNames.txt".split('/')), 'r') as f:
    names = f.readlines()
for name in names:
    name = name.strip()
    starNames.append(name)

################################################################################
#                               Begin Options
################################################################################
# starNames = starNames[41:]
# starNames = ['HD 106252']
# starNames = ['HD 204961', 'HD 105211', 'HD 106252']
# starNames = ['GJ 285', 'BD-19 1062']
# starNames = ['HD 105']
# starNames = ['HD 71722']
# starNames = ['HD 115600']
# starNames = ['HD 166', 'HD 10939', 'HD 71722']
# starNames = ['HD 15115']

# Set the upper and lower bounds for the fitting
# Upper bound for the cold belt radius:
maxRad = 500
# Sets the bound multiplier for the blowout size (these are default values):
# Inner belt
bosLowInner = 10
bosHighInner = 3
# Outer belt
bosHighOuter = 3
bosLowOuter = 10
# Sets the bound multiplier for the belt radius:
beltBound = 6
# Scalar multiplier for the warm belt minimum blowout size
bosScalar1 = 1
# Scalar multiplier for the cold belt minimum blowout size
bosScalar2 = 1


# Toggle plots showing after each fit:
showFigure = 1
# Toggle plots saving into their respective directories after each fit:
saveFigure = 0
# Toggle save into a results table
saveResults = 0
# Toggle discarding measurements beyond the saturation limit:
saturation_limits = 0
# Toggle using the resolved radii for the cold belt initial guess & bounds
useSpatialRadii = 1

# Show the resolved radius on the plot:
showResolved = 1
# Show the minimum blow out size on the plot:
showMinGrain = 1
# Show the IRS variance on the plot
showIRSVariance = 0
# Show the normalized belts
showNormedBelts = 0
# Show the luminosity ratios of the coldbelt/star and warmbelt/star
showLumRatios = 0
# Show the flux ratios of the coldbelt/star and warmbelt/star
showFluxRatios = 0
# Show the mass of each dust component in lunar masses
showDustMass = 0

# Fitting routines (only 1 can be active at once):
# One warm belt with a wandering grain size and one fixed grain size cold belt:
oneWander = 0
# One warm belt and one cold belt, both with wandering grain sizes:
twoWander = 0
# One warm belt and one cold belt (fixed blowout size):
noWander = 0
# Two co-located warm belts (one wandering & one fixed) and one fixed cold belt:
twoWarmBelts = 1

# Use AllSky or AllWise. You can only select ONE at a time. Or neither.
useAllSky = 0
useAllWise = 0

################################################################################
#                                 End Options
################################################################################
# Set the proper directories
if oneWander:
    IMG_DIR     = sconf.PLOTS_DIR + os.sep.join("Warm Wander/".split('/'))
    PARAMS_DIR  = sconf.PARAMS_DIR + os.sep.join("Warm Wander/".split('/'))
    FRATIOS_DIR = sconf.FRATIOS_DIR + os.sep.join("Warm Wander/".split('/'))
elif twoWander:
    IMG_DIR     = sconf.PLOTS_DIR + os.sep.join("Both Wander/".split('/'))
    PARAMS_DIR  = sconf.PARAMS_DIR + os.sep.join("Both Wander/".split('/'))
    FRATIOS_DIR = sconf.FRATIOS_DIR + os.sep.join("Both Wander/".split('/'))
elif noWander:
    IMG_DIR     = sconf.PLOTS_DIR + os.sep.join("No Wander/".split('/'))
    PARAMS_DIR  = sconf.PARAMS_DIR + os.sep.join("No Wander/".split('/'))
    FRATIOS_DIR = sconf.FRATIOS_DIR + os.sep.join("No Wander/".split('/'))
elif twoWarmBelts:
    IMG_DIR     = sconf.PLOTS_DIR + os.sep.join("Two Warm Belts/".split('/'))
    PARAMS_DIR  = sconf.PARAMS_DIR + os.sep.join("Two Warm Belts/".split('/'))
    FRATIOS_DIR = sconf.FRATIOS_DIR + os.sep.join("Two Warm Belts/".split('/'))
