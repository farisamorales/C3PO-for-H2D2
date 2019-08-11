'''
Determine the general options for the fitting routine. Can later be ignored in
the github repository, so that the client's code doesn't change when they
pull any updates to their machine.
'''


import os
import config_sed as conf_s

# Grab star names
starNames = []
with open(os.sep.join("Data/starNames.txt".split('/')), 'r') as f:
# with open(os.sep.join("Data/totalStarNames.txt".split('/')), 'r') as f:
    names = f.readlines()
for name in names:
    name = name.strip()
    starNames.append(name)

################################################################################
#                               Begin Options
################################################################################

# WISE data is crazy
# HD 191174

# cccc = starNames.index('HD 71722')
# starNames = ['HD 32297']
# starNames = ['HD 71722']
# starNames = ['HD 128207']
# starNames = ['HD 138965']
# starNames = starNames[cccc:]
# starNames = ['BD-19 1062']
# starNames = ['HD 10939']
# starNames = starNames[11:]

# Set the upper and lower bounds for the fitting
# Upper bound for the cold belt radius:
maxRad = 500
# Sets the bound multiplier for the blowout size (these are default values):
# Inner belt
bosLowInner = 0.1
bosHighInner = 3
# Outer belt
bosLowOuter = 0.1
bosHighOuter = 3

# Sets the bound multiplier for the belt radius:
beltBound = 10
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
# Stacking mode just saves more arrays that can be used for stacking analysis
stackSave = 0
# Save the spitzer data to a file
saveSpitzerData = 0
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
twoWander = 1
# One warm belt and one cold belt (fixed blowout size):
noWander = 0
# Two co-located warm belts (one wandering & one fixed) and one fixed cold belt:
twoWarmBelts = 0
# One cold belt with a wandering grainsize
oneBeltWander = 0

# Use AllSky or AllWise. Selecting both uses both, selecting none uses none
useAllSky = 0
useAllWise = 0

# If for some reason you want to recalculate all of the grain temps
recalcTemps = 0

################################################################################
#                                 End Options
################################################################################
# Set the proper directories
if oneWander:
    IMG_DIR     = conf_s.PLOTS_DIR + os.sep.join("Warm Wander/".split('/'))
    PARAMS_DIR  = conf_s.PARAMS_DIR + os.sep.join("Warm Wander/".split('/'))
    FRATIOS_DIR = conf_s.FRATIOS_DIR + os.sep.join("Warm Wander/".split('/'))
elif twoWander:
    IMG_DIR     = conf_s.PLOTS_DIR + os.sep.join("Both Wander/".split('/'))
    PARAMS_DIR  = conf_s.PARAMS_DIR + os.sep.join("Both Wander/".split('/'))
    FRATIOS_DIR = conf_s.FRATIOS_DIR + os.sep.join("Both Wander/".split('/'))
elif noWander:
    IMG_DIR     = conf_s.PLOTS_DIR + os.sep.join("No Wander/".split('/'))
    PARAMS_DIR  = conf_s.PARAMS_DIR + os.sep.join("No Wander/".split('/'))
    FRATIOS_DIR = conf_s.FRATIOS_DIR + os.sep.join("No Wander/".split('/'))
elif twoWarmBelts:
    IMG_DIR     = conf_s.PLOTS_DIR + os.sep.join("Two Warm Belts/".split('/'))
    PARAMS_DIR  = conf_s.PARAMS_DIR + os.sep.join("Two Warm Belts/".split('/'))
    FRATIOS_DIR = conf_s.FRATIOS_DIR + os.sep.join("Two Warm Belts/".split('/'))
elif oneBeltWander:
    IMG_DIR     = conf_s.PLOTS_DIR + os.sep.join("One Cold Wander/".split('/'))
    PARAMS_DIR  = conf_s.PARAMS_DIR + os.sep.join("One Cold Wander/".split('/'))
    FRATIOS_DIR = conf_s.FRATIOS_DIR + os.sep.join("One Cold Wander/".split('/'))
