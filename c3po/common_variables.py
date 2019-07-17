# Compositions of the grains.
grainComps = ['AstroSil', 'DirtyIceAstroSil']

# Grain densities. Important for calculating the blowout size of the grains.
densities = {
    # For the density of the grains:
    # Mixture of H20:NH3 is 3:1 and 10% vlfr of aC to make dirty ice.
    # then dirty ice is 50% vlfr of AstroSil to make DIAS
    # Assuming:
    # 1.0 g/cm^3 for H20
    # 0.817 g/cm^3 for NH3
    # 3.0 g/cm^3 for AstroSil
    'DirtyIce': 1.07,
    'AstroSil': 2.7,
    'DirtyIceAstroSil': 1.885, # 1:1 DirtyIce:AstroSil
    'WaterIce': 1,
    # 'AstroSil': 3.0,           # Correlates to DIAS -> 2.034 g/cm^3
    # 'DirtyIceAstroSil': 2.034, # Correlates to AS -> 3.0 g/cm^3
    }

# Colors for plotting accessed by instrument name
plotColors = {
    'B_Tycho': 'brown', 'V_Tycho': 'orange',

    '2MASSJ': 'red', '2MASSH': 'red', '2MASSK': 'red',

    'WISE1': 'blue', 'WISE2': 'blue', 'WISE3': 'blue', 'WISE4': 'blue',
    'AllSkyW1': 'blue', 'AllSkyW2': 'blue', 'AllSkyW3': 'blue',
    'AllSkyW4': 'blue', 'AllWiseW1': 'blue', 'AllWiseW2': 'blue',
    'AllWiseW3': 'blue', 'AllWiseW4': 'blue',
    'MIPS24': 'green', 'MIPS70': 'green',
    'MIPS24UL':'green', 'MIPS70UL':'green',

    'HerschelPACS70': 'purple', 'HerschelPACS100': 'purple',
    'HerschelPACS160': 'purple', 'HerschelPACS70UL': 'purple',
    'HerschelPACS100UL': 'purple', 'HerschelPACS160UL': 'purple',
    }

# Saturation limits for each instrument, in Jy
satLims = {
    '2MASSJ': 10.057, '2MASSH': 10.24, '2MASSK': 10.566,
    'WISE1':0.18, 'WISE2': 0.36, 'WISE3': 0.88, 'WISE4': 12.0,
    'HerschelPACS70': 220., 'HerschelPACS100': 510., 'HerschelPACS160':1125.,
    'MIPS24': float('inf'), 'MIPS70': float('inf')
    }

# Conversion factors
moon_mass = 7.34767309e22 # kg
AU_m = 1.4959787066e11 # Astronomical Units -> Meters
PC_m = 3.086e33  # Parsecs -> Meters
solL_m = 3.828e26 # Solar Luminosities -> W/m^2/s/m
