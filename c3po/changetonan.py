import os
import numpy as np
spatialRadii = []
rp_fname = os.sep.join('Data/spatialRadii.csv'.split('/'))
with open(rp_fname, 'r') as f:
    text = f.readlines()
for te in text[1:]:
    line = te.split(',')
    radial_name = line[0].strip()
    radial_hipname = line[1].strip()
    radial_loc = line[2].strip()
    radial_err = line[3].strip()
    radial_loc = np.nan if radial_loc == '-' else float(radial_loc)
    radial_err = np.nan if radial_err == '-' else float(radial_err)
    spatialRadii.append( [radial_name, radial_hipname, radial_loc, radial_err] )

with open(rp_fname, 'w+') as f:
    f.write('name,hipname,radius,radius_uncertainty\n')
    for rad in spatialRadii:
        rad = [str(r) for r in rad]
        f.write(','.join(rad)+'\n')
