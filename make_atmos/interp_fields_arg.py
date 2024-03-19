#!/usr/bin/env python

import numpy as np
import netCDF4 as nc
import iris, umfile
from um_fileheaders import *
import xesmf as xe
import xarray as xr
import os
import argparse

parser = argparse.ArgumentParser(description='Interpolate fields in UM ancillary file for idealised runs')
parser.add_argument('infile', help='input ancillary file')
parser.add_argument('outfile', help='output ancillary file for idealised run')

args = parser.parse_args()

lsm_pi_file = 'lsm_esm1.5.nc'
lsm_ideal_file = 'lsm_new.nc'

f = nc.Dataset(lsm_pi_file, 'r')
lsm_pi = f.variables['lsm'][:]
f.close()
lsm_pi = lsm_pi.astype('i4')

f = nc.Dataset(lsm_ideal_file, 'r')
lsm_ideal = f.variables['lsm'][:]
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]
f.close()
lsm_ideal = lsm_ideal.astype('i4')

nlat, nlon = lsm_pi.shape

new_ocean = np.logical_and(lsm_pi==1, lsm_ideal==0)
new_land = np.logical_and(lsm_pi==0, lsm_ideal==1)

lon_m, lat_m = np.meshgrid(lon, lat)

# Set up regridder for nearest neighbour interpolations:
old_grid = xr.Dataset(coords={
    "lon": (("y","x"), lon_m),
    "lat": (("y","x"), lat_m)},
    data_vars={
    "mask": (("y","x"), lsm_pi)
    })

new_grid = xr.Dataset(coords={
    "lon": (("y","x"), lon_m),
    "lat": (("y","x"), lat_m)
    })

regridder = xe.Regridder(old_grid, new_grid, method='nearest_s2d', periodic=True)

os.system('cp {} {}'.format(args.infile, args.outfile))

# Coordinates of a point in the Pacific Ocean:
oc_y = 50
oc_x = 118

fin = umfile.UMFile(args.infile, 'r')
f = umfile.UMFile(args.outfile, 'r+')
nvars = f.fixhd[FH_LookupSize2]

for k in range(nvars):
    ilookup = fin.ilookup[k]
    lbegin = ilookup[LBEGIN] 
    if lbegin == -99:
        break
    a = fin.readfld(k)
    a_new = regridder(a)
    a[new_land] = a_new[new_land]
    a[new_ocean] = a[oc_y, oc_x]
    f.writefld(a[:], k)
f.close()
