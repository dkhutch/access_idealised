import numpy as np
import netCDF4 as nc
import os
from scipy.interpolate import griddata

infile = 'basin_mask_orig.nc'
outfile = 'basin_mask.nc'
old_topogfile = 'topog_orig.nc'
new_topogfile = 'topog.nc'
oceangrid = 'ocean_hgrid.nc'

f = nc.Dataset(old_topogfile,'r')
depth_old = f.variables['depth'][:]
f.close()

f = nc.Dataset(new_topogfile,'r')
depth_new = f.variables['depth'][:]
f.close()

ny, nx = depth_old.shape

f = nc.Dataset(oceangrid,'r')
x = f.variables['x'][:]
y = f.variables['y'][:]
f.close()

geolon = x[1::2, 1::2]
geolat = y[1::2, 1::2]

oc_old = depth_old > 0
oc_new = depth_new > 0

xmask = geolon[oc_old].ravel()
ymask = geolat[oc_old].ravel()

cmd = f'cp {infile} {outfile}'
os.system(cmd)

f = nc.Dataset(outfile,'r+')
if 'history' in f.ncattrs(): 
    hist = f.history[:]
else:
    hist = ''
hist = 'basin_mask_fix.py \n' + hist
f.setncattr('history', hist)

varlist = ['BASIN_MASK']

for vname in varlist:
    vv = f.variables[vname]
    data = vv[:]
    dslice = data[0,:,:]

    dmask = dslice[oc_old].ravel()
    d_itp = griddata((xmask, ymask), dmask, (geolon, geolat), method='nearest')

    dslice.mask[oc_new] = False
    dslice.mask[~oc_new] = True
    dslice[oc_new] = d_itp[oc_new]
    data[0,:,:] = dslice
    vv[:] = data[:]
f.close()