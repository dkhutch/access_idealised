import numpy as np
import netCDF4 as nc
import os
from scipy.interpolate import griddata

infile = 'csiro_bgc_sediment.res.orig.nc'
outfile = 'csiro_bgc_sediment.res.nc'
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
nz = 50

f = nc.Dataset(oceangrid,'r')
x = f.variables['x'][:]
y = f.variables['y'][:]
f.close()

geolon = x[1::2, 1::2]
geolat = y[1::2, 1::2]

oc_old = depth_old > 0
oc_new = depth_new > 0

cmd = f'cp {infile} {outfile}'
os.system(cmd)
cmd = f'ncatted -a checksum,.,d,, {outfile}'
os.system(cmd)


f = nc.Dataset(outfile,'r+')
if 'history' in f.ncattrs(): 
    hist = f.history[:]
else:
    hist = ''
hist = 'bgc_sediment.py \n' + hist
f.setncattr('history', hist)

varlist = ['caco3_sediment',
           'det_sediment']

for vname in varlist:
    vv = f.variables[vname]
    data = vv[:]
    dslice = data[0,:,:]
    mask = dslice != 0

    dmask = dslice[mask].ravel()
    xmask = geolon[mask].ravel()
    ymask = geolat[mask].ravel()
    d_itp = griddata((xmask, ymask), dmask, (geolon, geolat), method='nearest')

    dslice[oc_new] = d_itp[oc_new]
    dslice[~oc_new] = 0.
    data[0,:,:] = dslice
    vv[:] = data[:]
f.close()