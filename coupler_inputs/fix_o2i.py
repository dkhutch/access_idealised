import numpy as np
import netCDF4 as nc
import os
from scipy.interpolate import griddata

o2ifile = 'o2i_orig.nc'
o2inew = 'o2i.nc'
old_topogfile = '../make_ocean/topog_orig.nc'
new_topogfile = '../make_ocean/topog.nc'
oceangrid = '../make_ocean/ocean_hgrid.nc'

os.system('cp {} {}'.format(o2ifile, o2inew))

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

varlist = ['sst_i','sss_i','ssu_i','ssv_i','sslx_i',
    'ssly_i','pfmice_i','co2_oi','co2fx_oi']

f = nc.Dataset(o2inew, 'r+')
if 'history' in f.ncattrs(): 
    hist = f.history[:]
else:
    hist = ''
hist = 'fix_o2i_mio.py \n ' + hist
f.setncattr('history',hist)
f.setncattr('title', o2inew)

for vname in varlist:
    vv = f.variables[vname]
    data = vv[:]
    dslice = data[0,:,:]
    dmask = dslice[oc_old].ravel()
    d_itp = griddata((xmask, ymask), dmask, (geolon, geolat), method='nearest')

    dslice[oc_new] = d_itp[oc_new]
    dslice[~oc_new] = 0.
    if vname not in ['ssu_i', 'ssv_i']:
        data[0,:,:] = dslice
    vv[:] = data[:]
f.close()