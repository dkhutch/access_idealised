import numpy as np
import netCDF4 as nc
import os
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

i2afile = 'i2a_orig.nc'
i2anew = 'i2a.nc'
old_topogfile = '../make_ocean/topog_orig.nc'
new_topogfile = '../make_ocean/topog.nc'
oceangrid = '../make_ocean/ocean_hgrid.nc'

os.system('cp {} {}'.format(i2afile, i2anew))

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

varlist = ['isst_ia_degc',
    'uvel_ia',
    'vvel_ia',
    'icecon01',
    'icecon02',
    'icecon03',
    'icecon04',
    'icecon05',
    'snwthk01',
    'snwthk02',
    'snwthk03',
    'snwthk04',
    'snwthk05',
    'icethk01',
    'icethk02',
    'icethk03',
    'icethk04',
    'icethk05',
    'isss_ia',
    'isst_ia',
    'co2_i2',
    'co2fx_i2']

f = nc.Dataset(i2anew, 'r+')
if 'history' in f.ncattrs(): 
    hist = f.history[:]
else:
    hist = ''
hist = 'fix_i2a_mio.py \n ' + hist
f.setncattr('history',hist)
f.setncattr('title',i2anew)

for vname in varlist:
    vv = f.variables[vname]
    data = vv[:]
    dslice = data[:]
    dmask = dslice[oc_old].ravel()
    d_itp = griddata((xmask, ymask), dmask, (geolon, geolat), method='nearest')

    dslice[oc_new] = d_itp[oc_new]
    dslice[~oc_new] = 0.
    if vname not in ['uvel_ia', 'vvel_ia']:
        data[:] = dslice
    vv[:] = data[:]
f.close()