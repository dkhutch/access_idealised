import numpy as np
import netCDF4 as nc
import os

mask_file = '../make_coupler_grids/masks.nc'
mask_file_orig = '../make_coupler_grids/masks_orig.nc'
a2ifile = 'a2i_orig.nc'
a2inew = 'a2i.nc'

os.system('cp {} {}'.format(a2ifile, a2inew))

f = nc.Dataset(mask_file,'r')
at_mask_cpl = f.variables['um1t.msk'][:]
f.close()

f = nc.Dataset(mask_file_orig,'r')
at_mask_cpl_orig = f.variables['um1t.msk'][:]
f.close()

nlat, nlon = at_mask_cpl.shape

new_ocean_um1t = np.logical_and(at_mask_cpl==0, at_mask_cpl_orig==1)
new_land_um1t = np.logical_and(at_mask_cpl==1, at_mask_cpl_orig==0)

varlist = ['bmlt01',
'bmlt02',
'bmlt03',
'bmlt04',
'bmlt05',
'evap2d',
'heatflux',
'lhflx',
'lwflx',
'pen_sol',
'press',
'runoff',
'shflx',
'swflx',
'taux',
'tauy',
'tmlt01',
'tmlt02',
'tmlt03',
'tmlt04',
'tmlt05',
'train',
'tsnow',
'wme',
'co2_a',
'wnd_a']

f = nc.Dataset(a2inew, 'r+')
if 'history' in f.ncattrs(): 
    hist = f.history[:]
else:
    hist = ''
hist = 'fix_a2i_mio.py \n ' + hist
f.setncattr('history',hist)
f.setncattr('title',a2inew)

src_y = 76
src_x = 145

for vname in varlist:
    vv = f.variables[vname]
    data = vv[:]
    if vname in ['wnd_a','tauy','co2_a']:
        continue
    elif vname[:4] == 'bmlt' or vname=='tsnow':
        data[:] = 0.
    else:
        data[new_land_um1t] = 0.
        if vname == 'press':
            data[new_ocean_um1t] = data[src_y, src_x]
    vv[:] = data[:]
f.close()