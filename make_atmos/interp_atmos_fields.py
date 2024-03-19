import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import iris, umfile
from um_fileheaders import *
from stashvar import atm_stashvar
import os
import xesmf as xe
import xarray as xr

restartfile = 'restart.subset'
new_restart = 'restart.idealised'
lsm_pi_file = 'lsm_esm1.5.nc'
lsm_new_file = 'lsm_new.nc'
maskvar = 'lsm'
netcdf_landfrac = '../make_coupler_grids/landfrac_um1t.nc'
mask_file = '../make_coupler_grids/masks.nc'
mask_file_orig = '../make_coupler_grids/masks_orig.nc'
orog_file = 'atmos_topog_idealised.nc'
stddev_file = 'stddev_idealised.nc'
gradient_file = 'xx_yy_scaled_idealised.nc'

f = nc.Dataset(lsm_pi_file, 'r')
lsm_pi = f.variables['lsm'][:]
f.close()
lsm_pi = lsm_pi.astype('i4')

f = nc.Dataset(lsm_new_file, 'r')
lsm_new = f.variables['lsm'][:]
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]
f.close()
lsm_new = lsm_new.astype('i4')

nlat, nlon = lsm_pi.shape

f = nc.Dataset(mask_file,'r')
at_mask_cpl = f.variables['um1t.msk'][:]
f.close()

f = nc.Dataset(mask_file_orig,'r')
at_mask_cpl_orig = f.variables['um1t.msk'][:]
f.close()

f = nc.Dataset(netcdf_landfrac,'r')
landfrac_new = f.variables['landfrac'][:]
f.close()

f = nc.Dataset(orog_file, 'r')
topo = f.variables['topog'][:]
f.close()

f = nc.Dataset(stddev_file,'r')
stddev = f.variables['stddev'][:]
f.close()

f = nc.Dataset(gradient_file,'r')
grad_xx = f.variables['grad_xx'][:]
grad_yy = f.variables['grad_yy'][:]
silhouette = f.variables['silhouette'][:]
peak_trough = f.variables['peak_trough'][:]
f.close()

new_ocean = np.logical_and(lsm_pi==1, lsm_new==0)
new_land = np.logical_and(lsm_pi==0, lsm_new==1)
land_pts = lsm_new==1

seaice = at_mask_cpl==0
no_seaice = at_mask_cpl==1

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

os.system('cp {} {}'.format(restartfile, new_restart))
os.system('um_replace_field.py -v 30 -n {} -V {} {}'.format(
    lsm_new_file, maskvar, new_restart))

fin = umfile.UMFile(restartfile, 'r')
f = umfile.UMFile(new_restart, 'r+')
nvars = f.fixhd[FH_LookupSize2]

lsm_code = 30
landfrac_code = 505
snow_codes = [23, 95, 416]
ice_reset = [49, 415]
ice_zeros = [31, 32, 413, 414, 416, 509]
ice_temp = 508
oc_curr = [28, 29, 269, 270]
orog_code = 33
stddev_code = 34
grad_xx_code = 35
grad_xy_code = 36
grad_yy_code = 37
sil_code = 17
peak_code = 18

ice_reset_val = 271.35

# orog_new = 10.
# orog_std = 10.
# orog_grad = 1.e-5
# orog_sil = 1.e-2
# orog_vars = [17, 18, 33, 34, 35, 36, 37]

for k in range(nvars):
    ilookup = fin.ilookup[k]
    lbegin = ilookup[LBEGIN] 
    if lbegin == -99:
        break
    a = fin.readfld(k)

    if ilookup[ITEM_CODE] in oc_curr:
        a[:] = 0.
    if a.shape == (nlat, nlon):
        if ilookup[ITEM_CODE] == lsm_code:
            continue
        elif ilookup[ITEM_CODE] == landfrac_code:
            a[:] = landfrac_new[:]
        elif ilookup[ITEM_CODE] in snow_codes:
            a[:] = 0.
        elif ilookup[ITEM_CODE] in ice_zeros:
            a[:] = 0.
        elif ilookup[ITEM_CODE] in ice_reset:
            a[seaice] = ice_reset_val
            a[no_seaice] = fin.missval_r
        elif ilookup[ITEM_CODE] == ice_temp:
            a[:] = ice_reset_val
        elif ilookup[ITEM_CODE] == orog_code:
            a[:] = topo
        elif ilookup[ITEM_CODE] == stddev_code:
            a[:] = stddev
        elif ilookup[ITEM_CODE] == grad_xx_code:
            a[:] = grad_xx
        elif ilookup[ITEM_CODE] == grad_yy_code:
            a[:] = grad_yy
        elif ilookup[ITEM_CODE] == grad_xy_code:
            a[:] = 0.
        elif ilookup[ITEM_CODE] == sil_code:
            a[:] = silhouette
        elif ilookup[ITEM_CODE] == peak_code:
            a[:] = peak_trough
        elif ilookup[LBPACK]==120:
            a_new = regridder(a)
            a[new_land] = a_new[new_land]
            a[new_ocean] = fin.missval_r
    f.writefld(a[:], k)
f.close()

