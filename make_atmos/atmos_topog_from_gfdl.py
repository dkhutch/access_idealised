import numpy as np
import netCDF4 as nc 
import xarray as xr 
import xesmf as xe
import argparse 

parser = argparse.ArgumentParser(description='Generate atmosphere topography on UM grid from GFDL model')
parser.add_argument("input", help="3x3.75 deg restart file from GFDL model")
parser.add_argument("output", help="topography on N96 grid for ACCESS-ESM1.5")
args = parser.parse_args()

newlsmfile = 'lsm_new.nc'

# run atmos_topog_from_gfdl.py /g/data/w40/dxd565/gfdl-model/experiments/solo150/input/fv_rst.res.nc atmos_topog_idealised.nc

f = nc.Dataset(args.input, 'r')
geopot = f.variables['Surface_geopotential'][:]
lat_in = f.variables['lat'][:]
lon_in = f.variables['lon'][:]
f.close()

ht_old = geopot / 9.8
ht_old = ht_old.data
ht_old = np.squeeze(ht_old)

f = nc.Dataset(newlsmfile,'r')
lat_um = f.variables['latitude'][:]
lon_um = f.variables['longitude'][:]
lsm = f.variables['lsm'][:]
f.close()

lon2d_in, lat2d_in = np.meshgrid(lon_in, lat_in)
lon2d_um, lat2d_um = np.meshgrid(lon_um, lat_um)

old_grid = xr.Dataset(coords={
    "lat": (("y","x"), lat2d_in),
    "lon": (("y","x"), lon2d_in)
})

new_grid = xr.Dataset(coords={
    "lat": (("y","x"), lat2d_um),
    "lon": (("y","x"), lon2d_um)
})

remap = xe.Regridder(old_grid, new_grid, method='bilinear', periodic=True)

topog = remap(ht_old)
# clean up outside of LSM
topog[lsm == 0] = 0.
# Smooth Antarctica polar cells:
topo_val = np.mean(topog[5,:])
topog[:4,:] = topo_val

ds = xr.Dataset(coords={
    "latitude": (("latitude"), lat_um),
    "longitude": (("longitude"), lon_um)
    },
    data_vars={
    "topog": (("latitude","longitude"), topog)
    })
encoding = {}
for key in ds.coords:
    encoding[key] = {"_FillValue": None}
for key in ds.data_vars:
    encoding[key] = {"_FillValue": -1.0e20}
ds.attrs['history'] = 'atmos_topog_from_gfdl.py'
ds.to_netcdf(args.output, encoding=encoding)
