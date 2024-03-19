import numpy as np
import netCDF4 as nc
import os
import argparse
import xarray as xr
import xesmf as xe

parser = argparse.ArgumentParser(description='Interpolate zonal mean temperature-salinity restart file from GFDL CM2.1')
parser.add_argument("input", help="TS file from GFDL")
parser.add_argument("output", help="TS file on ACCESS-ESM1.5 ocean grid")
args = parser.parse_args()

gfdl_gridfile = 'ocean_hgrid_gfdl.nc'
access_gridfile = 'ocean_hgrid.nc'

f = nc.Dataset(gfdl_gridfile,'r')
x_g = f.variables['x'][:]
y_g = f.variables['y'][:]
f.close()

f = nc.Dataset(access_gridfile,'r')
x = f.variables['x'][:]
y = f.variables['y'][:]
f.close()

gfdl_grid = xr.Dataset(coords={
    "lat": (("y","x"), y_g[1::2,1::2]),
    "lon": (("y","x"), x_g[1::2,1::2])
    })

access_grid = xr.Dataset(coords={
    "lat": (("y","x"), y[1::2,1::2]),
    "lon": (("y","x"), x[1::2,1::2])    
    })

remap = xe.Regridder(gfdl_grid, access_grid, method='bilinear', periodic=True)

f = nc.Dataset(args.input,'r')
temp_g = f.variables['temp'][:]
salt_g = f.variables['salt'][:]
f.close()

junk, nz, ny_g, nx_g = temp_g.shape

ny, nx = access_grid.lat.shape

temp = np.zeros((1,nz,ny,nx), 'f8')
salt = np.zeros((1,nz,ny,nx), 'f8')

for k in range(nz):
    temp_in = temp_g.data[0,k,:,:]
    salt_in = salt_g.data[0,k,:,:]
    temp_out = remap(temp_in)
    salt_out = remap(salt_in)
    temp[0,k,:,:] = temp_out
    salt[0,k,:,:] = salt_out

fo = nc.Dataset(args.output,'w')
fo.history = f'interp_temp_salt_from_gfdl.py {args.input} {args.output} \n '
fo.title = args.output

fo.createDimension('xaxis_1',nx)
fo.createDimension('yaxis_1',ny)
fo.createDimension('zaxis_1',nz)
fo.createDimension('Time',0)

x_o = fo.createVariable('xaxis_1','f4',('xaxis_1'))
x_o.cartesian_axis = 'X'
x_o[:] = np.arange(nx)

y_o = fo.createVariable('yaxis_1','f4',('yaxis_1'))
y_o.cartesian_axis = 'Y'
y_o[:] = np.arange(ny)

z_o = fo.createVariable('zaxis_1','f4',('zaxis_1'))
z_o.cartesian_axis = 'Z'
z_o[:] = np.arange(nz)

t_o = fo.createVariable('Time','f8',('Time'))
t_o.cartesian_axis = 'T'
t_o[:] = 1

t_o = fo.createVariable('temp','f8',('Time','zaxis_1','yaxis_1','xaxis_1'))
t_o[:] = temp[:]

s_o = fo.createVariable('salt','f8',('Time','zaxis_1','yaxis_1','xaxis_1'))
s_o[:] = salt[:]

fo.close()