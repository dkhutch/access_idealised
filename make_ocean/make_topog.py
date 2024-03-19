import numpy as np
import netCDF4 as nc
import xarray as xr
import xesmf as xe
import argparse

parser = argparse.ArgumentParser(description='Generate bathymetry on ocean grid using 1x1 deg input file')
parser.add_argument("input", help="1x1 deg input topography file")
parser.add_argument("output", help="bathymetry on ACCESS-ESM1.5 ocean grid")
args = parser.parse_args()

ocean_gridfile = 'ocean_hgrid.nc'

f = nc.Dataset(ocean_gridfile,'r')
oc_x = f.variables['x'][:]
oc_y = f.variables['y'][:]
f.close()

f = nc.Dataset(args.input,'r')
topo = f.variables['topo'][:]
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]
f.close()

topo = topo.astype('f8')
topo = topo.data
topo = -1. * topo

topo_grid = xe.util.grid_global(1., 1., lon1=360.)

oc_grid = xr.Dataset(coords={
    "lon": (("y","x"), oc_x[1::2,1::2]),
    "lat": (("y","x"), oc_y[1::2,1::2]),
    "lon_b": (("y_b","x_b"), oc_x[0::2,0::2]),
    "lat_b": (("y_b","x_b"), oc_y[0::2,0::2])  
    })

regridder = xe.Regridder(topo_grid, oc_grid, method="conservative_normed")

topo_oc = regridder(topo)
topo_oc[topo_oc < 0.] = 0.
min_depth = 40
ind = np.logical_and(topo_oc > 0, topo_oc < min_depth)
topo_oc[ind] = min_depth

ds_o = xr.Dataset(coords={
    "lon": (("y","x"), oc_x[1::2,1::2]),
    "lat": (("y","x"), oc_y[1::2,1::2])    
    },
    data_vars={
    "depth": (("y","x"), topo_oc)
    })
ds_o.attrs['history'] = 'make_topog.py on %s' % args.input
ds_o.to_netcdf(args.output)