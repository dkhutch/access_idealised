import numpy as np
import netCDF4 as nc
import xarray as xr
import xesmf as xe
import matplotlib.pyplot as plt
import os

oasis_atmosfile = 'grids.nc'
masks_infile = 'masks_orig.nc'
oceangridfile = '../make_ocean/ocean_hgrid.nc'
topofile = '../make_ocean/topog.nc'
masks_new = 'masks.nc'

landfracfile = 'landfrac_um1t.nc'

f = nc.Dataset(topofile,'r')
depth = f.variables['depth'][:]
f.close()

omask = depth > 0 
omask = omask.astype('f8')
momt_mask = depth==0
momt_mask = momt_mask.astype('i4')


f = nc.Dataset(oceangridfile,'r')
oc_x = f.variables['x'][:]
oc_y = f.variables['y'][:]
f.close()

# Create the ocean grid with mask variable:
oc_grid = xr.Dataset(coords={
    "lon": (("y","x"), oc_x[1::2,1::2]),
    "lat": (("y","x"), oc_y[1::2,1::2]),
    "lon_b": (("y_b","x_b"), oc_x[0::2,0::2]),
    "lat_b": (("y_b","x_b"), oc_y[0::2,0::2])    
    },
    data_vars={
    "omask": (("y","x"), omask)
    })

os.system('cp {} {}'.format(masks_infile, masks_new))

f_masks = nc.Dataset(masks_new, 'r+')
hist = 'regrid_um_mask.py on {} \n'.format(topofile)
f_masks.setncattr('history', hist)

vnames = ['um1t','um1v','um1u']
for vname in vnames:

    f = nc.Dataset(oasis_atmosfile,'r')
    at_lon_c = f.variables[f'{vname}.clo'][:]
    at_lat_c = f.variables[f'{vname}.cla'][:]
    at_lon = f.variables[f'{vname}.lon'][:]
    at_lat = f.variables[f'{vname}.lat'][:]
    f.close()

    nlat_a, nlon_a = at_lon.shape

    at_lon_b = np.zeros((nlat_a+1, nlon_a+1), 'f8')
    at_lat_b = np.zeros((nlat_a+1, nlon_a+1), 'f8')

    at_lon_b[0:-1,0:-1] = at_lon_c[0,:,:]
    at_lon_b[1:,1:] = at_lon_c[2,:,:]
    at_lon_b[-1,0] = at_lon_c[3,-1,0]
    at_lon_b[0,-1] = at_lon_c[1,0,-1]

    at_lat_b[0:-1,0:-1] = at_lat_c[0,:,:]
    at_lat_b[1:,1:] = at_lat_c[2,:,:]
    at_lat_b[-1,0] = at_lat_c[3,-1,0]
    at_lat_b[0,-1] = at_lat_c[1,0,-1]

    # Create the atmos grid for interpolation:
    at_grid = xr.Dataset(coords={
        "lon": (("y","x"), at_lon),
        "lat": (("y","x"), at_lat),
        "lon_b": (("y_b","x_b"), at_lon_b),
        "lat_b": (("y_b","x_b"), at_lat_b)    
        })

    # Generate the regridder function
    regridder = xe.Regridder(oc_grid, at_grid, method="conservative_normed")

    oc_in = oc_grid["omask"]
    at_mask = regridder(oc_in)
    at_mask = 1 - at_mask
    at_mask = at_mask.rename("landfrac")
    at_mask = at_mask.rename({'x': 'longitude','y': 'latitude'})
    at_mask = at_mask.assign_coords({
        "latitude": (("latitude"), at_mask.lat[:,0].data),
        "longitude": (("longitude"), at_mask.lon[0,:].data)
        })
    at_mask.attrs['history'] = 'regrid_um_mask.py \n '

    thresh = 0.01 # Minimum threshold for land fraction (a somewhat arbitrary number!)
    at_mask = at_mask.where(at_mask >= thresh, 0.)
    # at_mask = at_mask.where(at_mask < 0.99, 1.)

    cutoff = 0.99

    msk = at_mask.where(at_mask < cutoff, 1.)
    msk = msk.where(at_mask > cutoff, 0.)
    msk = msk.astype('i4')
    msk = msk.rename("msk")

    varstring = '{}.msk'.format(vname)
    outvar = f_masks.variables[varstring]
    outvar[:] = msk[:]

    if vname=='um1t':
        at_mask.to_netcdf(landfracfile)

outvar = f_masks.variables['cice.msk']
outvar[:] = momt_mask[:]

f_masks.close()

