import numpy as np
import netCDF4 as nc
import os
from scipy.interpolate import griddata

infile = 'ocean_temp_salt.res.orig.nc'
oceangrid = 'ocean_hgrid.nc'
outfile = 'ocean_temp_salt.res.nc'
topogfile = 'topog.nc'

f = nc.Dataset(oceangrid,'r')
x = f.variables['x'][:]
y = f.variables['y'][:]
f.close()

xc = x[1::2,1::2]
yc = y[1::2,1::2]

f = nc.Dataset(topogfile,'r')
depth = f.variables['depth'][:]
f.close()

land_mask = depth == 0

cmd = f'cp {infile} {outfile}'
os.system(cmd)

f = nc.Dataset(outfile,'r+')
hist = f.history[:]
hist = 'fill_restart_gaps.py \n' + hist
f.setncattr('history', hist)

temp = f.variables['temp']
salt = f.variables['salt']

tdata = temp[:]
sdata = salt[:]

salt_thresh = 20. # minimum value of salt restart

junk, nz, ny, nx = temp.shape

for k in range(nz):
    tslice = temp[0,k,:,:]
    sslice = salt[0,k,:,:]
    mask = sslice!=0

    tmask = tslice[mask].ravel()
    smask = sslice[mask].ravel()
    xmask = xc[mask].ravel()
    ymask = yc[mask].ravel()

    t_itp = griddata((xmask, ymask), tmask, (xc, yc), method='nearest')
    s_itp = griddata((xmask, ymask), smask, (xc, yc), method='nearest')

    s_itp[s_itp < salt_thresh] = salt_thresh

    t_itp[land_mask] = 0.
    s_itp[land_mask] = 0.

    tdata[0,k,:,:] = t_itp
    sdata[0,k,:,:] = s_itp

temp[:] = tdata[:]
salt[:] = sdata[:]

f.close()