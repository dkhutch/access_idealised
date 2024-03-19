import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, splrep, BSpline

infile = 'orog.esm1.5.nc'
maskfile = 'lsm_new.nc'
topofile = 'atmos_topog_idealised.nc'
outfile = 'xx_yy_scaled_idealised.nc'

ds_in = xr.open_dataset(infile)
grad_xx = ds_in['field152'].data
grad_yy = ds_in['field154'].data
silho = ds_in['field174'].data
peak_tr = ds_in['field175'].data
ht = ds_in['ht'].data

grad_xx = np.ma.masked_array(grad_xx, mask=np.isnan(grad_xx))
grad_yy = np.ma.masked_array(grad_yy, mask=np.isnan(grad_yy))
silho = np.ma.masked_array(silho, mask=np.isnan(silho))
peak_tr = np.ma.masked_array(peak_tr, mask=np.isnan(peak_tr))

ds_ideal = xr.open_dataset(topofile)
topo_ideal = ds_ideal['topog'].data
lat = ds_ideal['latitude'].data
lon = ds_ideal['longitude'].data

ds_mask = xr.open_dataset(maskfile)
lsm = ds_mask['lsm'].data
lsm = ~lsm.astype('bool')

# This applies if lon and lat are 1-Dimensional, otherwise comment out!
lon, lat = np.meshgrid(lon, lat)

wts = np.cos(np.deg2rad(lat))

ht_min = 0.
ht_max = 5500.

hist, bin_e = np.histogram(ht, bins=55, range=(ht_min, ht_max))
nbin = hist.shape[0]

xx_bin = np.zeros(nbin)
yy_bin = np.zeros(nbin)
silho_bin = np.zeros(nbin)
peak_tr_bin = np.zeros(nbin)

for i in range(nbin):
    ind = np.logical_and(ht >= bin_e[i], ht < bin_e[i+1])
    xx_bin[i] = np.ma.average(grad_xx[ind], weights=wts[ind])
    yy_bin[i] = np.ma.average(grad_yy[ind], weights=wts[ind])
    silho_bin[i] = np.ma.average(silho[ind], weights=wts[ind])
    peak_tr_bin[i] = np.ma.average(peak_tr[ind], weights=wts[ind])

ht_bin = (bin_e[1:] + bin_e[:-1]) * 0.5

sparm = 1e-4

tck_xx = splrep(ht_bin, xx_bin, s=sparm)
tck_yy = splrep(ht_bin, yy_bin, s=sparm)
tck_si = splrep(ht_bin, silho_bin, s=sparm)
tck_pe = splrep(ht_bin, peak_tr_bin, s=sparm)

ideal_xx = BSpline(*tck_xx)(topo_ideal)
ideal_yy = BSpline(*tck_yy)(topo_ideal)
ideal_si = BSpline(*tck_si)(topo_ideal)
ideal_pe = BSpline(*tck_pe)(topo_ideal)

ideal_xx[ideal_xx < 0.] = 0.
ideal_yy[ideal_yy < 0.] = 0.
ideal_si[ideal_si < 0.] = 0.
ideal_pe[ideal_pe < 0.] = 0.

ideal_xx = np.ma.masked_array(ideal_xx, mask=lsm)
ideal_yy = np.ma.masked_array(ideal_yy, mask=lsm)
ideal_si = np.ma.masked_array(ideal_si, mask=lsm)
ideal_pe = np.ma.masked_array(ideal_pe, mask=lsm)

ds_o = xr.Dataset(coords={
    "lat": ("lat", lat[:,0]),
    "lon": ("lon", lon[0,:])
    },
    data_vars={
    "grad_xx" : (("lat","lon"), ideal_xx),
    "grad_yy" : (("lat","lon"), ideal_yy),
    "silhouette" : (("lat","lon"), ideal_si),
    "peak_trough" : (("lat","lon"), ideal_pe)
    })
ds_o.lat.attrs['units'] = 'degrees_north'
ds_o.lon.attrs['units'] = 'degrees_east'
ds_o.attrs['history'] = f'scale_xx_yy_grads.py on {infile} and {topofile} \n'
ds_o.to_netcdf(outfile)


test_plot = False
if test_plot:
    ht_new = np.arange(0., 5501., 10.)

    spl_xx = BSpline(*tck_xx)(ht_new)
    spl_yy = BSpline(*tck_yy)(ht_new)
    spl_si = BSpline(*tck_si)(ht_new)
    spl_pe = BSpline(*tck_pe)(ht_new)

    clist = ['red','green','black','cyan']

    plt.plot(ht_new, spl_xx, color=clist[0])
    plt.plot(ht_bin, xx_bin, 'o', ms=5, color=clist[0])

    plt.plot(ht_new, spl_yy, color=clist[1])
    plt.plot(ht_bin, yy_bin, 'o', ms=5, color=clist[1])

    plt.plot(ht_new, spl_si, color=clist[2])
    plt.plot(ht_bin, silho_bin, 'o', ms=5, color=clist[2])

    # plt.plot(ht_new, spl_pe, color=clist[3])
    # plt.plot(ht_bin, peak_tr_bin, 'o', ms=5, color=clist[3])

    plt.legend(['xx spline','xx bins','yy spline','yy bins',
                 'sil spline','sil bins']) #,'pt spline','pt bins'])

    plt.xlabel('Altitude (m)')
    plt.ylabel('Gradient magnitude')
    plt.savefig('spline_fit_gradients.pdf')