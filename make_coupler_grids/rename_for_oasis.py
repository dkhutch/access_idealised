#!/usr/bin/env python

import os
import numpy as np
import netCDF4 as nc
import argparse

parser = argparse.ArgumentParser(description='Update variable names from ESMF_RegridWeightGen for oasis conventions')
parser.add_argument('infile', help='input netcdf file to be updated')
args = parser.parse_args()

cmdstring = ('ncrename -O -d n_a,src_grid_size -d n_b,dst_grid_size -d n_s,'
             'num_links -d nv_a,src_grid_corners -d nv_b,dst_grid_corner'
             's -v yc_a,src_grid_center_lat -v yc_b,dst_grid_center_lat '
             '-v xc_a,src_grid_center_lon -v xc_b,dst_grid_center_lon -v'
             ' yv_a,src_grid_corner_lat -v xv_a,src_grid_corner_lon -v y'
             'v_b,dst_grid_corner_lat -v xv_b,dst_grid_corner_lon -v mas'
             'k_a,src_grid_imask -v mask_b,dst_grid_imask -v area_a,src_'
             'grid_area -v area_b,dst_grid_area -v frac_a,src_grid_frac '
             '-v frac_b,dst_grid_frac -v col,src_address -v row,dst_addr'
             'ess {}')
cmd = cmdstring.format(args.infile)
print(cmd)
os.system(cmd)

f = nc.Dataset(args.infile, 'r+')
hist = f.history[:]
hist = 'rename_for_oasis.py on {} \n'.format(args.infile) + hist
f.setncattr('history', hist)

remap_old = f.variables['S'][:]
remap_new = f.createVariable('remap_matrix', 'f8', ('num_links', 'num_wgts'))
remap_new[:,0] = remap_old[:]

f.close()
