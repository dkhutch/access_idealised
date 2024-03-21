#!/bin/bash

EXP='solo150'
INPUTDIR='/g/data/w40/dxd565/gfdl-model/experiments/'${EXP}'/input/'
CODEDIR='/g/data/w40/dkh157/access_idealised/'
ANCILDIR='/g/data/w40/dkh157/access_idealised/ancils/'${EXP}'/'

export PYTHONPATH=/g/data/access/projects/access/apps/pythonlib/umfile_utils/:$PYTHONPATH

#--------------------
# Step 1: Make ocean
#--------------------
cd $CODEDIR/make_ocean
python make_topog.py $INPUTDIR/topo.nc ./topog.nc
python fix_kmt.py
python interp_temp_salt_from_gfdl.py $INPUTDIR/ocean_temp_salt.res.nc ocean_temp_salt.res.nc
python bgc_restart.py
python bgc_sediment.py
python ssw_fix.py
python basin_mask_fix.py

#---------------------------
# Step 2: Make coupler grids
#---------------------------
cd $CODEDIR/make_coupler_grids
python regrid_um_mask.py
ncl make_scrip_files.ncl
./run_esmf_regrid.sh
./run_rename_oasis.sh

#--------------------------------
# Step 3: Make atmosphere restart
#--------------------------------
cd $CODEDIR/make_atmos
python make_lsm.py
python atmos_topog_from_gfdl.py $INPUTDIR/fv_rst.res.nc atmos_topog_idealised.nc
python remake_stddev.py
python scale_xx_yy_grads.py
python interp_atmos_fields.py
python interp_fields_arg.py qrclim.slt.orig qrclim.slt
python interp_fields_arg.py qrclim.smow.orig qrclim.smow
python interp_fields_arg.py cable_vegfunc_N96.anc.orig cable_vegfunc_N96.anc
python interp_fields_arg.py qrparm.soil_igbp_vg.orig qrparm.soil_igbp_vg
python interp_fields_arg.py qrparm.mask.orig qrparm.mask

#-----------------------------------
# Step 4: Make coupler restart files
#-----------------------------------
cd $CODEDIR/coupler_inputs
python fix_a2i.py
python fix_i2a.py
python fix_o2i.py

#--------------------------------------
# Step 5: Make additional ocean forcing
#--------------------------------------
cd $CODEDIR/oc_common
python dust.py
python fice.py
python press.py
python xkw.py

#--------------------------------------
# Step 6: Collect files into $ANCILDIR
#--------------------------------------
mkdir -p $ANCILDIR
cd $CODEDIR
cp -p \
    ./make_ocean/topog.nc \
    ./make_ocean/kmt.nc \
    ./make_ocean/ocean_temp_salt.res.nc \
    ./make_ocean/ssw_atten_depth.nc \
    ./make_ocean/basin_mask.nc \
    ./make_ocean/csiro_bgc.res.nc \
    ./make_ocean/csiro_bgc_sediment.res.nc \
    ./make_coupler_grids/landfrac_um1t.nc \
    ./make_coupler_grids/masks.nc \
    ./make_coupler_grids/um1v_scrip.nc \
    ./make_coupler_grids/um1u_scrip.nc \
    ./make_coupler_grids/um1t_scrip.nc \
    ./make_coupler_grids/cice_scrip.nc \
    ./make_coupler_grids/rmp_cice_to_um1t_CONSERV_FRACNNEI.nc \
    ./make_coupler_grids/rmp_cice_to_um1v_CONSERV_FRACNNEI.nc \
    ./make_coupler_grids/rmp_cice_to_um1u_CONSERV_FRACNNEI.nc \
    ./make_coupler_grids/rmp_um1t_to_cice_CONSERV_DESTAREA.nc \
    ./make_coupler_grids/rmp_um1u_to_cice_CONSERV_FRACNNEI.nc \
    ./make_coupler_grids/rmp_um1t_to_cice_CONSERV_FRACNNEI.nc \
    ./make_coupler_grids/rmp_um1v_to_cice_CONSERV_FRACNNEI.nc \
    ./make_atmos/lsm_new.nc \
    ./make_atmos/atmos_topog_idealised.nc \
    ./make_atmos/stddev_idealised.nc \
    ./make_atmos/xx_yy_scaled_idealised.nc \
    ./make_atmos/land_field.txt \
    ./make_atmos/restart.idealised \
    ./make_atmos/qrclim.slt \
    ./make_atmos/qrclim.smow \
    ./make_atmos/cable_vegfunc_N96.anc \
    ./make_atmos/qrparm.soil_igbp_vg \
    ./make_atmos/qrparm.mask \
    ./coupler_inputs/a2i.nc \
    ./coupler_inputs/i2a.nc \
    ./coupler_inputs/o2i.nc \
    ./oc_common/dust.nc \
    ./oc_common/ocmip2_fice_monthly_om1p5_bc.nc \
    ./oc_common/ocmip2_press_monthly_om1p5_bc.nc \
    ./oc_common/ocmip2_xkw_monthly_om1p5_bc.nc \
    $ANCILDIR
