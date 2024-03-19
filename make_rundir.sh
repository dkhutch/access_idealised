#!/bin/bash

# PROJECT='y99' # Here you can specify a project or else it uses your default project
EXP='solo150'
ANCILDIR='/g/data/w40/dkh157/access_idealised/ancils/'${EXP}'/'
INPUTDIR="/g/data/$PROJECT/$USER/access_esm/input/${EXP}/"
RESTARTDIR="/g/data/$PROJECT/$USER/access_esm/restart/${EXP}/"
RUNDIR='/home/157/dkh157/ACCESS/'${EXP}'/'

BASE_INPUT='/g/data/w40/dkh157/access_esm_pi/input/'
BASE_RESTART='/g/data/w40/dkh157/access_esm_pi/restart/'
BASE_RUNDIR='/home/157/dkh157/ACCESS/PI/'

mkdir -p $INPUTDIR $RESTARTDIR $RUNDIR

cd $BASE_INPUT
cp -rp atmosphere coupler ice ocean start_dump $INPUTDIR
cd $BASE_RESTART 
cp -rp atmosphere coupler ice ocean $RESTARTDIR

cd $ANCILDIR
# RESTART FILES
cp -p restart.idealised $RESTARTDIR/atmosphere/restart_dump.astart
cp -p a2i.nc o2i.nc i2a.nc $RESTARTDIR/coupler
cp -p ocean_temp_salt.res.nc csiro_bgc.res.nc csiro_bgc_sediment.res.nc $RESTARTDIR/ocean

# INPUT FILES
cp -p qrclim.slt qrclim.smow qrparm.mask cable_vegfunc_N96.anc qrparm.soil_igbp_vg $INPUTDIR/atmosphere
cp -p rmp_cice_to_um1t_CONSERV_FRACNNEI.nc \
    rmp_cice_to_um1u_CONSERV_FRACNNEI.nc \
    rmp_cice_to_um1v_CONSERV_FRACNNEI.nc \
    rmp_um1t_to_cice_CONSERV_DESTAREA.nc \
    rmp_um1t_to_cice_CONSERV_FRACNNEI.nc \
    rmp_um1u_to_cice_CONSERV_FRACNNEI.nc \
    rmp_um1v_to_cice_CONSERV_FRACNNEI.nc \
    masks.nc \
    $INPUTDIR/coupler
cp -p kmt.nc $INPUTDIR/ice
cp -p basin_mask.nc ssw_atten_depth.nc topog.nc $INPUTDIR/ocean/common
cp -p dust.nc ocmip2_fice_monthly_om1p5_bc.nc ocmip2_press_monthly_om1p5_bc.nc ocmip2_xkw_monthly_om1p5_bc.nc $INPUTDIR/ocean/pre-industrial
cp -p restart.idealised $INPUTDIR/start_dump/pre-industrial.astart

# RUNDIR
cd $BASE_RUNDIR
cp -rp atmosphere coupler ice ocean config.yaml $RUNDIR
cd $RUNDIR
landpts=`cat $ANCILDIR/land_field.txt`
# Now we edit the number of land points in namelist files
sed -i "s/LAND_FIELD=10865/LAND_FIELD=${landpts}/g" $RUNDIR/atmosphere/namelists
sed -i "s/LAND_FIELD=10865/LAND_FIELD=${landpts}/g" $RUNDIR/atmosphere/SIZES

# Now change paths in config.yaml (CHECK RESULT BEFORE RUNNING)
sed -i "s@jobname: esm-pi@jobname: ${EXP}@g" $RUNDIR/config.yaml
sed -i "s@/g/data/access/payu/access-esm/input/pre-industrial@/g/data/$PROJECT/$USER/access_esm/input/${EXP}@g" $RUNDIR/config.yaml
sed -i "s@/g/data/access/payu/access-esm/restart/pre-industrial@/g/data/$PROJECT/$USER/access_esm/restart/${EXP}@g" $RUNDIR/config.yaml
