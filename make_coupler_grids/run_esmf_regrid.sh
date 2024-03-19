#!/bin/bash

ESMF_RegridWeightGen -s um1t_scrip.nc -d cice_scrip.nc -w rmp_um1t_to_cice_CONSERV_DESTAREA.nc -m conserve --norm_type dstarea 

ESMF_RegridWeightGen -s cice_scrip.nc -d um1t_scrip.nc -w rmp_cice_to_um1t_CONSERV_FRACNNEI.nc -m conserve --norm_type fracarea 
ESMF_RegridWeightGen -s cice_scrip.nc -d um1u_scrip.nc -w rmp_cice_to_um1u_CONSERV_FRACNNEI.nc -m conserve --norm_type fracarea 
ESMF_RegridWeightGen -s cice_scrip.nc -d um1v_scrip.nc -w rmp_cice_to_um1v_CONSERV_FRACNNEI.nc -m conserve --norm_type fracarea 
ESMF_RegridWeightGen -s um1t_scrip.nc -d cice_scrip.nc -w rmp_um1t_to_cice_CONSERV_FRACNNEI.nc -m conserve --norm_type fracarea 
ESMF_RegridWeightGen -s um1u_scrip.nc -d cice_scrip.nc -w rmp_um1u_to_cice_CONSERV_FRACNNEI.nc -m conserve --norm_type fracarea 
ESMF_RegridWeightGen -s um1v_scrip.nc -d cice_scrip.nc -w rmp_um1v_to_cice_CONSERV_FRACNNEI.nc -m conserve --norm_type fracarea 
