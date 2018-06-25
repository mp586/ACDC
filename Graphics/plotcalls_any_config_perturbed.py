from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7 import *
import stats as st

GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
import cell_area as ca

ctl_model = input('Enter model name as string ')
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'


control_dir= control_model + '/' + input('Enter control directory name as string ')
print control_dir
ctl_runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
ctl_runmax=input('Enter runmax number for comparison ')
ctl_timeseries_max = input('Enter end of ctl timeseries month ')

model = input('Enter model ')
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir = ''

testdir= input('Enter data directory name as string ')
runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number ')


outdir = output_dir + '/' + testdir
testdir = model_data + '/' + testdir

land = input('Which landmask? ')
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

area_array = ca.cell_area(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/')
area_array = xr.DataArray(area_array)

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)



total_sfc_area = np.sum(area_array)
print ('total sfc area (*10^14) = '+str(np.sum(area_array/(10**14)))) # -- test: correct, equals sfc area of earth (5.1*10**14 m^2)
land_sfc_area = np.sum(area_array.where(landmask==1.))
print ('land sfc area (*10^14) = '+str(land_sfc_area/(10**14)))
ocean_sfc_area = np.sum(area_array.where(landmask!=1.))
print ('ocean sfc area (*10^14) = '+str(ocean_sfc_area/(10**14)))

# for plotting a spin up run ('control') timeseries followed by the timeseries from the perturbed experiment
globavg_var_timeseries_total_and_land_perturbed(testdir,model,area_array,'t_surf',1,runmax,1.,landmask,control_dir,ctl_model,1,ctl_timeseries_max)
# globavg_var_timeseries_total_and_land_perturbed(testdir,model,area_array,'bucket_depth',1,runmax,1.,landmask,control_dir,ctl_model,1,ctl_timeseries_max,select='land')


[msf,msf_avg,msf_seasonal_avg,msf_month_avg] = mass_streamfunction(testdir,model,runmin,runmax) # 
plot_streamfunction_seasonal(msf_seasonal_avg)

# globavg_var_timeseries(testdir,'co2',1,runmax)


[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','W/m^2',factor = 1/28.) # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','kg/m2s', factor=86400)
#[bucket_depth,bucket_depth_avg,bucket_depth_seasonal_avg,bucket_depth_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'bucket_depth','m')
[rh,rh_avg,rh_seasonal_avg,rh_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'rh','%',level=39)
#[CIWV,CIWV_avg,CIWV_seasonal_avg,CIWV_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level='all')

[net_sw,net_sw_avg,net_sw_seasonal_avg,net_sw_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_sw','W/m^2',factor = 1.) # 
[net_lw,net_lw_avg,net_lw_seasonal_avg,net_lw_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lw','W/m^2',factor = 1.) # 
[net_t,net_t_avg,net_t_seasonal_avg,net_t_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_t','W/m^2',factor = 1.) # 

#[convection_rain,convection_rain_avg,convection_rain_seasonal_avg,convection_rain_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'convection_rain','kg/m2s')
#[condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,condensation_rain_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'condensation_rain','kg/m2s')



[tsurf_ctl,tsurf_avg_ctl,tsurf_seasonal_avg_ctl,tsurf_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[net_lhe_ctl,net_lhe_avg_ctl,net_lhe_seasonal_avg_ctl,net_lhe_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','W/m^2',factor = 1/28.) # latent heat flux at surface (UP)
[precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','kg/m2s', factor=86400)
#[bucket_depth_ctl,bucket_depth_avg_ctl,bucket_depth_seasonal_avg_ctl,bucket_depth_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'bucket_depth','m')
[rh_ctl,rh_avg_ctl,rh_seasonal_avg_ctl,rh_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%',level=39)
#[CIWV_ctl,CIWV_avg_ctl,CIWV_seasonal_avg_ctl,CIWV_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level='all')
#[flux_oceanq_ctl,flux_oceanq_avg_ctl,flux_oceanq_seasonal_avg_ctl,flux_oceanq_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_oceanq','W/m^2')
[net_sw_ctl,net_sw_avg_ctl,net_sw_seasonal_avg_ctl,net_sw_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_sw','W/m^2',factor = 1.) # 
[net_lw_ctl,net_lw_avg_ctl,net_lw_seasonal_avg_ctl,net_lw_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lw','W/m^2',factor = 1.) # 
[net_t_ctl,net_t_avg_ctl,net_t_seasonal_avg_ctl,net_t_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_t','W/m^2',factor = 1.) # 


# ctl_runmin = 421
# ctl_runmax = 481
# runmin = 1
# runmax = 60
# [precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','kg/m2s', factor=86400)
# [precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','kg/m2s', factor=86400)


# for i in range(ctl_runmin-ctl_runmin+1,ctl_runmax-ctl_runmin):
#     month_plot = any_configuration_plot(outdir,ctl_runmin,ctl_runmax,-90.,90.,precipitation_ctl[i,:,:],
#                                         area_array,'mm/day',
#                                         'P_avg_ctl','fromwhite',landmaskxr,landlats,landlons,nmb_contours=10,
#                                         minval = 0., maxval = 8.,month_annotate=i)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/P_month'+str(1000+i)+'.png',bbox_inches='tight')
# for i in range(runmin-runmin+1,runmax-runmin):
#     month_plot = any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation[i,:,:],area_array,'mm/day',
#                                         'P_avg','fromwhite',landmaskxr,landlats,landlons,nmb_contours=10,
#                                         minval = 0., maxval = 8.,month_annotate=i)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/P_month'+str(2000+i)+'.png',bbox_inches='tight')

# os.system('convert -delay 50 /scratch/mp586/Code/Graphics/'+outdir+'/P_month*.png /scratch/mp586/Code/Graphics/'+outdir+'/P_month_animation.gif')




rh_P_E_change(outdir,runmin,runmax,rh_avg,rh_avg_ctl,precipitation_avg,precipitation_avg_ctl,net_lhe_avg,net_lhe_avg_ctl,tsurf_avg,tsurf_avg_ctl,landmask,sfc='all')
rh_P_E_change(outdir,runmin,runmax,rh_avg,rh_avg_ctl,precipitation_avg,precipitation_avg_ctl,net_lhe_avg,net_lhe_avg_ctl,tsurf_avg,tsurf_avg_ctl,landmask,sfc='land')
rh_P_E_change(outdir,runmin,runmax,rh_avg,rh_avg_ctl,precipitation_avg,precipitation_avg_ctl,net_lhe_avg,net_lhe_avg_ctl,tsurf_avg,tsurf_avg_ctl,landmask,sfc='ocean')


##### DECOMPOSITION OF PRECIPITATION CHANGE ######################


[sphum,sphum_avg,sphum_seasonal_avg,sphum_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level=39)
[sphum_ctl,sphum_avg_ctl,sphum_seasonal_avg_ctl,sphum_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level=39)

M_ctl = precipitation_ctl / sphum_ctl
M_new = precipitation / sphum
Delta_M = M_new - M_ctl
Delta_P_circ = sphum_ctl * Delta_M
Delta_P_q = M_ctl * (sphum - sphum_ctl)
Delta_P_cross = (sphum - sphum_ctl) * (M_new - M_ctl)
Delta_P_sum = Delta_P_cross + Delta_P_q + Delta_P_circ
Delta_P_sum_avg = Delta_P_sum.mean('time')
Diff_DeltaPs = (precipitation_avg - precipitation_avg_ctl) - Delta_P_sum_avg
print('Diff_DeltaPs.max()'+str(Diff_DeltaPs.max()))
print('Diff_DeltaPs.min()'+str(Diff_DeltaPs.min()))
print('Diff_DeltaPs.mean()'+str(Diff_DeltaPs.mean()))

any_configuration_plot(outdir,runmin,runmax,-90.,90.,Delta_P_circ.mean('time'),area_array,'mm/day','Delta_P_circ','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,Delta_P_q.mean('time'),area_array,'mm/day','Delta_P_q','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,Delta_P_cross.mean('time'),area_array,'mm/day','Delta_P_cross','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,Diff_DeltaPs,area_array,'mm/day','Diff_Delta_Ps','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,Delta_P_circ.mean('time'),area_array,'mm/day','Delta_P_circ_narrowcbar','rainnorm',landmaskxr,landlats,landlons,minval=-1.,maxval=1.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,Delta_P_q.mean('time'),area_array,'mm/day','Delta_P_q_narrowcbar','rainnorm',landmaskxr,landlats,landlons,minval=-1.,maxval=1.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,Delta_P_cross.mean('time'),area_array,'mm/day','Delta_P_cross_narrowcbar','rainnorm',landmaskxr,landlats,landlons,minval=-1.,maxval=1.)

globavg_deltaT = area_weighted_avg((tsurf_avg - tsurf_avg_ctl),area_array,landmask,'all_sfcs')

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(Delta_M/globavg_deltaT).mean('time'),area_array,'mm/day/K','Delta_M_normalized','rainnorm',landmaskxr,landlats,landlons,minval=-40.,maxval=40) # change in circulation normalized by global mean T change # following Chadwick et al., 2013 (Fig. 8)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,sphum_avg - sphum_avg_ctl,area_array,'kg/kg','q_avg_minus_ctl_lev39','rainnorm',landmaskxr,landlats,landlons) 
any_configuration_plot(outdir,runmin,runmax,-90.,90.,sphum_avg_ctl,area_array,'kg/kg','q_ctl_lev39','fromwhite',landmaskxr,landlats,landlons,nmb_contours=10) 
any_configuration_plot(outdir,runmin,runmax,-90.,90.,sphum_avg,area_array,'kg/kg','q_avg_lev39','fromwhite',landmaskxr,landlats,landlons,nmb_contours=10) 
any_configuration_plot(outdir,runmin,runmax,-90.,90.,M_ctl,area_array,'mm/day','M_ctl','fromwhite',landmaskxr,landlats,landlons)  
any_configuration_plot(outdir,runmin,runmax,-90.,90.,M_new,area_array,'mm/day','M_avg','fromwhite',landmaskxr,landlats,landlons) 
any_configuration_plot(outdir,runmin,runmax,-90.,90.,(Delta_M).mean('time')/M_ctl.mean('time'),area_array,' ','DeltaM_avg_div_M_ctl_avg','rainnorm',landmaskxr,landlats,landlons,minval=-1.0,maxval=1.0)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,((sphum_avg - sphum_avg_ctl)/sphum_avg_ctl),area_array,' ','Deltaq_div_q_ctl','rainnorm',landmaskxr,landlats,landlons)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(Delta_M).mean('time')*(sphum_avg - sphum_avg_ctl).mean('time'),area_array,'mm/day','DeltaM_avg_dot_Deltaq_avg','rainnorm',landmaskxr,landlats,landlons,minval=-2.0,maxval=2.0)
                                            
                     
###################################################################

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(net_lhe_avg - net_lhe_avg_ctl),area_array,'mm/day','E_avg_minus_ctl_oceanscale','rainnorm',landmaskxr,landlats,landlons,minval=0.,maxval=0.5)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(tsurf_avg-tsurf_avg_ctl),area_array,'K','$T_S$_avg_minus_ctl_narrowcbar','tempdiff',landmaskxr,landlats,landlons, minval = -4., maxval = 4.)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(tsurf_avg-tsurf_avg_ctl),area_array,'K','$T_S$_avg_minus_ctl','tempdiff',landmaskxr,landlats,landlons, minval = -6., maxval = 6.)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(net_sw_avg - net_sw_avg_ctl),area_array,'W/m2','SW_avg_minus_ctl','rainnorm',landmaskxr,landlats,landlons,minval=-4.,maxval=4.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,(net_lw_avg - net_lw_avg_ctl),area_array,'W/m2','LW_avg_minus_ctl','rainnorm',landmaskxr,landlats,landlons,minval=-40.,maxval=40.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,(net_t_avg - net_t_avg_ctl),area_array,'W/m2','SH_avg_minus_ctl','rainnorm',landmaskxr,landlats,landlons,minval=-10.,maxval=10.)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,((tsurf_avg-tsurf_avg_ctl)/tsurf_avg_ctl)*100.,area_array,'%','$T_S$_avg_minus_ctl_relativechange','tempdiff',landmaskxr,landlats,landlons)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,(((net_lhe_avg - net_lhe_avg_ctl)/net_lhe_avg_ctl)*100.),area_array,'%','lhe_flux_minus_ctl_oceanscale_relativechange','rainnorm',landmaskxr,landlats,landlons)

PE_avg=precipitation_avg-net_lhe_avg # 28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# # # see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC

PE_avg_sum = area_integral(PE_avg,area_array,landmaskxr,'all_sfcs',factor = 10**(-3)) # factor to 
#convert from mm/d to m/d
print('P avg - E avg global integral / total sfc area'+str(PE_avg_sum/total_sfc_area))

PE_avg_ctl=precipitation_avg_ctl-net_lhe_avg_ctl # 28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg

PE_seasonal_avg = precipitation_seasonal_avg - net_lhe_seasonal_avg
PE_seasonal_avg_ctl = precipitation_seasonal_avg_ctl - net_lhe_seasonal_avg_ctl


# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(bucket_depth_avg - bucket_depth_avg_ctl),area_array,'mm/day','bucket depth avg minus ctl','rainnorm',landmaskxr,landlats,landlons)
#any_configuration_plot(outdir,runmin,runmax,-90.,90.,(bucket_depth_avg - bucket_depth_avg_ctl).where(landmask==1.),area_array,'mm/day','bucket_depth_avg_minus_ctl_land','rainnorm',landmaskxr,landlats,landlons,minval = -0.1, maxval = 0.1)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,rh_avg - rh_avg_ctl,area_array,'%','surface_rh_avg_minus_ctl','rainnorm',landmaskxr,landlats,landlons,nmb_contours=5, minval=-7., maxval = 7.)

#any_configuration_plot(outdir,runmin,runmax,-90.,90.,CIWV_avg - CIWV_avg_ctl,area_array,'kg/kg','column_IWV_minus_ctl','rainnorm',landmaskxr,landlats,landlons,nmb_contours=5, minval = -0.05, maxval = 0.05)

# # degrees C symbol : ...,u"\u00b0"+'C',..

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,flux_oceanq_avg_ctl,area_array,'W/m^2','ocean_heat_transport_ctl','tempdiff',landmaskxr,landlats,landlons,minval=-200,maxval=200)

# # # #any_configuration_plot_minuszonavg(-90.,90.,tsurf_avg,'K','tsurf avg minus zonavg','temp','T avg')
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(PE_avg - PE_avg_ctl).where(landmask==1.),area_array,'mm/day','P-E avg minus ctl','rainnorm',landmaskxr,landlats,landlons,minval=-0.5,maxval=0.5
any_configuration_plot(outdir,runmin,runmax,-90.,90.,(PE_avg - PE_avg_ctl),area_array,'mm/day','P-E_avg_minus_ctl','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(PE_avg),area_array,'mm/day','P-E avg','rainnorm',landmaskxr,landlats,landlons,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(PE_avg_ctl),area_array,'mm/day','P-E avg ctl','rainnorm',landmaskxr,landlats,landlons,nmb_contours=4)


any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_avg,area_array,'mm/day','P_avg','fromwhite',landmaskxr,landlats,landlons,nmb_contours=10,minval = 0., maxval = 8.)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_avg.where(landmask==1.),area_array,'mm/day','P avg','fromwhite',landmaskxr,landlats,landlons,nmb_contours=4)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_avg_ctl,area_array,'mm/day','P_avg_ctl','fromwhite',landmaskxr,landlats,landlons,nmb_contours=10,minval = 0., maxval = 8.)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(precipitation_avg - precipitation_avg_ctl),area_array,'mm/day','P_avg_minus_ctl','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,((precipitation_avg - precipitation_avg_ctl)/precipitation_avg_ctl)*100.,area_array,'%','P_avg_minus_ctl_relativechange_','rainnorm',landmaskxr,landlats,landlons,minval=-100.,maxval=100.)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(precipitation_avg - precipitation_avg_ctl).where(landmask==1.),area_array,'mm/day','P_avg_minus_ctl_land','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,(net_lhe_avg - net_lhe_avg_ctl).where(landmask==1.),area_array,'mm/day','E_avg_minus_ctl_land','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,(net_lhe_avg - net_lhe_avg_ctl),area_array,'mm/day','E_avg_minus_ctl','rainnorm',landmaskxr,landlats,landlons,minval=0.,maxval=2.)

# [rh_ctl,rh_ctl_avg,rh_ctl_seasonal_avg,rh_ctl_month_avg,rh_ctl_annual_avg,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
# [rh,rh_avg,rh_seasonal_avg,rh_month_avg,rh_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'rh','%')

# for i in range(20,40):
#     any_configuration_plot(outdir,runmin,runmax,-90.,90.,rh_avg[i,:,:] - rh_ctl_avg[i,:,:],area_array,'%','lev'+str(i)+'_rh_avg_minus_ctl','rainnorm',landmaskxr,landlats,landlons,nmb_contours=5, minval=-7., maxval = 7.)

ctl_runmin_6hrly = 493
ctl_runmax_6hrly = 601

runmin_6hrly = 492
runmax_6hrly = 600

[sphum_6hrly,sphum_6hrly_avg,time]=seasonal_surface_variable_6hrly(testdir,model,runmin_6hrly,runmax_6hrly,'sphum','kg/kg',level=39)
[sphum_6hrly_ctl,sphum_6hrly_avg_ctl,time]=seasonal_surface_variable_6hrly(control_dir,ctl_model,ctl_runmin_6hrly,ctl_runmax_6hrly,'sphum','kg/kg',level=39)

[ucomp_6hrly,ucomp_6hrly_avg,time]=seasonal_surface_variable_6hrly(testdir,model,runmin_6hrly,runmax_6hrly,'ucomp','m/s',level=39)
[vcomp_6hrly,vcomp_6hrly_avg,time]=seasonal_surface_variable_6hrly(testdir,model,runmin_6hrly,runmax_6hrly,'vcomp','m/s',level=39)
[ucomp_6hrly_ctl,ucomp_6hrly_avg_ctl,time]=seasonal_surface_variable_6hrly(control_dir,ctl_model,ctl_runmin_6hrly,ctl_runmax_6hrly,'ucomp','m/s',level=39)
[vcomp_6hrly_ctl,vcomp_6hrly_avg_ctl,time]=seasonal_surface_variable_6hrly(control_dir,ctl_model,ctl_runmin_6hrly,ctl_runmax_6hrly,'vcomp','m/s',level=39)

[precipitation_6hrly,precipitation_6hrly_avg,time]=seasonal_surface_variable_6hrly(testdir,model,runmin_6hrly,runmax_6hrly,'precipitation','mm/d',factor=86400.)
[precipitation_6hrly_ctl,precipitation_6hrly_avg_ctl,time]=seasonal_surface_variable_6hrly(control_dir,ctl_model,ctl_runmin_6hrly,ctl_runmax_6hrly,'precipitation','mm/d',factor=86400.)


winds_one_level(outdir,runmin_6hrly,runmax_6hrly,'moisture_flux_total_',(ucomp_6hrly*sphum_6hrly).mean('time'),(vcomp_6hrly*sphum_6hrly).mean('time'),precipitation_6hrly_avg,'fromwhite','mm/d',0.,8.,landmaskxr,landlats,landlons,veclen=0.1,level=39)
winds_one_level(outdir,ctl_runmin_6hrly,ctl_runmax_6hrly,'moisture_flux_total_ctl_',(ucomp_6hrly_ctl*sphum_6hrly_ctl).mean('time'),(vcomp_6hrly_ctl*sphum_6hrly_ctl).mean('time'),precipitation_6hrly_avg_ctl,'fromwhite','mm/d',0.,8.,landmaskxr,landlats,landlons,veclen=0.1,level=39)

winds_one_level(outdir,runmin_6hrly,runmax_6hrly,'moisture_flux_stationary_',(ucomp_6hrly_avg*sphum_6hrly_avg),(vcomp_6hrly_avg*sphum_6hrly_avg),precipitation_6hrly_avg,'fromwhite','mm/d',0.,8.,landmaskxr,landlats,landlons,veclen=0.1,level=39)
winds_one_level(outdir,ctl_runmin_6hrly,ctl_runmax_6hrly,'moisture_flux_stationary_ctl_',(ucomp_6hrly_avg_ctl*sphum_6hrly_avg_ctl),(vcomp_6hrly_avg_ctl*sphum_6hrly_avg_ctl),precipitation_6hrly_avg_ctl,'fromwhite','mm/d',0.,8.,landmaskxr,landlats,landlons,veclen=0.1,level=39)





[ucomp,ucomp_avg,ucomp_seasonal_avg,ucomp_month_avg,ucomp_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[vcomp,vcomp_avg,vcomp_seasonal_avg,vcomp_month_avg,vcomp_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'vcomp','m/s')
[ucomp_ctl,ucomp_avg_ctl,ucomp_seasonal_avg_ctl,ucomp_month_avg_ctl,ucomp_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')
[vcomp_ctl,vcomp_avg_ctl,vcomp_seasonal_avg_ctl,vcomp_month_avg_ctl,vcomp_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'vcomp','m/s')

winds_seasons(ucomp_seasonal_avg,vcomp_seasonal_avg,39,(PE_seasonal_avg - PE_seasonal_avg_ctl),'rainnorm','mm/d',-2.,2.,landmaskxr,landlats,landlons,outdir,runmin,runmax)

winds_seasons((ucomp_seasonal_avg - ucomp_seasonal_avg_ctl),(vcomp_seasonal_avg - vcomp_seasonal_avg_ctl),39,(precipitation_seasonal_avg - precipitation_seasonal_avg_ctl),'rainnorm','mm/d',-2.,2.,landmaskxr,landlats,landlons,outdir,runmin,runmax,veclen=1.)

winds_at_heightlevel((ucomp_avg - ucomp_avg_ctl),(vcomp_avg - vcomp_avg_ctl),39,(precipitation_avg - precipitation_avg_ctl),'rainnorm','mm/d',-2.,2.,landmaskxr,landlats,landlons,veclen=1.)



################### Finding the level 'CH' were delta T over land and ocean is the same, (see Joshi et al., 2008) ###

[temp,temp_avg,temp_seasonal_avg,temp_month_avg,temp_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'temp','K')

[temp_ctl,temp_ctl_avg,temp_ctl_seasonal_avg,temp_ctl_month_avg,temp_ctl_annual_avg,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')


delta_temp_avg = temp_avg - temp_ctl_avg
# global mean
delta_temp_awave_ocean = xr.DataArray(area_weighted_avg_4D(delta_temp_avg,area_array_3D,landmaskxr,'ocean',minlat=-90.,maxlat=90.,axis=(1,2)))
delta_temp_awave_land = xr.DataArray(area_weighted_avg_4D(delta_temp_avg,area_array_3D,landmaskxr,'land',minlat=-90.,maxlat=90.,axis=(1,2)))

# tropical mean 
delta_temp_awave_ocean_tropics = xr.DataArray(area_weighted_avg_4D(delta_temp_avg,area_array_3D,landmaskxr,'ocean',minlat=-30.,maxlat=30.,axis=(1,2)))
delta_temp_awave_land_tropics = xr.DataArray(area_weighted_avg_4D(delta_temp_avg,area_array_3D,landmaskxr,'land',minlat=-30.,maxlat=30.,axis=(1,2)))

levs = temp_avg.pres_lev

fig, ax = plt.subplots(1,2,sharex = True, sharey=True,figsize=(25,10))
ax[0].plot(delta_temp_awave_land_tropics,levs,'r',label = 'Delta T land (tropics)')
ax[0].plot(delta_temp_awave_ocean_tropics,levs,'b',label = 'Delta T ocean (tropics)')

ax[1].plot(delta_temp_awave_land,levs,'r',label = 'Delta T land (global)')
ax[1].plot(delta_temp_awave_ocean,levs,'b',label = 'Delta T ocean (global)')

ax[0].set_xlabel("Delta T (K)")
ax[1].set_xlabel("Delta T (K)")
ax[0].tick_params()
ax[1].tick_params()
ax[0].legend()
ax[1].legend()
fig.gca().invert_yaxis()
ax[0].set_ylabel('Pressure (hPa)')

################################################################################################


# #[omega,omega_avg,omega_seasonal_avg,omega_month_avg,time]=seasonal_4D_variable(testdir,runmin,runmax,'omega','Pa/s')

# #[omega_ctl,omega_avg_ctl,omega_seasonal_avg_ctl,omega_month_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_runmin,ctl_runmax,'omega','Pa/s')


# winds_seasons((ucomp_seasonal_avg-ucomp_seasonal_avg_ctl),(vcomp_seasonal_avg-vcomp_seasonal_avg_ctl),39,(precipitation_seasonal_avg-precipitation_seasonal_avg_ctl),'rainnorm','mm/d',-2.,2.,landmaskxr,landlats,landlons,outdir,runmin,runmax)




#[convection_rain_ctl,convection_rain_avg_ctl,convection_rain_seasonal_avg_ctl,convection_rain_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'convection_rain','kg/m2s')
#[condensation_rain_ctl,condensation_rain_avg_ctl,condensation_rain_seasonal_avg_ctl,condensation_rain_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'condensation_rain','kg/m2s')



# maxval_precip = np.absolute((precipitation_month_avg*86400).max())

# for i in range(1,13):

#     month_plot = winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,precipitation_month_avg.sel(month=i)*86400,'fromwhite','mm/day',0,maxval_precip)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_precip'+str(i)+'.png',bbox_inches='tight')
# os.system('convert -delay 100 /scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_precip*.png /scratch/mp586/Code/Graphics/'+testdir+'/precip_wind_monthly_clim.gif')

# maxval_omega_surf = np.absolute((omega_month_avg[:,39,:,:]).max())
# minval_omega_surf = np.absolute((omega_month_avg[:,39,:,:]).min())


# for i in range(1,13):

#     month_plot = winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,(omega_month_avg[:,39,:,:]).sel(month=i),'rainnorm','mm/day',minval_omega_surf,maxval_omega_surf)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_omega'+str(i)+'.png',bbox_inches='tight')
# os.system('convert -delay 100 /scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_omega*.png /scratch/mp586/Code/Graphics/'+testdir+'/omega_wind_monthly_clim.gif')


# JJA = 'JJA'
# DJF = 'DJF'
# MAM = 'MAM'
# SON = 'SON'

# summer_plot = winds_at_heightlevel(ucomp_seasonal_avg.sel(season=JJA),vcomp_seasonal_avg.sel(season=JJA),39,precipitation_seasonal_avg.sel(season=JJA)*86400,'fromwhite','mm/day')
# summer_plot.savefig('anim_plot1.png',bbox_inches='tight')
# fall_plot = winds_at_heightlevel(ucomp_seasonal_avg.sel(season=SON),vcomp_seasonal_avg.sel(season=SON),39,precipitation_seasonal_avg.sel(season=SON)*86400,'fromwhite','mm/day')
# fall_plot.savefig('anim_plot2.png',bbox_inches='tight')
# winter_plot = winds_at_heightlevel(ucomp_seasonal_avg.sel(season=DJF),vcomp_seasonal_avg.sel(season=DJF),39,precipitation_seasonal_avg.sel(season=DJF)*86400,'fromwhite','mm/day')
# winter_plot.savefig('anim_plot3.png',bbox_inches='tight')
# spring_plot = winds_at_heightlevel(ucomp_seasonal_avg.sel(season=MAM),vcomp_seasonal_avg.sel(season=MAM),39,precipitation_seasonal_avg.sel(season=MAM)*86400,'fromwhite','mm/day')
# spring_plot.savefig('anim_plot4.png',bbox_inches='tight')

# os.system('convert -delay 50 anim_plot*.png animation_seasons.gif')



#animated_map(testdir,outdir,bucket_depth.where(landmask==1.),'m','bucket depth','bucket_depth','fromwhite',0,140) # need runmin = 1!






# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(omega_avg[39,:,:] - omega_avg_ctl[39,:,:]),'Pa/s','Omega avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,contourson = False)



# land_temp_global=tsurf_avg.where(landmask==1.).mean()
# ocean_temp_global=tsurf_avg.where(landmask==0.).mean()
# print('Average temperature over land (global) = '+str(land_temp_global))
# print('Average temperature over ocean (global) = '+str(ocean_temp_global))


# minlat=-30. #input('Enter minimum latitude ')
# maxlat=30. #input('Enter maximum latitude ')

# # difference between exp and control
# land_temp=(tsurf_avg-tsurf_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_temp=(tsurf_avg-tsurf_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!

# print('Average temperature diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_temp.mean()))
# print('Average temperature diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_temp.mean()))

# # difference between exp and control
# land_precip_global=(precipitation_avg-precipitation_avg_ctl).where(landmask==1.).mean() 
# ocean_precip_global=(precipitation_avg-precipitation_avg_ctl).where(landmask==0.).mean()
# print('Average precipitation diff over land (global) = '+str(land_precip_global*86400)+' mm/day')
# print('Average precipitation diff over ocean (global) = '+str(ocean_precip_global*86400)+' mm/day')

# # difference between exp and control
# land_precip=(precipitation_avg-precipitation_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_precip=(precipitation_avg-precipitation_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
# print('Average precipitation diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.mean()*86400)+'+/-'+str(land_precip.std()*86400)+'mm/day') # spatial standard deviation
# print('Average precipitation diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.mean()*86400)+'+/-'+str(ocean_precip.std()*86400)+'mm/day')
# print('Min precipitation diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.min()*86400)+'mm/day') # spatial standard deviation
# print('Max precipitation diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.max()*86400)+'mm/day') # spatial standard deviation
# print('Min precipitation diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.min()*86400)+'mm/day') # spatial standard deviation
# print('Max precipitation diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.max()*86400)+'mm/day')

# # difference between exp and control
# land_lhe=(net_lhe_avg-net_lhe_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_lhe=(net_lhe_avg-net_lhe_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
# print('Average net_lhe diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_lhe.mean()/28.)+'+/-'+str(land_lhe.std()/28.)+'mm/day') # spatial standard deviation
# print('Average net_lhe diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_lhe.mean()/28.)+'+/-'+str(ocean_lhe.std()/28.)+'mm/day')
# print('Min net_lhe diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_lhe.min()/28.)+'mm/day') # spatial standard deviation
# print('Max net_lhe diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_lhe.max()/28.)+'mm/day') # spatial standard deviation
# print('Min net_lhe diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_lhe.min()/28.)+'mm/day') # spatial standard deviation
# print('Max net_lhe diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_lhe.max()/28.)+'mm/day')




# JJA = 'JJA'
# DJF = 'DJF'
# MAM = 'MAM'
# SON = 'SON'

# land_precip_JJA=precipitation_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_precip_JJA=precipitation_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
# land_temp_JJA=tsurf_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_temp_JJA=tsurf_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!

# several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=JJA),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=JJA)*86400,'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=JJA)/28.,'E (mm/day)','g','JJA T_S, P and E') # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=DJF),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=DJF)*86400,'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=DJF)/28.,'E (mm/day)','g','DJF T_S, P and E')
# several_vars_zonalavg2(tsurf_avg,'T_S (K)','Red',precipitation_avg*86400,'P (mm/day)','Blue',net_lhe_avg/28.,'E (mm/day)','g','avg T_S, P and E')


# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_month_avg.sel(month=7)*86400,'mm/day','P_July (mm/day)','rainnorm',landmask,landlats,landlons)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_month_avg.sel(month=1)*86400,'mm/day','P_January (mm/day)','rainnorm',landmask,landlats,landlons)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,tsurf_month_avg.sel(month=7),'K','tsurf_July (K)','temp',landmask,landlats,landlons)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,tsurf_month_avg.sel(month=1),'K','tsurf_January (K)','temp',landmask,landlats,landlons)

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_month_avg.sel(month=7)*86400-net_lhe_month_avg.sel(month=7)/28.,'mm/day','P-E_July (mm/day)','rainnorm',landmask,landlats,landlons)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_month_avg.sel(month=1)*86400-net_lhe_month_avg.sel(month=1)/28.,'mm/day','P-E_January (mm/day)','rainnorm',landmask,landlats,landlons)


# print('JJA tsurf_avg (global) = '+str(tsurf_seasonal_avg.sel(season=JJA).mean()))
# print('DJF tsurf_avg (global) = '+str(tsurf_seasonal_avg.sel(season=DJF).mean()))
