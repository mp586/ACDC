from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
import plotting_routines
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
total_sfc_area = np.sum(area_array)
print ('total sfc area (*10^14) = '+str(np.sum(area_array/(10**14)))) # -- test: correct, equals sfc area of earth (5.1*10**14 m^2)
land_sfc_area = np.sum(area_array.where(landmask==1.))
print ('land sfc area (*10^14) = '+str(land_sfc_area/(10**14)))
ocean_sfc_area = np.sum(area_array.where(landmask!=1.))
print ('ocean sfc area (*10^14) = '+str(ocean_sfc_area/(10**14)))

# for plotting a spin up run ('control') timeseries followed by the timeseries from the perturbed experiment
globavg_var_timeseries_total_and_land_perturbed(testdir,model,area_array,'t_surf',1,runmax,1.,landmask,control_dir,ctl_model,1,ctl_timeseries_max)
# globavg_var_timeseries_total_and_land_perturbed(testdir,model,area_array,'bucket_depth',1,runmax,1.,landmask,control_dir,ctl_model,1,ctl_timeseries_max,select='land')


#### not working atm, see plotcalls any config. py
# [msf,msf_avg,msf_seasonal_avg,msf_month_avg] = mass_streamfunction(testdir,model,runmin,runmax) # 
# plot_streamfunction_seasonal(msf_seasonal_avg)

# globavg_var_timeseries(testdir,'co2',1,runmax)


[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','W/m^2',factor = 1/28.) # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','kg/m2s', factor=86400)
#[bucket_depth,bucket_depth_avg,bucket_depth_seasonal_avg,bucket_depth_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'bucket_depth','m')
[rh,rh_avg,rh_seasonal_avg,rh_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'rh','%',level=39)
#[sphum,sphum_avg,sphum_seasonal_avg,sphum_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level='all')

#[convection_rain,convection_rain_avg,convection_rain_seasonal_avg,convection_rain_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'convection_rain','kg/m2s')
#[condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,condensation_rain_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'condensation_rain','kg/m2s')



[tsurf_ctl,tsurf_avg_ctl,tsurf_seasonal_avg_ctl,tsurf_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[net_lhe_ctl,net_lhe_avg_ctl,net_lhe_seasonal_avg_ctl,net_lhe_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','W/m^2',factor = 1/28.) # latent heat flux at surface (UP)
[precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','kg/m2s', factor=86400)
#[bucket_depth_ctl,bucket_depth_avg_ctl,bucket_depth_seasonal_avg_ctl,bucket_depth_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'bucket_depth','m')
[rh_ctl,rh_ctl_avg,rh_ctl_seasonal_avg,rh_ctl_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%',level=39)
#[sphum_ctl,sphum_ctl_avg,sphum_ctl_seasonal_avg,sphum_ctl_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level='all')
#[flux_oceanq_ctl,flux_oceanq_avg_ctl,flux_oceanq_seasonal_avg_ctl,flux_oceanq_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_oceanq','W/m^2')

rh_P_E_change(rh_avg,rh_ctl_avg,precipitation_avg,precipitation_avg_ctl,net_lhe_avg,net_lhe_avg_ctl,landmask)






#[convection_rain_ctl,convection_rain_avg_ctl,convection_rain_seasonal_avg_ctl,convection_rain_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'convection_rain','kg/m2s')
#[condensation_rain_ctl,condensation_rain_avg_ctl,condensation_rain_seasonal_avg_ctl,condensation_rain_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'condensation_rain','kg/m2s')


# [ucomp,ucomp_avg,ucomp_seasonal_avg,ucomp_month_avg,time]=seasonal_4D_variable(testdir,runmin,runmax,'ucomp','m/s')
# [vcomp,vcomp_avg,vcomp_seasonal_avg,vcomp_month_avg,time]=seasonal_4D_variable(testdir,runmin,runmax,'vcomp','m/s')
# [omega,omega_avg,omega_seasonal_avg,omega_month_avg,time]=seasonal_4D_variable(testdir,runmin,runmax,'omega','Pa/s')


# [ucomp_ctl,ucomp_avg_ctl,ucomp_seasonal_avg_ctl,ucomp_month_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_runmin,ctl_runmax,'ucomp','m/s')
# [vcomp_ctl,vcomp_avg_ctl,vcomp_seasonal_avg_ctl,vcomp_month_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_runmin,ctl_runmax,'vcomp','m/s')
# [omega_ctl,omega_avg_ctl,omega_seasonal_avg_ctl,omega_month_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_runmin,ctl_runmax,'omega','Pa/s')



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



#animated_map(testdir,bucket_depth.where(landmask==1.),'m','bucket depth','bucket_depth','fromwhite',0,140) # need runmin = 1!


PE_avg=precipitation_avg-net_lhe_avg # 28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# # # see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC

PE_avg_sum = area_integral(PE_avg,area_array,landmaskxr,'all_sfcs',factor = 10**(-3)) # factor to 
#convert from mm/d to m/d
print('P avg - E avg global integral / total sfc area'+str(PE_avg_sum/total_sfc_area))

PE_avg_ctl=precipitation_avg_ctl-net_lhe_avg_ctl # 28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg

# any_configuration_plot(-90.,90.,(bucket_depth_avg - bucket_depth_avg_ctl),area_array,'mm/day','bucket depth avg minus ctrl','rainnorm',landmaskxr,landlats,landlons)
#any_configuration_plot(-90.,90.,(bucket_depth_avg - bucket_depth_avg_ctl).where(landmask==1.),area_array,'mm/day','bucket depth avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,minval = -0.1, maxval = 0.1)
any_configuration_plot(-90.,90.,rh_avg - rh_ctl_avg,area_array,'%','surface rh avg minus ctl','rainnorm',landmaskxr,landlats,landlons,nmb_contours=5, minval=-7., maxval = 7.)

any_configuration_plot(-90.,90.,sphum_avg - sphum_ctl_avg,area_array,'kg/kg','column IWV minus ctl','rainnorm',landmaskxr,landlats,landlons,nmb_contours=5, minval = -0.05, maxval = 0.05)

# # degrees C symbol : ...,u"\u00b0"+'C',...
any_configuration_plot(-90.,90.,(tsurf_avg-tsurf_avg_ctl),area_array,'K','$T_S$ avg minus ctrl','tempdiff',landmaskxr,landlats,landlons, minval = -6., maxval = 6.)

any_configuration_plot(-90.,90.,flux_oceanq_avg_ctl,area_array,'W/m^2','ocean heat transport ctl','tempdiff',landmaskxr,landlats,landlons,minval=-200,maxval=200)

# # # #any_configuration_plot_minuszonavg(-90.,90.,tsurf_avg,'K','tsurf avg minus zonavg','temp','T avg')
# any_configuration_plot(-90.,90.,(PE_avg - PE_avg_ctl).where(landmask==1.),area_array,'mm/day','P-E avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,minval=-0.5,maxval=0.5)
# any_configuration_plot(-90.,90.,(PE_avg - PE_avg_ctl),area_array,'mm/day','P-E avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)

# any_configuration_plot(-90.,90.,(PE_avg),area_array,'mm/day','P-E avg','rainnorm',landmaskxr,landlats,landlons,nmb_contours=4)
# any_configuration_plot(-90.,90.,(PE_avg_ctl),area_array,'mm/day','P-E avg ctl','rainnorm',landmaskxr,landlats,landlons,nmb_contours=4)


any_configuration_plot(-90.,90.,precipitation_avg,area_array,'mm/day','P avg','fromwhite',landmaskxr,landlats,landlons,nmb_contours=10,minval = 0., maxval = 8.)
# any_configuration_plot(-90.,90.,precipitation_avg.where(landmask==1.),area_array,'mm/day','P avg','fromwhite',landmaskxr,landlats,landlons,nmb_contours=4)
any_configuration_plot(-90.,90.,precipitation_avg_ctl,area_array,'mm/day','P avg ctl','fromwhite',landmaskxr,landlats,landlons,nmb_contours=10,minval = 0., maxval = 8.)

any_configuration_plot(-90.,90.,(precipitation_avg - precipitation_avg_ctl),area_array,'mm/day','P avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)
any_configuration_plot(-90.,90.,(precipitation_avg - precipitation_avg_ctl).where(landmask==1.),area_array,'mm/day','P avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)
any_configuration_plot(-90.,90.,(net_lhe_avg - net_lhe_avg_ctl).where(landmask==1.),area_array,'mm/day','E avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)

any_configuration_plot(-90.,90.,(net_lhe_avg - net_lhe_avg_ctl),area_array,'mm/day','E avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,minval=-2.,maxval=2.)



# any_configuration_plot(-90.,90.,(omega_avg[39,:,:] - omega_avg_ctl[39,:,:]),'Pa/s','Omega avg minus ctrl','rainnorm',landmaskxr,landlats,landlons,contourson = False)



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


# any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=7)*86400,'mm/day','P_July (mm/day)','rainnorm',landmask,landlats,landlons)
# any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=1)*86400,'mm/day','P_January (mm/day)','rainnorm',landmask,landlats,landlons)
# any_configuration_plot(-90.,90.,tsurf_month_avg.sel(month=7),'K','tsurf_July (K)','temp',landmask,landlats,landlons)
# any_configuration_plot(-90.,90.,tsurf_month_avg.sel(month=1),'K','tsurf_January (K)','temp',landmask,landlats,landlons)

# any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=7)*86400-net_lhe_month_avg.sel(month=7)/28.,'mm/day','P-E_July (mm/day)','rainnorm',landmask,landlats,landlons)
# any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=1)*86400-net_lhe_month_avg.sel(month=1)/28.,'mm/day','P-E_January (mm/day)','rainnorm',landmask,landlats,landlons)


# print('JJA tsurf_avg (global) = '+str(tsurf_seasonal_avg.sel(season=JJA).mean()))
# print('DJF tsurf_avg (global) = '+str(tsurf_seasonal_avg.sel(season=DJF).mean()))
