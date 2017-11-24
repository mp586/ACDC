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
import plotting_routines_kav7
import stats as st
sys.path.insert(0, '/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/src/extra/python/scripts') 
import cell_area as ca

# NOTE: Don't use landmaskxr as argument for where(....==1.), use np array landmask!!!!


testdir=input('Enter data directory name as string ')
runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number ')

# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/two_continents/land_two_continents.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/squareland/land_square.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/sqland_plus_antarctica/land_sqland_plus_antarctica.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/sqland_plus_antarctica/land_sqland_plus_antarctica_to35S.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/aquaplanet/land_aquaplanet.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/square_South_America/land_square_South_America.nc')
landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/all_continents/land.nc')


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

# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,area_array,'t_surf',1,runmax,1.,landmaskxr)
# # plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,area_array,'t_surf',1,runmax,1.,landmaskxr,minlat=-30.,maxlat=30.)
# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,area_array,'bucket_depth',1,runmax,1.,landmaskxr)





# still need to change to area_weighted zonal averages
# plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','bucket_depth',1.,'k',0,1,runmax,landmaskxr)
# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'precipitation',1,runmax,86400,landmask)
# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'flux_lhe',1,runmax,1./28.,landmask)
# plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','rh',1.,'m',39,1,runmax,landmaskxr)
# plotting_routines_kav7.tropics_severalvars_timeseries_oceanonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','rh',1.,'m',39,1,runmax,landmaskxr)
# plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','bucket_depth',1.,'k',0,1,runmax,landmaskxr)
# plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','t_surf',1.,'r',0,1,runmax,landmaskxr)
# plotting_routines_kav7.tropics_severalvars_timeseries_oceanonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','t_surf',1.,'r',0,1,runmax,landmaskxr)

[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'t_surf','K')
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'flux_lhe','W/m^2',factor = 1/28.) # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'precipitation','kg/m2s', factor=86400)
# [bucket_depth,bucket_depth_avg,bucket_depth_seasonal_avg,bucket_depth_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'bucket_depth','m')
[flux_oceanq,flux_oceanq_avg,flux_oceanq_seasonal_avg,flux_oceanq_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'flux_oceanq','W/m^2')




#[convection_rain,convection_rain_avg,convection_rain_seasonal_avg,convection_rain_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'convection_rain','kg/m2s')
#[condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,condensation_rain_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'condensation_rain','kg/m2s')

# [ucomp,ucomp_avg,ucomp_seasonal_avg,ucomp_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'ucomp','m/s')
# [vcomp,vcomp_avg,vcomp_seasonal_avg,vcomp_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'vcomp','m/s')
# [omega,omega_avg,omega_seasonal_avg,omega_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'omega','Pa/s')


# plotting_routines_kav7.any_configuration_plot(-90.,90.,(bucket_depth_avg),'m','bucket_depth','fromwhite',landmask,landlats,landlons)

# plotting_routines_kav7.animated_map(testdir,flux_oceanq_month_avg,'W/m^2','resulting_q_flux','qflux_clim_animated','rainnorm',0,12,-300,300)



# maxval_omega_surf = np.absolute((omega_month_avg[:,39,:,:]).max())
# minval_omega_surf = np.absolute((omega_month_avg[:,39,:,:]).min())


# for i in range(1,13):

#     month_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,(omega_month_avg[:,39,:,:]).sel(month=i),'rainnorm','dp/dt',minval_omega_surf,maxval_omega_surf)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_omega_run_'+str(runmin)+'-'+str(runmax)+'month_'+str(i)+'.png',bbox_inches='tight')
# os.system('convert -delay 100 /scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_omega_run_'+str(runmin)+'-'+str(runmax)+'*.png /scratch/mp586/Code/Graphics/'+testdir+'/omega_wind_monthly_clim_run_'+str(runmin)+'-'+str(runmax)+'.gif')


maxval_precip = np.absolute((precipitation_month_avg).max())

# for i in range(1,13):

#     month_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,precipitation_month_avg.sel(month=i),'fromwhite','mm/day',0,maxval_precip)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_precip_run_'+str(runmin)+'-'+str(runmax)+'month_'+str(i)+'.png',bbox_inches='tight')
# os.system('convert -delay 100 /scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_precip_run_'+str(runmin)+'-'+str(runmax)+'*.png /scratch/mp586/Code/Graphics/'+testdir+'/precip_wind_monthly_clim_run_'+str(runmin)+'-'+str(runmax)+'.gif')

# if runmin == 1:
#    plotting_routines_kav7.animated_map(testdir,bucket_depth,'m','bucket depth','bucket_depth','fromwhite',0,runmax-2,0,20)
#    plotting_routines_kav7.animated_map(testdir,tsurf,'m','tsurf','tsurf','temp',0,runmax-2,240,310)

PE_avg=precipitation_avg-net_lhe_avg # 28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg, see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC

PE_avg_sum = plotting_routines_kav7.area_integral(PE_avg,area_array,landmaskxr,'all_sfcs',factor = 10**(-3)) # factor to 
#convert from mm/d to m/d
print('P avg - E avg global integral / total sfc area'+str(PE_avg_sum/total_sfc_area))


###just a test -- gives same result
# PE_avg_sum = plotting_routines_kav7.area_integral(precipitation_avg,area_array,landmaskxr,'all_sfcs',factor = 10**(-3)) - plotting_routines_kav7.area_integral(net_lhe_avg,area_array,landmaskxr,'all_sfcs',factor = 10**(-3)) # factor to 
# #convert from mm/d to m/d
# print('P avg - E avg global integral / total sfc area'+str(PE_avg_sum/total_sfc_area))

PE_avg_sum = plotting_routines_kav7.area_integral(PE_avg,area_array,landmaskxr,'land',factor = 10**(-3)) # factor to 
#convert from mm/d to m/d
print('P avg - E avg global integral land / landsfc area'+str(PE_avg_sum/land_sfc_area))

# Need runmin = 1 for all of those 
# bd_0 = plotting_routines_kav7.area_integral(bucket_depth[0,:,:],area_array,landmaskxr,'all_sfcs')
# print('bd_0 allsfcs= '+str(bd_0/total_sfc_area))
# bd_end = plotting_routines_kav7.area_integral(bucket_depth[runmax-2,:,:],area_array,landmaskxr,'all_sfcs')
# print('bd_end allsfcs= '+str(bd_end/total_sfc_area))

# bd_0 = plotting_routines_kav7.area_integral(bucket_depth[0,:,:],area_array,landmaskxr,'land')
# print('bd_0 land= '+str(bd_0/land_sfc_area))
# bd_end = plotting_routines_kav7.area_integral(bucket_depth[runmax-2,:,:],area_array,landmaskxr,'land')
# print('bd_end land= '+str(bd_end/land_sfc_area))

#plotting_routines_kav7.any_configuration_plot(-90.,90.,bucket_depth_avg.where(landmask==1.),'m','bucket_depth','fromwhite',landmask,landlats,landlons,minval=0.,maxval=2.)

# plotting_routines_kav7.any_configuration_plot(-100.,100.,(bucket_depth[runmax-2,:,:] - bucket_depth[11,:,:]),area_array,'m','bucket depth dec year 40 - year 1','rainnorm',landmaskxr,landlats,landlons)

# plotting_routines_kav7.any_configuration_plot(-100.,100.,(tsurf[runmax-2,:,:] - tsurf[11,:,:]),area_array,'m','tsurf dec year 40 - year 1','tempdiff',landmaskxr,landlats,landlons)





plotting_routines_kav7.any_configuration_plot(-100.,100.,PE_avg,area_array,'mm/day','P-E avg','rainnorm',landmaskxr,landlats,landlons,minval=-12.,maxval=12.)
plotting_routines_kav7.any_configuration_plot(-90.,90.,tsurf_avg,area_array,'K','T_S avg','temp',landmaskxr,landlats,landlons,minval=240.,maxval=310.)
plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_avg,area_array,'mm/day','P avg','fromwhite',landmaskxr,landlats,landlons,contourson=True,minval=0.,maxval=10.)
plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_avg.where(landmask==1.),area_array,'mm/day','P avg','fromwhite',landmaskxr,landlats,landlons,contourson=True,minval=0., maxval = 5.)
plotting_routines_kav7.any_configuration_plot(-90.,90.,net_lhe_avg,area_array,'mm/day','E avg','fromwhite',landmaskxr,landlats,landlons)

# # plotting_routines_kav7.any_configuration_plot(-90.,90.,convection_rain_avg,'mm/day','convection_rain avg','fromwhite')
# # plotting_routines_kav7.any_configuration_plot(-90.,90.,condensation_rain_avg,'mm/day','condensation_rain avg','fromwhite')
#plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_avg.where(landmaskxr==1.),'mm/day','P avg','fromwhite')
# # plotting_routines_kav7.any_configuration_plot(-90.,90.,net_lhe_avg.where(landmaskxr==1.),'mm/day','E avg','fromwhite')


globavg_tsurf_w = plotting_routines_kav7.area_weighted_avg(tsurf_avg,area_array,landmask,'all_sfcs')
print('global_temp_unweighted = '+str(tsurf_avg.mean()))
print('global_temp_weighted = '+str(globavg_tsurf_w))

land_temp_globavg_w = plotting_routines_kav7.area_weighted_avg(tsurf_avg,area_array,landmask,'land')
ocean_temp_globavg_w = plotting_routines_kav7.area_weighted_avg(tsurf_avg,area_array,landmask,'ocean')
print('Area weighted Average temperature over land (global) = '+str(land_temp_globavg_w))
print('Area weighted Average temperature over ocean (global) = '+str(ocean_temp_globavg_w))

globavg_precipitation_w = plotting_routines_kav7.area_weighted_avg(precipitation_avg,area_array,landmask,'all_sfcs')
print('global_precip_unweighted = '+str(precipitation_avg.mean()))
print('global_precip_weighted = '+str(globavg_precipitation_w))

land_precip_globavg_w = plotting_routines_kav7.area_weighted_avg(precipitation_avg,area_array,landmask,'land')
ocean_precip_globavg_w = plotting_routines_kav7.area_weighted_avg(precipitation_avg,area_array,landmask,'ocean')
print('Area weighted Average preciperature over land (global) = '+str(land_precip_globavg_w))
print('Area weighted Average preciperature over ocean (global) = '+str(ocean_precip_globavg_w))





# still need to make those the weighted averages 
minlat=-30. #input('Enter minimum latitude ')
maxlat=30. #input('Enter maximum latitude ')

# land_temp=tsurf_avg_w.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_temp=tsurf_avg_w.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
# print('Average temperature over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_temp.sum()))
# print('Average temperature over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_temp.sum()))


# land_precip=precipitation_avg_w.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_precip=precipitation_avg_w.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
# print('Average precipitation over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.sum())+'+/-'+str(land_precip.std())+'mm/day') # spatial standard deviation
# print('Average precipitation over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.sum())+'+/-'+str(ocean_precip.std())+'mm/day')
# print('Min precipitation over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.min())+'mm/day') # spatial standard deviation
# print('Max precipitation over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.max())+'mm/day') # spatial standard deviation
# print('Min precipitation over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.min())+'mm/day') # spatial standard deviation
# print('Max precipitation over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.max())+'mm/day')

# JJA = 'JJA'
# DJF = 'DJF'
# MAM = 'MAM'
# SON = 'SON'

# land_precip_JJA=precipitation_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_precip_JJA=precipitation_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
# land_temp_JJA=tsurf_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
# ocean_temp_JJA=tsurf_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!

# plotting_routines_kav7.several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=JJA),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=JJA),'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=JJA),'E (mm/day)','g','JJA T_S, P and E') # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# plotting_routines_kav7.several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=DJF),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=DJF),'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=DJF),'E (mm/day)','g','DJF T_S, P and E')
# plotting_routines_kav7.several_vars_zonalavg2(tsurf_avg,'T_S (K)','Red',precipitation_avg,'P (mm/day)','Blue',net_lhe_avg,'E (mm/day)','g','avg T_S, P and E')

# [r1,pval1]=st.pattern_corr_2d(np.nan_to_num(tsurf_avg.where(landmask==0.)),np.nan_to_num(precipitation_avg.where(landmask==0.))) # correlation can not deal with nans
# print('Spatial correlation of tsurf_avg_ocean and P_avg_ocean (global)='+str(r1))
# [r1,pval1]=st.pattern_corr_2d(np.nan_to_num(tsurf_avg.where(landmask==1.)),np.nan_to_num(precipitation_avg.where(landmask==1.)))
# print('Spatial correlation of tsurf_avg_land and P_avg_land (global)='+str(r1))
# [r1,pval1]=st.pattern_corr_2d(np.nan_to_num(land_precip),np.nan_to_num(land_temp))
# print('Spatial correlation between tsurf_avg_land and P_avg_land ('+str(minlat)+'N to '+str(maxlat)+'N) '+str(r1))
# [r1,pval1]=st.pattern_corr_2d(np.nan_to_num(ocean_precip),np.nan_to_num(ocean_temp))
# print('Spatial correlation between tsurf_avg_ocean and P_avg_ocean ('+str(minlat)+'N to '+str(maxlat)+'N) '+str(r1))


# [r1,pval1]=st.pattern_corr_2d(np.nan_to_num(tsurf_seasonal_avg.sel(season=JJA).where(landmask==0.)),np.nan_to_num(precipitation_seasonal_avg.sel(season=JJA).where(landmask==0.)))
# print('Spatial correlation of tsurf_ocean and P_ocean for JJA (global)='+str(r1))
# [r1,pval1]=st.pattern_corr_2d(np.nan_to_num(tsurf_seasonal_avg.sel(season=JJA).where(landmask==1.)),np.nan_to_num(precipitation_seasonal_avg.sel(season=JJA).where(landmask==1.)))
# print('Spatial correlation of tsurf_land and P_land for JJA (global)='+str(r1))
# [r1,pval1]=st.pattern_corr_2d(np.nan_to_num(land_precip_JJA),np.nan_to_num(land_temp_JJA))
# print('Spatial correlation between tsurf_JJA and P_JJA land ('+str(minlat)+'N to '+str(maxlat)+'N) '+str(r1))
# [r1,pval1]=st.pattern_corr_2d(np.nan_to_num(ocean_precip_JJA),np.nan_to_num(ocean_temp_JJA))
# print('Spatial correlation between tsurf_JJA and P_JJA ocean ('+str(minlat)+'N to '+str(maxlat)+'N) '+str(r1))

#plotting_routines_kav7.any_configuration_plot_correlation(-90.,90.,tsurf_avg,precipitation_avg,'tsurf vs precip')


plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=7),area_array,'mm/day','P_July (mm/day)','fromwhite',landmaskxr,landlats,landlons,contourson=True)
plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=1),area_array,'mm/day','P_January (mm/day)','fromwhite',landmaskxr,landlats,landlons,contourson=True)
#plotting_routines_kav7.any_configuration_plot(-90.,90.,tsurf_month_avg.sel(month=7),'K','tsurf_July (K)','temp',landmask,landlats,landlons)
#plotting_routines_kav7.any_configuration_plot(-90.,90.,tsurf_month_avg.sel(month=1),'K','tsurf_January (K)','temp',landmask,landlats,landlons)

#plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=7)-net_lhe_month_avg.sel(month=7),'mm/day','P-E_July (mm/day)','rainnorm',landmask,landlats,landlons)
#plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=1)-net_lhe_month_avg.sel(month=1),'mm/day','P-E_January (mm/day)','rainnorm',landmask,landlats,landlons)







