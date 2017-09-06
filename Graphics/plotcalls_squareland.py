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

testdir=input('Enter data directory name as string ')
runmin=97 #input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=241 #input('Enter runmax number ')

landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
landmask=landfile.variables['land_mask'][:]
lats=landfile.variables['lat'][:]
lons=landfile.variables['lon'][:]


# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'t_surf',1,runmax,1.,'true')
# # plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'coszen',1,runmax,1.,'true')

# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'precipitation',1,runmax,86400,'true')
# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'flux_lhe',1,runmax,1./28.,'true')
# plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'Orange','rh',1.,'m',39,1,runmax,'true')
# plotting_routines_kav7.tropics_severalvars_timeseries_oceanonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'Orange','rh',1.,'m',39,1,runmax,'true')
# #plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'Orange','bucket_depth',1.,'k',0,1,runmax,'true')
# plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'Orange','t_surf',1.,'r',0,1,runmax,'true')
# plotting_routines_kav7.tropics_severalvars_timeseries_oceanonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'Orange','t_surf',1.,'r',0,1,runmax,'true')
[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'t_surf','K')

[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'precipitation','kg/m2s')
[convection_rain,convection_rain_avg,convection_rain_seasonal_avg,convection_rain_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'convection_rain','kg/m2s')
[condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,condensation_rain_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'condensation_rain','kg/m2s')
#[bucket_depth,bucket_depth_avg,bucket_depth_seasonal_avg,bucket_depth_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'bucket_depth','m')

PE_avg=precipitation_avg*86400-net_lhe_avg/28. # 28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# # # see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC
# # plotting_routines_kav7.squareland_plot(-90.,90.,convection_rain_avg*86400,'mm/day','convection_rain avg','fromwhite')
# # plotting_routines_kav7.squareland_plot(-90.,90.,condensation_rain_avg*86400,'mm/day','condensation_rain avg','fromwhite')
plotting_routines_kav7.squareland_plot(-90.,90.,precipitation_avg.where(landmask==1.)*86400,'mm/day','P avg','fromwhite')
# # #plotting_routines_kav7.squareland_plot(-90.,90.,bucket_depth_avg.where(landmask==1.),'m','bucket_depth','fromwhite')
# # plotting_routines_kav7.squareland_plot(-90.,90.,net_lhe_avg.where(landmask==1.)/28.,'mm/day','E avg','fromwhite')
plotting_routines_kav7.squareland_plot(-100.,100.,PE_avg,'mm/day','P-E avg','rainnorm')
# # #plotting_routines_kav7.squareland_plot_minuszonavg(-90.,90.,PE_avg,'mm/day','P-E avg minus zonavg','rainnorm','P-E avg')
plotting_routines_kav7.squareland_plot(-90.,90.,tsurf_avg,'K','T_S avg','temp')
# # #plotting_routines_kav7.squareland_plot_minuszonavg(-90.,90.,tsurf_avg,'K','tsurf avg minus zonavg','temp','T avg')
plotting_routines_kav7.squareland_plot(-90.,90.,precipitation_avg*86400,'mm/day','P avg','fromwhite')
# # #plotting_routines_kav7.squareland_plot_minuszonavg(-90.,90.,precipitation_avg*86400,'mm/day','P avg minus zonavg','rainnorm','P avg')
# plotting_routines_kav7.squareland_plot(-90.,90.,net_lhe_avg/28.,'mm/day','E avg','fromwhite')

land_temp_global=tsurf_avg.where(landmask==1.).mean()
ocean_temp_global=tsurf_avg.where(landmask==0.).mean()
print('Average temperature over land (global) = '+str(land_temp_global))
print('Average temperature over ocean (global) = '+str(ocean_temp_global))

# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[lats,lons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

minlat=-30. #input('Enter minimum latitude ')
maxlat=30. #input('Enter maximum latitude ')

land_temp=tsurf_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
land_temp.plot()
plt.show()
ocean_temp=tsurf_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
ocean_temp.plot()
plt.show()
print('Average temperature over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_temp.mean()))
print('Average temperature over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_temp.mean()))


land_precip_global=precipitation_avg.where(landmask==1.).mean()
ocean_precip_global=precipitation_avg.where(landmask==0.).mean()
print('Average precipitation over land (global) = '+str(land_precip_global*86400)+' mm/day')
print('Average precipitation over ocean (global) = '+str(ocean_precip_global*86400)+' mm/day')


land_precip=precipitation_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_precip=precipitation_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
print('Average precipitation over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.mean()*86400)+'+/-'+str(land_precip.std()*86400)+'mm/day') # spatial standard deviation
print('Average precipitation over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.mean()*86400)+'+/-'+str(ocean_precip.std()*86400)+'mm/day')
print('Min precipitation over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.min()*86400)+'mm/day') # spatial standard deviation
print('Max precipitation over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.max()*86400)+'mm/day') # spatial standard deviation
print('Min precipitation over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.min()*86400)+'mm/day') # spatial standard deviation
print('Max precipitation over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.max()*86400)+'mm/day')

JJA = 'JJA'
DJF = 'DJF'
MAM = 'MAM'
SON = 'SON'

land_precip_JJA=precipitation_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_precip_JJA=precipitation_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
land_temp_JJA=tsurf_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_temp_JJA=tsurf_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!

plotting_routines_kav7.several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=JJA),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=JJA)*86400,'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=JJA)/28.,'E (mm/day)','Orange','JJA T_S, P and E') # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
plotting_routines_kav7.several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=DJF),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=DJF)*86400,'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=DJF)/28.,'E (mm/day)','Orange','DJF T_S, P and E')
plotting_routines_kav7.several_vars_zonalavg2(tsurf_avg,'T_S (K)','Red',precipitation_avg*86400,'P (mm/day)','Blue',net_lhe_avg/28.,'E (mm/day)','Orange','avg T_S, P and E')

[r1,pval1]=st.pattern_corr_2d(np.nan_to_num(tsurf_avg.where(landmask==0.)),np.nan_to_num(precipitation_avg.where(landmask==0.))) # correlation can not deal with nans
print('Spatial correlation of tsurf_avg_ocean and P_avg_ocean (global)='+str(r1))
[r1,pval1]=st.pattern_corr_2d(np.nan_to_num(tsurf_avg.where(landmask==1.)),np.nan_to_num(precipitation_avg.where(landmask==1.)))
print('Spatial correlation of tsurf_avg_land and P_avg_land (global)='+str(r1))
[r1,pval1]=st.pattern_corr_2d(np.nan_to_num(land_precip),np.nan_to_num(land_temp))
print('Spatial correlation between tsurf_avg_land and P_avg_land ('+str(minlat)+'N to '+str(maxlat)+'N) '+str(r1))
[r1,pval1]=st.pattern_corr_2d(np.nan_to_num(ocean_precip),np.nan_to_num(ocean_temp))
print('Spatial correlation between tsurf_avg_ocean and P_avg_ocean ('+str(minlat)+'N to '+str(maxlat)+'N) '+str(r1))


[r1,pval1]=st.pattern_corr_2d(np.nan_to_num(tsurf_seasonal_avg.sel(season=JJA).where(landmask==0.)),np.nan_to_num(precipitation_seasonal_avg.sel(season=JJA).where(landmask==0.)))
print('Spatial correlation of tsurf_ocean and P_ocean for JJA (global)='+str(r1))
[r1,pval1]=st.pattern_corr_2d(np.nan_to_num(tsurf_seasonal_avg.sel(season=JJA).where(landmask==1.)),np.nan_to_num(precipitation_seasonal_avg.sel(season=JJA).where(landmask==1.)))
print('Spatial correlation of tsurf_land and P_land for JJA (global)='+str(r1))
[r1,pval1]=st.pattern_corr_2d(np.nan_to_num(land_precip_JJA),np.nan_to_num(land_temp_JJA))
print('Spatial correlation between tsurf_JJA and P_JJA land ('+str(minlat)+'N to '+str(maxlat)+'N) '+str(r1))
[r1,pval1]=st.pattern_corr_2d(np.nan_to_num(ocean_precip_JJA),np.nan_to_num(ocean_temp_JJA))
print('Spatial correlation between tsurf_JJA and P_JJA ocean ('+str(minlat)+'N to '+str(maxlat)+'N) '+str(r1))

#plotting_routines_kav7.squareland_plot_correlation(-90.,90.,tsurf_avg,precipitation_avg,'tsurf vs precip')


plotting_routines_kav7.squareland_plot(-90.,90.,precipitation_month_avg.sel(month=7)*86400,'mm/day','P_July (mm/day)','rainnorm')
plotting_routines_kav7.squareland_plot(-90.,90.,precipitation_month_avg.sel(month=1)*86400,'mm/day','P_January (mm/day)','rainnorm')
plotting_routines_kav7.squareland_plot(-90.,90.,tsurf_month_avg.sel(month=7),'K','tsurf_July (K)','temp')
plotting_routines_kav7.squareland_plot(-90.,90.,tsurf_month_avg.sel(month=1),'K','tsurf_January (K)','temp')

plotting_routines_kav7.squareland_plot(-90.,90.,precipitation_month_avg.sel(month=7)*86400-net_lhe_month_avg.sel(month=7)/28.,'mm/day','P-E_July (mm/day)','rainnorm')
plotting_routines_kav7.squareland_plot(-90.,90.,precipitation_month_avg.sel(month=1)*86400-net_lhe_month_avg.sel(month=1)/28.,'mm/day','P-E_January (mm/day)','rainnorm')



print('JJA tsurf_avg (global) = '+str(tsurf_seasonal_avg.sel(season=JJA).mean()))
print('DJF tsurf_avg (global) = '+str(tsurf_seasonal_avg.sel(season=DJF).mean()))

# print('Antarctic January temp = '+str(tsurf_month_avg.sel(month=1,lat=slice(-90.,-60.)).mean()))
# print('Antarctic July temp = '+str(tsurf_month_avg.sel(month=7,lat=slice(-90.,-60.)).mean()))

# print('Arctic January temp = '+str(tsurf_month_avg.sel(month=1,lat=slice(60.,90.)).mean()))
# print('Arctic July temp = '+str(tsurf_month_avg.sel(month=7,lat=slice(60.,90.)).mean()))
