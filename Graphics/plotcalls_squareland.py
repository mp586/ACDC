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

testdir= input('Enter data directory name as string ')
runmin=97 #input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=241 #input('Enter runmax number ')

landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
landmask=landfile.variables['land_mask'][:]
lats=landfile.variables['lat'][:]
lons=landfile.variables['lon'][:]

#plotting_routines_kav7.globavg_var_timeseries(testdir,'t_surf',109,122)
#plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'t_surf',1,runmax,1.,'true')
# plotting_routines_kav7.globavg_var_timeseries(testdir,'co2',1,runmax)

# # plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'coszen',1,runmax,1.,'true')

# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'precipitation',1,runmax,86400,'true')
# plotting_routines_kav7.globavg_var_timeseries_total_and_land(testdir,'flux_lhe',1,runmax,1./28.,'true')
# plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','rh',1.,'m',39,1,runmax,'true')
# plotting_routines_kav7.tropics_severalvars_timeseries_oceanonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','rh',1.,'m',39,1,runmax,'true')
# #plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','bucket_depth',1.,'k',0,1,runmax,'true')
# plotting_routines_kav7.tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','t_surf',1.,'r',0,1,runmax,'true')
# plotting_routines_kav7.tropics_severalvars_timeseries_oceanonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','t_surf',1.,'r',0,1,runmax,'true')

[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'t_surf','K')
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'precipitation','kg/m2s')
[convection_rain,convection_rain_avg,convection_rain_seasonal_avg,convection_rain_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'convection_rain','kg/m2s')
[condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,condensation_rain_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'condensation_rain','kg/m2s')
[bucket_depth,bucket_depth_avg,bucket_depth_seasonal_avg,bucket_depth_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'bucket_depth','m')
[ucomp,ucomp_avg,ucomp_seasonal_avg,ucomp_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'ucomp','m/s')
[vcomp,vcomp_avg,vcomp_seasonal_avg,vcomp_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'vcomp','m/s')
[omega,omega_avg,omega_seasonal_avg,omega_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'omega','Pa/s')


# JJA = 'JJA'
# DJF = 'DJF'
# MAM = 'MAM'
# SON = 'SON'

# summer_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_seasonal_avg.sel(season=JJA),vcomp_seasonal_avg.sel(season=JJA),39,precipitation_seasonal_avg.sel(season=JJA)*86400,'fromwhite','mm/day')
# summer_plot.savefig('anim_plot1.png',bbox_inches='tight')
# fall_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_seasonal_avg.sel(season=SON),vcomp_seasonal_avg.sel(season=SON),39,precipitation_seasonal_avg.sel(season=SON)*86400,'fromwhite','mm/day')
# fall_plot.savefig('anim_plot2.png',bbox_inches='tight')
# winter_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_seasonal_avg.sel(season=DJF),vcomp_seasonal_avg.sel(season=DJF),39,precipitation_seasonal_avg.sel(season=DJF)*86400,'fromwhite','mm/day')
# winter_plot.savefig('anim_plot3.png',bbox_inches='tight')
# spring_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_seasonal_avg.sel(season=MAM),vcomp_seasonal_avg.sel(season=MAM),39,precipitation_seasonal_avg.sel(season=MAM)*86400,'fromwhite','mm/day')
# spring_plot.savefig('anim_plot4.png',bbox_inches='tight')

# os.system('convert -delay 50 anim_plot*.png animation_seasons.gif')

maxval_precip = np.absolute((precipitation_month_avg*86400).max())

for i in range(1,13):

    month_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,precipitation_month_avg.sel(month=i)*86400,'fromwhite','mm/day',0,maxval_precip)
    month_plot.savefig('anim_plot'+str(i)+'.png',bbox_inches='tight')
os.system('convert -delay 100 anim_plot*.png precip_wind_monthly_clim.gif')

maxval_omega_surf = np.absolute((omega_month_avg[:,39,:,:]).max())
minval_omega_surf = np.absolute((omega_month_avg[:,39,:,:]).min())


for i in range(1,13):

    month_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,(omega_month_avg[:,39,:,:]).sel(month=i),'rainnorm','mm/day',minval_omega_surf,maxval_omega_surf)
    month_plot.savefig('anim_plot_omega'+str(i)+'.png',bbox_inches='tight')
os.system('convert -delay 100 anim_plot_omega*.png omega_wind_monthly_clim.gif')

#precipitation = convection_rain+condensation_rain
#precipitation_avg=convection_rain_avg+condensation_rain_avg

#plotting_routines_kav7.animated_map(testdir,bucket_depth.where(landmask==1.),'m','bucket depth','bucket_depth','fromwhite',0,140) # need runmin = 1!


PE_avg=precipitation_avg*86400-net_lhe_avg/28. # 28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# # # see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC


# # plotting_routines_kav7.squareland_plot(-90.,90.,convection_rain_avg*86400,'mm/day','convection_rain avg','fromwhite')
# # plotting_routines_kav7.squareland_plot(-90.,90.,condensation_rain_avg*86400,'mm/day','condensation_rain avg','fromwhite')
#plotting_routines_kav7.squareland_plot(-90.,90.,precipitation_avg.where(landmask==1.)*86400,'mm/day','P avg','fromwhite')
# # #plotting_routines_kav7.squareland_plot(-90.,90.,bucket_depth_avg.where(landmask==1.),'m','bucket_depth','fromwhite')
# # plotting_routines_kav7.squareland_plot(-90.,90.,net_lhe_avg.where(landmask==1.)/28.,'mm/day','E avg','fromwhite')
plotting_routines_kav7.squareland_plot(-100.,100.,PE_avg,'mm/day','P-E avg','rainnorm')
# # #plotting_routines_kav7.squareland_plot_minuszonavg(-90.,90.,PE_avg,'mm/day','P-E avg minus zonavg','rainnorm','P-E avg')
plotting_routines_kav7.squareland_plot(-90.,90.,tsurf_avg,'K','$T_S$ avg','temp') # degrees C symbol : ...,u"\u00b0"+'C',...

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
land_PE=PE_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
land_PE = land_PE.sel(lon = slice(0.,50.))
print('Average P-E over land between '+str(minlat)+'N and '+str(maxlat)+'N , and 0. - 50. E= '+str(land_PE.mean())+'+/-'+str(land_PE.std())+'mm/day') # spatial standard deviation

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

plotting_routines_kav7.several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=JJA),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=JJA)*86400,'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=JJA)/28.,'E (mm/day)','g','JJA T_S, P and E') # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
plotting_routines_kav7.several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=DJF),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=DJF)*86400,'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=DJF)/28.,'E (mm/day)','g','DJF T_S, P and E')
plotting_routines_kav7.several_vars_zonalavg2(tsurf_avg,'T_S (K)','Red',precipitation_avg*86400,'P (mm/day)','Blue',net_lhe_avg/28.,'E (mm/day)','g','avg T_S, P and E')

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
