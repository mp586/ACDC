from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from mpl_toolkits.basemap import interp

import xarray as xr
import pandas as pd
import os
import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
# import plotting_routines
from plotting_routines_kav7 import * # isca and gfdl have 0:04 and 0:03 
# filename format, respectively --> choose correct plotting routines kav7 or kav7_isca .py

import stats as st

GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts')) 
import cell_area as ca
# NOTE: Don't use landmaskxr as argument for where(....==1.), use np array landmask!!!!

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


landfile=Dataset(os.path.join(GFDL_BASE,'input/'+input('Which landmask? ')+'/land.nc'),mode='r')


landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]

# for specified lats

landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

area_array = ca.cell_area(t_res=42,base_dir='/scratch/mp586/Isca/')
area_array = xr.DataArray(area_array)

area_array_1deg = ca.cell_area(t_res='1_deg',base_dir='/scratch/mp586/Isca/')
area_array_1deg = xr.DataArray(area_array_1deg)

total_sfc_area = np.sum(area_array)
print ('total sfc area (*10^14) = '+str(total_sfc_area/(10**14))) # -- test: correct, equals sfc area of earth (5.1*10**14 m^2)
land_sfc_area = np.sum(area_array.where(landmask==1.))
print ('land sfc area (*10^14) = '+str(land_sfc_area/(10**14)))
ocean_sfc_area = np.sum(area_array.where(landmask!=1.))
print ('ocean sfc area (*10^14) = '+str(ocean_sfc_area/(10**14)))

globavg_var_timeseries_total_and_land(outdir,testdir,model,area_array,'t_surf',runmin,runmax,1.,landmaskxr,select='all')


### only for 6 hourly data ###
# globavg_var_timeseries_total_and_land_6hrly_severalvars(outdir,testdir,model,area_array,'flux_lhe','precipitation','t_surf','bucket_depth',runmin,runmax,1/28.,86400.,1.,1.,landmaskxr,minlat=-10.,maxlat=10.,select='land')


########### Mass Stream function plots 
# not working atm: File "Graphics/plotcalls_any_config.py", line 72, in <module>
  #   [msf,msf_avg,msf_seasonal_avg,msf_month_avg] = mass_streamfunction(testdir,model,runmin,runmax) # 
  # File "/scratch/mp586/Code/PYCODES/py", line 2490, in mass_streamfunction
  #   msf[i-runmin] = (np.cumsum(vbar*dp, axis='pfull')*c)
# TypeError: 'str' object cannot be interpreted as an index

[msf,msf_avg,msf_seasonal_avg,msf_month_avg] = mass_streamfunction(testdir,model,runmin,runmax) # 
plot_streamfunction_seasonal(msf_seasonal_avg)
###########

[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','kg/m2s', factor=86400)
#[bucket_depth,bucket_depth_avg,bucket_depth_seasonal_avg,bucket_depth_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'bucket_depth','m')

# [flux_oceanq,flux_oceanq_avg,flux_oceanq_seasonal_avg,flux_oceanq_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_oceanq','W/m^2')
[net_sw,net_sw_avg,net_sw_seasonal_avg,net_sw_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_sw','W/m^2',factor = 1.) # 
[lw_down,lw_down_avg,lw_down_seasonal_avg,lw_down_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lw','W/m^2',factor = 1.) # 
[net_t,net_t_avg,net_t_seasonal_avg,net_t_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_t','W/m^2',factor = 1.) # 
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','W/m^2',factor = 1.) #for SEB in W/m²

sigma = 5.67*10**(-8)
net_lw_avg = sigma*(tsurf_avg**4) - lw_down_avg

SEB = area_weighted_avg(net_sw_avg,area_array,landmaskxr,option='all_sfcs') - area_weighted_avg(net_lw_avg,area_array,landmaskxr,option='all_sfcs')- area_weighted_avg(net_lhe_avg,area_array,landmaskxr,option='all_sfcs') - area_weighted_avg(net_t_avg,area_array,landmaskxr,option='all_sfcs') # optional, if there is a qflux + area_weighted_avg(flux_oceanq_avg,area_array,landmaskxr,option='ocean')

print('SEB in W/m2 = '+str(SEB))

[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','W/m^2',factor = 1/28.) # latent heat flux at surface (UP) in mm/day
# 1/28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg, see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC


# only possible for full continents newbucket fullnewbucketqflux 
# runmin=481
# runmax=493
# globavg_var_timeseries_selected_points__6hrly_severalvars(outdir,testdir,model,area_array,'flux_lhe','precipitation','t_surf','bucket_depth','rh',runmin,runmax,1/28.,86400.,1.,1.,1.,landmask,precipitation_avg,minlat=-10.,maxlat=10.,maxormin='max')

# globavg_var_timeseries_selected_points__6hrly_severalvars(outdir,testdir,model,area_array,'flux_lhe','precipitation','t_surf','bucket_depth','rh',runmin,runmax,1/28.,86400.,1.,1.,1.,landmask,precipitation_avg,minlat=-10.,maxlat=10.,maxormin='min')




[slp,slp_avg,slp_seasonal_avg,slp_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'slp','hPa', factor = 1/100.)

[rh,rh_avg,rh_seasonal_avg,rh_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'rh','%',level=39)
# [sphum,sphum_avg,sphum_seasonal_avg,sphum_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level='all')

rh_P_E_T(outdir,runmin,runmax,rh_avg,precipitation_avg,net_lhe_avg,tsurf_avg,landmask)




# animated_map(outdir,slp_month_avg,'hPa','slp','slp_clim_animated','slp',0,12)
# animated_map(outdir,bucket_depth_month_avg,'m','soil moisture','bucket_depth_animated','fromwhite',0,12,minval = 0., maxval = 0.15)
# animated_map(outdir,tsurf_month_avg,'K','T','tsurf_animated','temp',0,12,minval = 240.,maxval = 320.)
# animated_map(outdir,net_lhe_month_avg,'mm/d','E','E_animated','fromwhite',0,12,minval = 0.,maxval = 8.0)
# animated_map(outdir,precipitation_month_avg,'mm/d','P','P_animated','fromwhite',0,12,minval = 0.,maxval = 8.)



any_configuration_plot(outdir,runmin,runmax,-100.,100.,rh_avg,area_array,'%','rh_avg','fromwhite',landmaskxr,minval = 0, maxval = 100, nmb_contours=5)
# any_configuration_plot(outdir,runmin,runmax,-100.,100.,sphum_avg,area_array,'kg/kg','column int WV','fromwhite',landmaskxr,nmb_contours=10)



# any_configuration_plot(outdir,runmin,runmax,outdir,runmin,runmax,-90.,90.,bucket_depth_avg.where(landmask==1.),area_array,'m','bucket_depth','fromwhite',landmaskxr,minval=0.,maxval=.5)
# if runmin == 1:
#    animated_map(testdir,bucket_depth,'m','bucket depth','bucket_depth','fromwhite',0,runmax-2,0,2)
#    animated_map(testdir,tsurf,'m','tsurf','tsurf','temp',0,runmax-2,240,310)


globavg_tsurf_w = area_weighted_avg(tsurf_avg,area_array,landmask,'all_sfcs')
print('global_temp_unweighted = '+str(tsurf_avg.mean()))
print('global_temp_weighted = '+str(globavg_tsurf_w))

land_temp_globavg_w = area_weighted_avg(tsurf_avg,area_array,landmask,'land')
ocean_temp_globavg_w = area_weighted_avg(tsurf_avg,area_array,landmask,'ocean')
print('Area weighted Average temperature over land (global) = '+str(land_temp_globavg_w))
print('Area weighted Average temperature over ocean (global) = '+str(ocean_temp_globavg_w))

globavg_precipitation_w = area_weighted_avg(precipitation_avg,area_array,landmask,'all_sfcs')
print('global_precip_unweighted = '+str(precipitation_avg.mean()))
print('global_precip_weighted = '+str(globavg_precipitation_w))

land_precip_globavg_w = area_weighted_avg(precipitation_avg,area_array,landmask,'land')
ocean_precip_globavg_w = area_weighted_avg(precipitation_avg,area_array,landmask,'ocean')
print('Area weighted Average precip over land (global) = '+str(land_precip_globavg_w))
print('Area weighted Average precip over ocean (global) = '+str(ocean_precip_globavg_w))



PE_avg=precipitation_avg-net_lhe_avg 

# any_configuration_plot(outdir,runmin,runmax,-100.,100.,(PE_avg).where(landmask==1.),area_array,'mm/day','P-E avg','rainnorm',landmaskxr,minval=-2.,maxval=2.)

any_configuration_plot(outdir,runmin,runmax,-100.,100.,(PE_avg),area_array,'mm/day','P-E_avg','rainnorm',landmaskxr,minval=-3.,maxval=3.)

# # any_configuration_plot(outdir,runmin,runmax,-90.,90.,net_lhe_avg.where(landmask==1.),area_array,'mm/day','E avg','fromwhite',landmaskxr,nmb_contours=4, minval = 0., maxval = 8.)

#any_configuration_plot(outdir,runmin,runmax,-90.,90.,flux_oceanq_avg,area_array,'W/m^2','ocean_heat_transport','tempdiff',landmaskxr,minval=-200,maxval=200)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,net_lhe_avg,area_array,'mm/day','E_avg','fromwhite',landmaskxr,nmb_contours=4,minval = 0., maxval = 8.)

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,slp_avg,area_array,'hPa','slp avg','slp',landmaskxr,nmb_contours=10)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_avg,area_array,'mm/day','P_avg','fromwhite',landmaskxr,nmb_contours=8,minval=0.,maxval=8.)
any_configuration_plot(outdir,runmin,runmax,-90.,90.,tsurf_avg,area_array,'K','avg_surface_T','temp',landmaskxr,nmb_contours=5)

any_configuration_plot(outdir,runmin,runmax,-90.,90.,tsurf_avg.where(landmask==1.),area_array,'K','avg_surface_T_land','temp',landmaskxr,nmb_contours=5)


JJA = 'JJA'
DJF = 'DJF'
MAM = 'MAM'
SON = 'SON'

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_seasonal_avg.sel(season=MAM),area_array,'mm/day','P_MAM (mm/day)','fromwhite',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_seasonal_avg.sel(season=SON),area_array,'mm/day','P_SON (mm/day)','fromwhite',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,net_lhe_seasonal_avg.sel(season=MAM),area_array,'mm/day','E_MAM (mm/day)','fromwhite',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,net_lhe_seasonal_avg.sel(season=SON),area_array,'mm/day','E_SON (mm/day)','fromwhite',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,tsurf_seasonal_avg.sel(season=MAM),area_array,'K','tsurf_MAM (K)','temp',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,tsurf_seasonal_avg.sel(season=SON),area_array,'K','tsurf_SON (K)','temp',landmaskxr,nmb_contours=4)




# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_seasonal_avg.sel(season=JJA),area_array,'mm/day','P_JJA (mm/day)','fromwhite',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_seasonal_avg.sel(season=DJF),area_array,'mm/day','P_DJF (mm/day)','fromwhite',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,net_lhe_seasonal_avg.sel(season=JJA),area_array,'mm/day','E_JJA (mm/day)','fromwhite',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,net_lhe_seasonal_avg.sel(season=DJF),area_array,'mm/day','E_DJF (mm/day)','fromwhite',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,tsurf_seasonal_avg.sel(season=JJA),area_array,'K','tsurf_JJA (K)','temp',landmaskxr,nmb_contours=4)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,tsurf_seasonal_avg.sel(season=DJF),area_array,'K','tsurf_DJF (K)','temp',landmaskxr,nmb_contours=4)


[ucomp,ucomp_avg,ucomp_seasonal_avg,ucomp_month_avg,ucomp_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[vcomp,vcomp_avg,vcomp_seasonal_avg,vcomp_month_avg,vcomp_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'vcomp','m/s')


winds_at_heightlevel(ucomp_avg,vcomp_avg,39,precipitation_avg,'fromwhite','mm/d',0.,8.,landmaskxr)

#[omega,omega_avg,omega_seasonal_avg,omega_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'omega','Pa/s',level=39)

winds_seasons(ucomp_seasonal_avg,vcomp_seasonal_avg,39,precipitation_seasonal_avg,'fromwhite','mm/d',0.,8.,landmaskxr,outdir,runmin,runmax)

winds_seasons(ucomp_seasonal_avg,vcomp_seasonal_avg,39,(precipitation_seasonal_avg-net_lhe_seasonal_avg),'rainnorm','mm/d',-2.,2.,landmaskxr,outdir,runmin,runmax)


# winds_anomaly_uv_vectors(ucomp_avg,vcomp_avg,landmaskxr, level = 37)

# # for i in range(0,10):
# #      winds_anomaly_uv_vectors(ucomp_annual_avg[i,:,:,:],vcomp_annual_avg[i,:,:,:],landmaskxr)


# winds_anomaly(ucomp_avg,vcomp_avg,landmaskxr)

# winds_anomaly(ucomp_seasonal_avg.sel(season='JJA'),vcomp_seasonal_avg.sel(season='JJA'),landmaskxr)
# winds_anomaly(ucomp_seasonal_avg.sel(season='DJF'),vcomp_seasonal_avg.sel(season='DJF'),landmaskxr)



# animated_map(testdir,slp_month_avg,'hPa','slp','slp_clim_animated','slp',0,12)

#any_configuration_plot(outdir,runmin,runmax,-90.,90.,flux_oceanq_avg,area_array,'mm/d','annual mean qflux','tempdiff',landmaskxr,contourson=False)

# print(area_integral((flux_oceanq_avg*0 + 1),area_array,landmaskxr,'all_sfcs')) # to test that integration
# # function works --> gives same result as total area 


# input_qflux_file = Dataset(os.path.join(GFDL_BASE,'input/all_continents/ocean_qflux.nc'))
# input_qflux = xr.DataArray(input_qflux_file.variables['ocean_qflux'][:], coords = [input_qflux_file.variables['time'][:],input_qflux_file.variables['lat'][:],input_qflux_file.variables['lon'][:]], dims = ['time','lat','lon'])
# input_qflux_avg = input_qflux.mean(dim = 'time')
# input_qflux_int = area_integral(input_qflux_avg,area_array,landmaskxr,'ocean')
# print('Global integrated q flux area weighted = '+str(input_qflux_int)) # should be almost identical
# integr_qflux = area_integral(flux_oceanq_avg,area_array,landmaskxr,'ocean')
# print('Global integrated q flux area weighted = '+str(integr_qflux))

# [convection_rain,convection_rain_avg,convection_rain_seasonal_avg,convection_rain_month_avg,time]=seasonal_surface_variable(testdir,runmin,runmax,'convection_rain','kg/m2s')
# [condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,condensation_rain_month_avg,time]=seasonal_surface_variable(testdir,runmin,runmax,'condensation_rain','kg/m2s')

# precipitation_avg = (condensation_rain + convection_rain).mean('time')*86400


# compare precip against gpcp data 

nc_gpcp = Dataset('/scratch/mp586/gpcp_detrendedpentads.nc')

gpcp_P = nc_gpcp.variables['precip_clim'][:]
gpcp_P = xr.DataArray(gpcp_P, coords = [nc_gpcp.variables['xofyear'][:],nc_gpcp.variables['lat'][:],nc_gpcp.variables['lon'][:]], dims = ['time','lat','lon'])

gpcp_P_avg = gpcp_P.mean('time')

# # interpolate gpcp data to model grid # not working 

# arlat = precipitation_avg.lat
# arlon = precipitation_avg.lon

# m = Basemap(projection='kav7',lon_0=0.,resolution='c')
# x_out, y_out = m(arlat,arlon)

# gpcp_precip = interp(gpcp_P,gpcp_P.lat,gpcp_P.lon,x_out,y_out) 


# worldmap_variable(gpcp_P_avg,'mm/day','GPCP precip annual mean','fromwhite',0.,8.,contourson=True)

gpcp_avg = area_weighted_avg(gpcp_P_avg,area_array_1deg,landmaskxr,'all_sfcs',minlat=-90., maxlat=90.,axis=None)
print('GPCP mean precip = '+str(gpcp_avg))

# if runmin == 1: 
#     for i in range(1,runmax):
#         month_plot = any_configuration_plot(outdir,runmin,runmax,-90.,90.,sphum[i,:,:],area_array,'kg/kg * dp','IWV','fromwhite',landmaskxr,contourson=False,month_annotate=i)
#         month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/IWV_month'+str(1000+i)+'.png',bbox_inches='tight')
#         os.system('convert -delay 50 /scratch/mp586/Code/Graphics/'+testdir+'/IWV_month*.png /scratch/mp586/Code/Graphics/'+testdir+'/IWV.gif')



# any_configuration_plot(outdir,runmin,runmax,-90.,90.,omega_avg,area_array,'hPa/s','omega','rainnorm',landmaskxr)


# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(bucket_depth_avg),area_array,'m','bucket_depth','fromwhite',landmaskxr)

# animated_map(testdir,flux_oceanq_month_avg,area_array,'W/m^2','resulting_q_flux','qflux_clim_animated','rainnorm',0,12,-300,300)


maxval_omega_surf = np.absolute(omega_month_avg.max())
minval_omega_surf = np.absolute(omega_month_avg.min())


# for i in range(1,13):

#     month_plot = winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,(omega_month_avg).sel(month=i),'rainnorm','dp/dt',minval_omega_surf,maxval_omega_surf)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_omega_run_'+str(runmin)+'-'+str(runmax)+'month_'+str(i)+'.png',bbox_inches='tight')
# os.system('convert -delay 100 /scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_omega_run_'+str(runmin)+'-'+str(runmax)+'*.png /scratch/mp586/Code/Graphics/'+testdir+'/omega_wind_monthly_clim_run_'+str(runmin)+'-'+str(runmax)+'.gif')


maxval_precip = np.absolute((precipitation_month_avg).max())
# maxval_sphum = np.absolute((sphum_avg).max())

for i in range(1,13):

    month_plot = winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,precipitation_month_avg.sel(month=i),'fromwhite','mm/day',0,8,landmaskxr)
    month_plot.savefig('/scratch/mp586/Code/Graphics/anim_plot_precip_run_'+str(runmin)+'-'+str(runmax)+'month_'+str(i)+'.png',bbox_inches='tight')
os.system('convert -delay 100 /scratch/mp586/Code/Graphics/anim_plot_precip_run_'+str(runmin)+'-'+str(runmax)+'*.png /scratch/mp586/Code/Graphics/precip_wind_monthly_clim_run_'+str(runmin)+'-'+str(runmax)+'.gif')


# for i in range(1,13):

#     month_plot = winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,tsurf_month_avg.sel(month=i),'temp','mm/day',240,310,landmaskxr)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/anim_plot_tsurf_run_'+str(runmin)+'-'+str(runmax)+'month_'+str(i)+'.png',bbox_inches='tight')
# os.system('convert -delay 100 /scratch/mp586/Code/Graphics/anim_plot_tsurf_run_'+str(runmin)+'-'+str(runmax)+'*.png /scratch/mp586/Code/Graphics/tsurf_wind_monthly_clim_run_'+str(runmin)+'-'+str(runmax)+'.gif')



PE_avg_sum = area_integral(PE_avg,area_array,landmaskxr,'all_sfcs',factor = 10**(-3)) # factor to convert from mm/d to m/d
print('P avg - E avg global integral / total sfc area'+str(PE_avg_sum/total_sfc_area))

###just a test -- gives same result
# PE_avg_sum = area_integral(precipitation_avg,area_array,landmaskxr,'all_sfcs',factor = 10**(-3)) - area_integral(net_lhe_avg,area_array,landmaskxr,'all_sfcs',factor = 10**(-3)) # factor to 
# #convert from mm/d to m/d
# print('P avg - E avg global integral / total sfc area'+str(PE_avg_sum/total_sfc_area))

#PE_avg_sum = area_integral(PE_avg,area_array,landmaskxr,'land',factor = 10**(-3)) # factor to 
#convert from mm/d to m/d
#print('P avg - E avg global integral land / landsfc area'+str(PE_avg_sum/land_sfc_area))

# Need runmin = 1 for all of those 
# bd_0 = area_integral(bucket_depth[0,:,:],area_array,landmaskxr,'all_sfcs')
# print('bd_0 allsfcs= '+str(bd_0/total_sfc_area))
# bd_end = area_integral(bucket_depth[runmax-2,:,:],area_array,landmaskxr,'all_sfcs')
# print('bd_end allsfcs= '+str(bd_end/total_sfc_area))

# bd_0 = area_integral(bucket_depth[0,:,:],area_array,landmaskxr,'land')
# print('bd_0 land= '+str(bd_0/land_sfc_area))
# bd_end = area_integral(bucket_depth[runmax-2,:,:],area_array,landmaskxr,'land')
# print('bd_end land= '+str(bd_end/land_sfc_area))

#any_configuration_plot(outdir,runmin,runmax,-90.,90.,bucket_depth_avg.where(landmask==1.),'m','bucket_depth','fromwhite',landmaskxr,minval=0.,maxval=2.)

# any_configuration_plot(outdir,runmin,runmax,-100.,100.,(bucket_depth[runmax-2,:,:] - bucket_depth[11,:,:]),area_array,'m','bucket depth dec year 40 - year 1','rainnorm',landmaskxr)

# any_configuration_plot(outdir,runmin,runmax,-100.,100.,(sphum[runmax-2,:,:] - sphum[11,:,:]),area_array,'m','IWV dec year 40 - year 1','rainnorm',landmaskxr)

#any_configuration_plot(outdir,runmin,runmax,-90.,90.,net_lhe_avg,area_array,'mm/day','E avg','fromwhite',landmaskxr)

# # any_configuration_plot(outdir,runmin,runmax,-90.,90.,convection_rain_avg,'mm/day','convection_rain avg','fromwhite')
# # any_configuration_plot(outdir,runmin,runmax,-90.,90.,condensation_rain_avg,'mm/day','condensation_rain avg','fromwhite')
#any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_avg.where(landmaskxr==1.),'mm/day','P avg','fromwhite')

#any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_month_avg.sel(month=7)-net_lhe_month_avg.sel(month=7),'mm/day','P-E_July (mm/day)','rainnorm',landmaskxr)
#any_configuration_plot(outdir,runmin,runmax,-90.,90.,precipitation_month_avg.sel(month=1)-net_lhe_month_avg.sel(month=1),'mm/day','P-E_January (mm/day)','rainnorm',landmaskxr)





####################### For comparing my run against Ruth's run ################

# # # ruth climatology for new bucket model formulation with E = E_0 when bucket is full 
# nc_test = Dataset('/scratch/mp586/GFDL_DATA/bucket_amendment_test.nc') 
# bucket_amended = xr.DataArray(nc_test.variables['bucket_depth'][:],coords=[precipitation_month_avg.month,precipitation.lat,precipitation.lon],dims=['month','lat','lon'])
# precip_amended = xr.DataArray(nc_test.variables['precipitation'][:],coords=[precipitation_month_avg.month,precipitation.lat,precipitation.lon],dims=['month','lat','lon'])*86400
# precip_amended_annavg = precip_amended.mean(dim='month')
# bucket_amended_annavg = bucket_amended.mean(dim='month')

# precip_diff_month_avg = precip_amended - precipitation_month_avg
# # bucket_diff_month_avg = bucket_amended - bucket_depth_month_avg

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(precip_amended_annavg - precipitation_avg),area_array,'mm/d','precip_rg_newbucket - precip this run avg','rainnorm',landmaskxr,contourson=False)
# # any_configuration_plot(outdir,runmin,runmax,-90.,90.,(bucket_amended_annavg - bucket_depth_avg),area_array,'m','bucket_new - bucket_old annual avg','rainnorm',landmaskxr,contourson=False)

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(precip_amended_annavg - precipitation_avg).where(landmask == 1.),area_array,'mm/d','precip_rg_newbucket - precip this run avg','rainnorm',landmaskxr,contourson=False)
# # any_configuration_plot(outdir,runmin,runmax,-90.,90.,(bucket_amended_annavg - bucket_depth_avg).where(landmask == 1.),area_array,'m','bucket_new - bucket_old annual avg','rainnorm',landmaskxr,contourson=False)
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precip_amended_annavg,area_array,'mm/d','precip_rg_newbucket','fromwhite',landmaskxr,nmb_contours=4,minval = 0., maxval = 8.)

# # for i in range(1,13):
# #     month_plot = any_configuration_plot(outdir,runmin,runmax,-90.,90.,precip_diff_month_avg.sel(month=i),area_array,'mm/d','precip_new - precip_old','rainnorm',landmaskxr,contourson=False,month_annotate=i)
# #     month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/precip_diff_clim_month'+str(i)+'.png',bbox_inches='tight')
# # os.system('convert -delay 100 /scratch/mp586/Code/Graphics/'+testdir+'/precip_diff_clim_month*.png /scratch/mp586/Code/Graphics/'+testdir+'/precip_diff_clim.gif')



# # ruth climatology for full qflux lhpref = 0.7 run lhpref outside and 5-daily data
# nc_ruth = Dataset('/scratch/mp586/GFDL_DATA/full_qflux.nc') 
# precip_ruth = xr.DataArray(nc_ruth.variables['convection_rain'][:] + nc_ruth.variables['condensation_rain'][:], coords = [nc_ruth.variables['xofyear'][:],precipitation_avg.lat,precipitation_avg.lon],dims = ['time','lat','lon'])*86400
# precip_ruth_avg = precip_ruth.mean(dim='time')

# precip_diff_annavg = precipitation_avg - precip_ruth_avg
# #bucket_diff_month_avg = bucket_amended - bucket_depth_month_avg

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precip_diff_annavg.where(landmask == 1.),area_array,'mm/d','precip_lhpref07_me - precip_ruth ann avg','rainnorm',landmaskxr,contourson=False)
# #any_configuration_plot(outdir,runmin,runmax,-90.,90.,(bucket_amended_annavg - bucket_depth_avg).where(landmask == 1.),area_array,'m','bucket_new - bucket_old annual avg','rainnorm',landmaskxr,contourson=False)


# precip_newbucket_minus_ruthlhpref = precip_amended_annavg - precip_ruth_avg
# #bucket_diff_month_avg = bucket_amended - bucket_depth_month_avg

# any_configuration_plot(outdir,runmin,runmax,-90.,90.,precip_newbucket_minus_ruthlhpref,area_array,'mm/d','precip_lhprefinside_07_ruth - precip_newbucket ann avg','rainnorm',landmaskxr,contourson=False)
# #any_configuration_plot(outdir,runmin,runmax,-90.,90.,(bucket_amended_annavg - bucket_depth_avg).where(landmask == 1.),area_array,'m','bucket_new - bucket_old annual avg','rainnorm',landmaskxr,contourson=False)



#############################################################################



# still need to change to area_weighted zonal averages
# tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','bucket_depth',1.,'k',0,1,runmax,landmaskxr)
# globavg_var_timeseries_total_and_land(testdir,'precipitation',1,runmax,86400,landmask)
# globavg_var_timeseries_total_and_land(testdir,'flux_lhe',1,runmax,1./28.,landmask)
# tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','rh',1.,'m',39,1,runmax,landmaskxr)
# tropics_severalvars_timeseries_oceanonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','rh',1.,'m',39,1,runmax,landmaskxr)
# tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','bucket_depth',1.,'k',0,1,runmax,landmaskxr)
# tropics_severalvars_timeseries_landonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','t_surf',1.,'r',0,1,runmax,landmaskxr)
# tropics_severalvars_timeseries_oceanonly(testdir,'precipitation',86400,'Blue','flux_lhe',1./28.,'g','t_surf',1.,'r',0,1,runmax,landmaskxr)




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

#any_configuration_plot_correlation(-90.,90.,tsurf_avg,precipitation_avg,'tsurf vs precip')




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

# several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=JJA),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=JJA),'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=JJA),'E (mm/day)','g','JJA T_S, P and E') # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=DJF),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=DJF),'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=DJF),'E (mm/day)','g','DJF T_S, P and E')
# several_vars_zonalavg2(tsurf_avg,'T_S (K)','Red',precipitation_avg,'P (mm/day)','Blue',net_lhe_avg,'E (mm/day)','g','avg T_S, P and E')
