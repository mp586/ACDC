# not sure about this, especially linear regression part not working (sometimes taking whole globe sometimes only tropics -- need to correct this script!

from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
from scipy.stats import linregress as linreg


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

landfile=Dataset(os.path.join(GFDL_BASE,'input/'+input('Which landmask ? ')+'/land.nc'),mode='r') # landmask both continents

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

landfileWC=Dataset(os.path.join(GFDL_BASE,'input/square_South_America/land.nc'),mode='r')

landmaskWC=landfileWC.variables['land_mask'][:] # land mask west continent
landmaskWCxr=xr.DataArray(landmaskWC,coords=[landlats,landlons],dims=['lat','lon']) 


landfileEC=Dataset(os.path.join(GFDL_BASE,'input/square_Africa/land.nc'),mode='r')

landmaskEC=landfileEC.variables['land_mask'][:] # land mask east continent
landmaskECxr=xr.DataArray(landmaskEC,coords=[landlats,landlons],dims=['lat','lon']) 


area_array = ca.cell_area(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/')
area_array = xr.DataArray(area_array)

sigma = 5.67*10**(-8)

[tper,tpermean,tper_seasonal_avg,tper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[lhper,lhpermean,lhper_seasonal_avg,lhper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','W/m^2',factor = 1.) # latent heat flux at surface (UP) in Wm-2
[prper,prpermean,prper_seasonal_avg,prper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','kg/m2s', factor=86400)
[swper,swpermean,swper_seasonal_avg,swper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_sw','W/m^2',factor = 1.) # 
[lwdper,lwdpermean,lwdper_seasonal_avg,lwdper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lw','W/m^2',factor = 1.) # 
[shper,shpermean,shper_seasonal_avg,shper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_t','W/m^2',factor = 1.) # 
[relper,relpermean,relper_seasonal_avg,relper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'rh','%',level=39)
[specper,specpermean,specper_seasonal_avg,specper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level=39)
lwupermean = sigma*(tpermean**4)
lwuper_month_avg = sigma*(tper_month_avg**4)

[tcon,tconmean,tcon_seasonal_avg,tcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[lhcon,lhconmean,lhcon_seasonal_avg,lhcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','W/m^2',factor = 1.) # latent heat flux at surface (UP) in Wm-2
[prcon,prconmean,prcon_seasonal_avg,prcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','kg/m2s', factor=86400)
[swcon,swconmean,swcon_seasonal_avg,swcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_sw','W/m^2',factor = 1.) # 
[lwdcon,lwdconmean,lwdcon_seasonal_avg,lwdcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lw','W/m^2',factor = 1.) # 
[shcon,shconmean,shcon_seasonal_avg,shcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_t','W/m^2',factor = 1.) # 
[relcon,relconmean,relcon_seasonal_avg,relcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%',level=39)
[speccon,specconmean,speccon_seasonal_avg,speccon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level=39)
lwuconmean = sigma*(tconmean**4)
lwucon_month_avg = sigma*(tcon_month_avg**4)


# annual means 

tconl = tconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
prconl = prconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
specconl = specconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
relconl = relconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
swconl = swconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
lwdconl = lwdconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
lwuconl = lwuconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
lhconl = lhconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
shconl = shconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]

tperl = tpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
prperl = prpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
specperl = specpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
relperl = relpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
swperl = swpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
lwdperl = lwdpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
lwuperl = lwupermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
lhperl = lhpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]
shperl = shpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]

relpero = relconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==0.)]
relcono = relconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==0.)]


# Finding the constants a, b, c, k, f, and dq


# a and dq
drysurfenb = np.asarray(((shperl+lwuperl-lwdperl-swperl) - (shconl+lwuconl-lwdconl-swconl))).flatten()
dTl = np.asarray(tperl - tconl).flatten()
[aa,dq,r,p,stderr] = linreg(dTl,drysurfenb) # aa = 8.4, dq = -32
dq = -dq
# plt.plot(dT_all,drysurfenb,'.b')


# b

prconmean_awavg = area_weighted_avg(prconmean,area_array,landmaskxr,'all_sfcs',-30.,30.)
prpermean_awavg = area_weighted_avg(prpermean,area_array,landmaskxr,'all_sfcs',-30.,30.)
tropprecipmean = (prpermean_awavg - prconmean_awavg) / prconmean_awavg  # mean relative precip change in the tropics

dTt_awavg = area_weighted_avg(tpermean - tconmean, area_array, landmaskxr,'all_sfcs',-30.,30.)

bb = tropprecipmean / dTt_awavg # %/K interest in tropical mean precipitation. bb = 0.021

# L*c

[lcc,intercept,r,p,stderr] = linreg(np.asarray(prconl).flatten(),np.asarray(lhconl).flatten()) #lcc = 20.97

# d

dd = 0.07 # assume C-C scaling

# f

curly_r = prconl / prconmean_awavg # Function that links local precipitation to tropical mean precipitation change. (Note we want to consider land data only for the individual points.) 
[ff,intercept,r,p,stderr] = linreg(np.asarray(relconl/100.).flatten(),np.asarray(curly_r).flatten()) #ff = 2.13

# k

kk = (area_weighted_avg(relpermean,area_array,landmaskxr,'ocean',-30.,30.) - area_weighted_avg(relconmean,area_array,landmaskxr,'ocean',-30.,30.)) / 100. # correction that comes from the fact that oceanic RH actually increases # kk = 0.021


dttrop = area_weighted_avg(tpermean - tconmean, area_array, landmaskxr, 'all_sfcs',-30.,30.) # tropical mean temperature change
dprtrop = prpermean_awavg - prconmean_awavg # tropical mean precipitation change 


dt = tpermean.sel(lat=slice(-30.,30.)) - tconmean.sel(lat=slice(-30.,30.)) # tropical temperature change
dpr = prpermean.sel(lat=slice(-30.,30.)) - prconmean.sel(lat=slice(-30.,30.)) # tropical precipitation change



# 1) Crude prediction: relies only on tropical mean P and T change 

dtpred = (dq - (((lcc * bb * prconmean.sel(lat=slice(-30.,30.)) + aa) * dttrop)) / aa # Note: these predictions will only make sense on land. Temperature change wrt tropical mean change.

dprpred = (bb * prconmean.sel(lat=slice(-30.,30.) * dttrop) - dprtrop # Precipitation change wrt tropical mean. (Chadwick prediction, effectively.)


any_configuration_plot(outdir,runmin,runmax,-100.,100.,dprpred,area_array,'mm/day','Pchange_crude_prediction_corrected','rainnorm',landmaskxr,minval=-.5,maxval=.5)

any_configuration_plot(outdir,runmin,runmax,-100.,100.,dtpred,area_array,'K','Tchange_crude_prediction_corrected','tempdiff',landmaskxr,minval=-1.,maxval=1.)

# 2) More advanced prediction that allows for shifts in precipitation caused by changes in land RH (that is used as a proxy for RH rank, curly_R).


dtpred2 = (dq - (aa * dttrop) - (lcc * ff * kk * prconmean_awavg) - (
        lcc * bb * prconmean.sel(lat=slice(-30.,30.) * dttrop)) / (
    aa - (lcc * ff * (relconmean.sel(lat=slice(-30.,30.)/100.)* prconmean_awavg * dd))
                                   
dprpred2 = (bb * prconmean.sel(lat=slice(-30.,30.) * dttrop) - (
    ff * (relconmean.sel(lat=slice(-30.,30.)/100.)* prconmean_awavg * dd * dtpred2) + (
    ff * kk * prconmean_awavg) - dprtrop

dprcheck2 = (bb * prconmean.sel(lat=slice(-30.,30.) * dttrop) - (
    ff * (relconmean.sel(lat=slice(-30.,30.)/100.)* prconmean_awavg * dd * (dt - dttrop)) + (
    ff * kk * prconmean_awavg) - dprtrop # If we give our precipitation prediction the right temperature change, how does it do? A bit better but not much better. However, the consistency on the delta p versus delta t plot is much better.

any_configuration_plot(outdir,runmin,runmax,-100.,100.,dprpred2,area_array,'mm/day','Pchange_advanced_prediction_corrected','rainnorm',landmaskxr,minval=-.5,maxval=.5)

any_configuration_plot(outdir,runmin,runmax,-100.,100.,dtpred2,area_array,'K','Tchange_advanced_prediction_corrected','tempdiff',landmaskxr,minval=-6.,maxval=6.)



# monthly means
# same constants as for annual mean 


area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 12, axis = 0) # to make area_array 3D (months, lat, lon)

# vectors of len 12 containing area weighted avg for each month
prcon_month_avg_awavg = area_weighted_avg_4D(prcon_month_avg,area_array_3D,landmaskxr,'all_sfcs',-30.,30.)
prper_month_avg_awavg = area_weighted_avg_4D(prper_month_avg,area_array_3D,landmaskxr,'all_sfcs',-30.,30.)
dttrop_months = area_weighted_avg_4D(tper_month_avg - tcon_month_avg, area_array_3D, landmaskxr, 'all_sfcs',-30.,30.) # tropical _month_avg temperature change

# arrays of dims 12xlatxlon
dprtrop_months = prper_month_avg_awavg - prcon_month_avg_awavg # tropical _month_avg[i,:,:] precipitation change 
dt_months = tper_month_avg.sel(lat=slice(-30.,30.) - tcon_month_avg.sel(lat=slice(-30.,30.) # tropical temperature change
dpr_months = prper_month_avg.sel(lat=slice(-30.,30.) - prcon_month_avg.sel(lat=slice(-30.,30.) # tropical precipitation change



# 1) Crude prediction: relies only on tropical _month_avg[i,:,:] P and T change 

dtpred_months = (dq - (((lcc * bb * prcon_month_avg.sel(lat=slice(-30.,30.) + aa) * dttrop_months))) / aa # Note: these predictions will only make sense on land. Temperature change wrt tropical _month_avg[i,:,:] change.

dprpred_months = (bb * prcon_month_avg.sel(lat=slice(-30.,30.) * dttrop_months) - dprtrop_months # Precipitation change wrt tropical _month_avg[i,:,:]. (Chadwick prediction, effectively.)


# 2) More advanced prediction that allows for shifts in precipitation caused by changes in land RH (that is used as a proxy for RH rank, curly_R).

dtpred2_months = (dq - (aa * dttrop_months) - (lcc * ff * kk * prcon_month_avg_awavg) - (lcc * bb * prcon_month_avg.sel(lat=slice(-30.,30.) * dttrop_months)) / (aa - (lcc * ff * (relcon_month_avg.sel(lat=slice(-30.,30.)/100.)* prcon_month_avg_awavg * dd))

dprpred2_months = (bb * prcon_month_avg.sel(lat=slice(-30.,30.) * dttrop_months) - (ff * (relcon_month_avg.sel(lat=slice(-30.,30.)/100.)* prcon_month_avg_awavg * dd * dtpred2_months) + (ff * kk * prcon_month_avg_awavg) - dprtrop_months

dprcheck2_months = (bb * prcon_month_avg.sel(lat=slice(-30.,30.) * dttrop_months) - (ff * (relcon_month_avg.sel(lat=slice(-30.,30.)/100.)* prcon_month_avg_awavg * dd * (dt_months.sel(lat=slice(-30.,30.) - dttrop_months)) + (ff * kk * prcon_month_avg_awavg) - dprtrop_months

# If we give our precipitation prediction the right temperature change, how does it do? A bit better but not much better. However, the consistency on the delta p versus delta t plot is much better.


for i in range(0,12): # plots for each month
    any_configuration_plot(outdir,runmin,runmax,-100.,100.,dprpred_months[i,:,:],area_array,'mm/day','Pchange_crude_prediction_month'+str(i+1),'rainnorm',landmaskxr,minval=-.5,maxval=.5)

    any_configuration_plot(outdir,runmin,runmax,-100.,100.,dtpred_months[i,:,:],area_array,'K','Tchange_crude_prediction_month'+str(i+1),'tempdiff',landmaskxr,minval=-2.,maxval=2.)

    any_configuration_plot(outdir,runmin,runmax,-100.,100.,dt_months[i,:,:],area_array,'K','Tchange_actual_month'+str(i+1),'tempdiff',landmaskxr,minval=-6.,maxval=6.)

    
    any_configuration_plot(outdir,runmin,runmax,-100.,100.,tcon_month_avg[i,:,:],area_array,'K','Tctl_actual_month'+str(i+1),'temp',landmaskxr)


for i in range(0,12): # plots for each month

    any_configuration_plot(outdir,runmin,runmax,-100.,100.,(dtpred2_months[i,:,:]).where(landmask==1.),area_array,'K','Tchange_advanced_prediction_month'+str(i+1),'tempdiff',landmaskxr,minval=-6.,maxval=6.)


    any_configuration_plot(outdir,runmin,runmax,-100.,100.,(dprpred2_months[i,:,:]).where(landmask==1.),area_array,'mm/day','Pchange_advanced_prediction_month'+str(i+1),'rainnorm',landmaskxr,minval=-1.,maxval=1.)




any_configuration_plot(outdir,runmin,runmax,-100.,100.,(dprpred2_months.mean('month')).where(landmask == 1.),area_array,'mm/day','Pchange_advanced_prediction_avg_allmonths','rainnorm',landmaskxr,minval=-1.,maxval=1.)

any_configuration_plot(outdir,runmin,runmax,-100.,100.,(dprcheck2_months.mean('month')).where(landmask == 1.),area_array,'mm/day','Pchange_check_avg_allmonths','rainnorm',landmaskxr,minval=-1.,maxval=1.)



any_configuration_plot(outdir,runmin,runmax,-100.,100.,(dprpred_months.mean('month')),area_array,'mm/day','Pchange_crude_prediction_avg_allmonths','rainnorm',landmaskxr,minval=-.5,maxval=.5)

any_configuration_plot(outdir,runmin,runmax,-100.,100.,(dtpred2_months.mean('month')),area_array,'K','Tchange_advanced_prediction_avg_allmonths','tempdiff',landmaskxr,minval=-6.,maxval=6.)


any_configuration_plot(outdir,runmin,runmax,-100.,100.,(dtpred_months.mean('month')),area_array,'K','Tchange_crude_prediction_avg_allmonths','tempdiff',landmaskxr,minval=-6.,maxval=6.)




#### x-y plots


##### temperature predictions

fig, ax = plt.subplots(1,2,sharey=True,figsize=(25,15))

ax[0].plot(np.asarray(dtpred.where(landmaskWC == 1.)).flatten(),np.asarray(dt.where(landmaskWC == 1.)-dttrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[0].plot(np.asarray(dtpred.where(landmaskEC == 1.)).flatten(),np.asarray(dt.where(landmaskEC == 1.)-dttrop).flatten(),'r.',label = 'East C') # East Island
ax[0].set_xlabel('dt_simple_pred')
ax[0].set_ylabel('dt_actual - dt_actual_tropmean')
ax[0].set_xlim(-1.0,2.5)
ax[0].legend()


ax[1].plot(np.asarray(dtpred2.where(landmaskWC == 1.)).flatten(),np.asarray(dt.where(landmaskWC == 1.)-dttrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[1].plot(np.asarray(dtpred2.where(landmaskEC == 1.)).flatten(),np.asarray(dt.where(landmaskEC == 1.)-dttrop).flatten(),'r.',label='East C') # East Island
ax[1].set_xlabel('dt_advanced_pred')
ax[1].legend()
ax[1].set_xlim(-1.0,2.5)

fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/dT_simple_vs_advanced_pred_'+str(runmin)+'-'+str(runmax), bbox_inches='tight')




fig, ax = plt.subplots(1,2,sharey=True,figsize=(25,15))

ax[0].plot(np.asarray((dtpred_months.mean('month')).where(landmaskWC == 1.)).flatten(),np.asarray(dt.where(landmaskWC == 1.)-dttrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[0].plot(np.asarray((dtpred_months.mean('month')).where(landmaskEC == 1.)).flatten(),np.asarray(dt.where(landmaskEC == 1.)-dttrop).flatten(),'r.',label = 'East C') # East Island
ax[0].set_xlabel('dt_simple_pred_months_mean')
ax[0].set_ylabel('dt_actual - dt_actual_tropmean')
ax[0].legend()
ax[0].set_xlim(-1.0,2.5)


ax[1].plot(np.asarray((dtpred2_months.mean('month')).where(landmaskWC == 1.)).flatten(),np.asarray(dt.where(landmaskWC == 1.)-dttrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[1].plot(np.asarray((dtpred2_months.mean('month')).where(landmaskEC == 1.)).flatten(),np.asarray(dt.where(landmaskEC == 1.)-dttrop).flatten(),'r.',label='East C') # East Island
ax[1].set_xlabel('dt_advanced_pred_months_mean')
ax[1].legend()
ax[1].set_xlim(-1.0,2.5)


##### precipitation predictions

fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/dT_simple_vs_advanced_pred_months_mean_'+str(runmin)+'-'+str(runmax), bbox_inches='tight')



fig, ax = plt.subplots(1,3,sharey=True,figsize=(40,15))

ax[0].plot(np.asarray(dprpred.where(landmaskWC == 1.)).flatten(),np.asarray(dpr.where(landmaskWC == 1.)-dprtrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[0].plot(np.asarray(dprpred.where(landmaskEC == 1.)).flatten(),np.asarray(dpr.where(landmaskEC == 1.)-dprtrop).flatten(),'r.',label = 'East C') # East Island
ax[0].set_xlabel('dpr_simple_pred')
ax[0].set_ylabel('dpr_actual - dpr_actual_tropmean')
ax[0].legend()
ax[0].set_xlim(-0.5,0.5)


ax[1].plot(np.asarray(dprpred2.where(landmaskWC == 1.)).flatten(),np.asarray(dpr.where(landmaskWC == 1.)-dprtrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[1].plot(np.asarray(dprpred2.where(landmaskEC == 1.)).flatten(),np.asarray(dpr.where(landmaskEC == 1.)-dprtrop).flatten(),'r.',label='East C') # East Island
ax[1].set_xlabel('dpr_advanced_pred')
ax[1].legend()
ax[1].set_xlim(-0.5,0.5)


ax[2].plot(np.asarray(dprcheck2.where(landmaskWC == 1.)).flatten(),np.asarray(dpr.where(landmaskWC == 1.)-dprtrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[2].plot(np.asarray(dprcheck2.where(landmaskEC == 1.)).flatten(),np.asarray(dpr.where(landmaskEC == 1.)-dprtrop).flatten(),'r.',label='East C') # East Island
ax[2].set_xlabel('dpr_check')
ax[2].legend()
ax[2].set_xlim(-0.5,0.5)

fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/dpr_simple_vs_advanced_pred_vs_check_'+str(runmin)+'-'+str(runmax), bbox_inches='tight')



fig, ax = plt.subplots(1,3,sharey=True,figsize=(40,15))

ax[0].plot(np.asarray((dprpred_months.mean('month')).where(landmaskWC == 1.)).flatten(),np.asarray(dpr.where(landmaskWC == 1.)-dprtrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[0].plot(np.asarray((dprpred_months.mean('month')).where(landmaskEC == 1.)).flatten(),np.asarray(dpr.where(landmaskEC == 1.)-dprtrop).flatten(),'r.',label = 'East C') # East Island
ax[0].set_xlabel('dpr_simple_pred_months_mean')
ax[0].set_ylabel('dpr_actual - dpr_actual_tropmean')
ax[0].legend()
ax[0].set_xlim(-0.5,0.5)


ax[1].plot(np.asarray((dprpred2_months.mean('month')).where(landmaskWC == 1.)).flatten(),np.asarray(dpr.where(landmaskWC == 1.)-dprtrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[1].plot(np.asarray((dprpred2_months.mean('month')).where(landmaskEC == 1.)).flatten(),np.asarray(dpr.where(landmaskEC == 1.)-dprtrop).flatten(),'r.',label='East C') # East Island
ax[1].set_xlabel('dpr_advanced_pred_months_mean')
ax[1].legend()
ax[1].set_xlim(-0.5,0.5)


ax[2].plot(np.asarray(dprcheck2_months.mean('month').where(landmaskWC == 1.)).flatten(),np.asarray(dpr.where(landmaskWC == 1.)-dprtrop).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[2].plot(np.asarray(dprcheck2_months.mean('month').where(landmaskEC == 1.)).flatten(),np.asarray(dpr.where(landmaskEC == 1.)-dprtrop).flatten(),'r.',label='East C') # East Island
ax[2].set_xlabel('dpr_check_months_mean')
ax[2].legend()
ax[2].set_xlim(-0.5,0.5)

fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/dpr_simple_vs_advanced_pred_months_mean_vs_check_'+str(runmin)+'-'+str(runmax), bbox_inches='tight')

fig, ax = plt.subplots(1,2,sharey=True,figsize=(25,15))

ax[0].plot(np.asarray(dt.where(landmaskWC == 1.)).flatten(),np.asarray(dpr.where(landmaskWC == 1.)).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[0].plot(np.asarray(dt.where(landmaskEC == 1.)).flatten(),np.asarray(dpr.where(landmaskEC == 1.)).flatten(),'r.',label = 'East C') # East Island
ax[0].set_xlabel('dt_actual')
ax[0].set_ylabel('dpr_actual')
#ax[0].set_xlim(-1.0,2.5)
ax[0].legend()

ax[1].plot(np.asarray(dt.where(landmask == 0.)).flatten(),np.asarray(dpr.where(landmask == 0.)).flatten(),'g.',label='ocean') # West Island. Temperature change.
ax[1].set_xlabel('dt_actual')
ax[1].legend()
#ax[1].set_xlim(-1.0,2.5)

fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/dP_vs_dT_each_continent_and_ocean_'+str(runmin)+'-'+str(runmax), bbox_inches='tight')

fig, ax = plt.subplots(1,2,figsize=(25,15))

ax[0].plot(np.asarray(dtpred2_months.mean('month').where(landmaskWC == 1.)).flatten(),np.asarray(dprpred2_months.mean('month').where(landmaskWC == 1.)).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[0].plot(np.asarray(dtpred2_months.mean('month').where(landmaskEC == 1.)).flatten(),np.asarray(dprpred2_months.mean('month').where(landmaskEC == 1.)).flatten(),'r.',label = 'East C') # East Island
ax[0].set_xlabel('dt_advanced_pred_month_mean')
ax[0].set_ylabel('dpr_advanced_pred_month_mean')
#ax[0].set_xlim(-1.0,2.5)
ax[0].legend()

ax[1].plot(np.asarray(dtpred2_months.mean('month').where(landmask == 0.)).flatten(),np.asarray(dprpred2_months.mean('month').where(landmask == 0.)).flatten(),'g.',label='ocean') # West Island. Temperature change.
ax[1].set_xlabel('dt_advanced_pred_month_mean')
ax[1].legend()
#ax[1].set_xlim(-1.0,2.5)

fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/dP_vs_dT_each_continent_and_ocean_advanced_month_mean_predictions_'+str(runmin)+'-'+str(runmax), bbox_inches='tight')

#############################################################################################################################
### LINEAR REGRESIONS ###

maskWC = ~np.isnan(np.asarray(dt[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]).flatten())
maskEC = ~np.isnan(np.asarray(dt[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)]).flatten())
maskOC = ~np.isnan(np.asarray(dt[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==0.)]).flatten())

fig, ax = plt.subplots(1,2,sharey=True,figsize=(25,15))

ax[0].plot(np.asarray(dt.where(landmaskWC == 1.)).flatten(),np.asarray(dpr.where(landmaskWC == 1.)).flatten(),'b.',label='West C') # West Island. Temperature change.
ax[0].plot(np.asarray(dt.where(landmaskEC == 1.)).flatten(),np.asarray(dpr.where(landmaskEC == 1.)).flatten(),'r.',label = 'East C') # East Island
ax[0].set_xlabel('dt_actual')
ax[0].set_ylabel('dpr_actual')
#ax[0].set_xlim(-1.0,2.5)
ax[0].legend()

[k,dy,r,p,stderr] = linreg(np.asarray(dt.where(landmaskWC == 1.)).flatten()[maskWC],np.asarray(dpr.where(landmaskWC == 1.)).flatten()[maskWC]) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(np.asarray(dt.where(landmaskWC == 1.)).flatten()[maskWC]),np.max(np.asarray(dt.where(landmaskWC == 1.)).flatten()[maskWC]),500)
y = k*x1 + dy
ax[0].plot(x1,y,'b-')
ax[0].annotate('West C: k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r_val = '+str("%.2f" % r), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)

[k,dy,r,p,stderr] = linreg(np.asarray(dt.where(landmaskWC == 1.)).flatten()[maskEC],np.asarray(dpr.where(landmaskWC == 1.)).flatten()[maskEC]) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(np.asarray(dt.where(landmaskWC == 1.)).flatten()[maskEC]),np.max(np.asarray(dt.where(landmaskWC == 1.)).flatten()[maskEC]),500)
y = k*x1 + dy
ax[0].plot(x1,y,'r-')
ax[0].annotate('East C: k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r_val = '+str("%.2f" % r), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)

ax[1].plot(np.asarray(dt.where(landmask == 0.)).flatten(),np.asarray(dpr.where(landmask == 0.)).flatten(),'g.',label='ocean') # West Island. Temperature change.
ax[1].set_xlabel('dt_actual')
ax[1].legend()
#ax[1].set_xlim(-1.0,2.5)


[k,dy,r,p,stderr] = linreg(np.asarray(dt.where(landmask == 0.)).flatten()[maskOC],np.asarray(dpr.where(landmask == 0.)).flatten()[maskOC]) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(np.asarray(dt.where(landmask == 0.)).flatten()[maskOC]),np.max(np.asarray(dt.where(landmask == 0.)).flatten()[maskOC]),500)
y = k*x1 + dy
ax[1].plot(x1,y,'g-')
ax[1].annotate('ocean: k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r_val = '+str("%.2f" % r), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)


fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/dP_vs_dT_each_continent_and_ocean_'+str(runmin)+'-'+str(runmax), bbox_inches='tight')

