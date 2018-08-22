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

[qfluxcon,qfluxconmean,qfluxcon_seasonal_avg,qfluxcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_oceanq','W/m2')

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
tcono = tconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==0.)]
tpero = tpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==0.)]
prcono = prconmean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==0.)]
prpero = prpermean.sel(lat=slice(-30.,30.))[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==0.)]

# Finding m


# a and dq
# dTo = np.asarray(tpero - tcono).flatten()

# [mm,intercept,r,p,stderr] = linreg(dTl,drysurfenb) # aa = 8.4, dq = -32
# dq = -dq
# # plt.plot(dT_all,drysurfenb,'.b')
dTo_flat = np.asarray((tpero - tcono)).flatten()
qflux_flat = np.asarray(qfluxconmean[np.where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==0.)]).flatten()

plt.plot(qflux_flat,dTo_flat,'k.')

dPo_flat = np.asarray((prpero - prcono)).flatten()

plt.plot(dTo_flat,dPo_flat,'b.')
