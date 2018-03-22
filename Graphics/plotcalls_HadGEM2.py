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


nc = Dataset('/scratch/mp586/HadGEM2_picontrol/hurs_Amon_HadGEM2-ES_piControl_r1i1p1_240512-243011.nc')
rh_sfc = nc.variables['hurs'][:]
time = nc.variables['time'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
rh_sfc = xr.DataArray(rh_sfc,coords=[time,lat,lon],dims=['time','lat','lon'])
rh_avg = rh_sfc.mean('time')
worldmap_variable(rh_avg,'%','avg rh_sfc','fromwhite',0,100)

nc = Dataset('/scratch/mp586/HadGEM2_picontrol/pr_Amon_HadGEM2-ES_piControl_r1i1p1_240512-243011.nc')
P = nc.variables['pr'][:]
P = xr.DataArray(P*86400,coords=[time,lat,lon],dims=['time','lat','lon'])
P_avg = P.mean('time')
worldmap_variable(P_avg,'mm/d','avg_P','fromwhite',0,8)

nc = Dataset('/scratch/mp586/HadGEM2_picontrol/hfls_Amon_HadGEM2-ES_piControl_r1i1p1_240512-243011.nc')
E = nc.variables['hfls'][:]
E = xr.DataArray(E/28.,coords=[time,lat,lon],dims=['time','lat','lon'])
E_avg = E.mean('time')
worldmap_variable(E_avg,'mm/d','avg_E','fromwhite',0,8)


rh_avg_tropics = rh_avg.sel(lat=slice(-30.,30.))
precip_avg_tropics = P_avg.sel(lat=slice(-30.,30.))
evap_avg_tropics = E_avg.sel(lat=slice(-30.,30.))

rh_tropics_1d = np.asarray(rh_avg_tropics).flatten()
P_tropics_1d = np.asarray(precip_avg_tropics).flatten()
E_tropics_1d = np.asarray(evap_avg_tropics).flatten()

PE_tropics_1d = P_tropics_1d - E_tropics_1d

plt.plot(rh_tropics_1d,PE_tropics_1d,'k*',label = 'P-E tropics (all sfcs)')
plt.legend()
plt.xlabel('RH %')
plt.ylabel('P-E (mm/d)')
plt.title('P-E versus RH, annual mean (tropics)')
plt.show()
