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

outdir = 'HadGEM2_picontrol'

nc = Dataset('/scratch/mp586/HadGEM2_picontrol/hurs_Amon_HadGEM2-ES_piControl_r1i1p1_240512-243011.nc')
rh_sfc = nc.variables['hurs'][:]
time = nc.variables['time'][:]
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
rh_sfc = xr.DataArray(rh_sfc,coords=[time,lat,lon],dims=['time','lat','lon'])
rh_avg = rh_sfc.mean('time')
worldmap_variable(outdir,rh_avg,'%','avg_rh_sfc','fromwhite',0,100)

nc = Dataset('/scratch/mp586/HadGEM2_picontrol/pr_Amon_HadGEM2-ES_piControl_r1i1p1_240512-243011.nc')
P = nc.variables['pr'][:]
P = xr.DataArray(P*86400,coords=[time,lat,lon],dims=['time','lat','lon'])
P_avg = P.mean('time')
worldmap_variable(outdir,P_avg,'mm/d','avg_P','fromwhite',0,8,nmb_contours=6)

nc = Dataset('/scratch/mp586/HadGEM2_picontrol/hfls_Amon_HadGEM2-ES_piControl_r1i1p1_240512-243011.nc')
E = nc.variables['hfls'][:]
E = xr.DataArray(E/28.,coords=[time,lat,lon],dims=['time','lat','lon'])
E_avg = E.mean('time')
worldmap_variable(outdir,E_avg,'mm/d','avg_E','fromwhite',0,8)
PE_avg = P_avg - E_avg
worldmap_variable(outdir,PE_avg,'mm/d','avg_P-E','rainnorm',-3.,3.)


rh_avg_tropics = rh_avg.sel(lat=slice(-30.,30.))
precip_avg_tropics = P_avg.sel(lat=slice(-30.,30.))
evap_avg_tropics = E_avg.sel(lat=slice(-30.,30.))

rh_tropics_1d = np.asarray(rh_avg_tropics).flatten()
P_tropics_1d = np.asarray(precip_avg_tropics).flatten()
E_tropics_1d = np.asarray(evap_avg_tropics).flatten()

PE_tropics_1d = P_tropics_1d - E_tropics_1d

fig2, ax2 = plt.subplots(3,1,sharex = True,figsize = (25,10))
P_all_1d = np.asarray(P_avg.sel(lat=slice(-30.,30.))).flatten()
E_all_1d = np.asarray(E_avg.sel(lat=slice(-30.,30.))).flatten()
ax2[0].plot(rh_tropics_1d,P_all_1d,'b.', label='P tropics')
ax2[1].plot(rh_tropics_1d,E_all_1d,'g.', label='E tropics')
ax2[2].plot(rh_tropics_1d,PE_tropics_1d,'k.', label='P-E tropics')
ax2[0].legend()
ax2[1].legend()
ax2[2].legend()
ax2[2].set_xlabel('RH %')
ax2[0].set_ylabel('P (mm/d)')
ax2[1].set_ylabel('E (mm/d)')
ax2[2].set_ylabel('P-E (mm/d)')
fig2.suptitle('P, E and P-E versus RH (tropics)')
fig2.savefig('/scratch/mp586/Code/Graphics/HadGEM2_picontrol/RH_PE_land_oc_all_highres.png', bbox_inches='tight', dpi=400)

#all and tropics is equivalent, both means both land and ocean in the tropics (all sfcs)
