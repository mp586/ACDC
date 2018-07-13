# basically copy and pasted the necessary stuff from src/extra/python/scripts/calculate_qglux/land_qflux_zero_correction_zero_integral.py

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
model_data = 'Isca_DATA'
output_dir = 'Isca'

lmask = input('Which landmask? ')
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+lmask+'/land.nc'),mode='r')


landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]

# for specified lats

landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

area_array = ca.cell_area(t_res=42,base_dir='/scratch/mp586/Isca/')
area_array = xr.DataArray(area_array)

APqflux = Dataset(os.path.join(GFDL_BASE,'input/aquaplanet/isca_qflux/ocean_qflux.nc'))


lats = APqflux.variables['lat'][:]
lons = APqflux.variables['lon'][:]
time = APqflux.variables['time'][:]

APqflux = APqflux.variables['ocean_qflux'][:]

Cqflux = np.empty_like(APqflux)
# make corrected idealised continents qflux 
for i in range(0,12):
            qflux_i_init = APqflux[i,:,:]
            qflux_i_init[landmask == 1.] = 0. # 2.0
            qflux_i_init = xr.DataArray(qflux_i_init, coords = [lats,lons], dims = ['lat','lon'])
            Cqflux[i,:,:] = qflux_i_init

Cqflux = xr.DataArray(Cqflux,coords=[time,lats,lons],dims=['time','lat','lon'])
APqflux = xr.DataArray(APqflux,coords=[time,lats,lons],dims=['time','lat','lon'])

Cmean = Cqflux.mean('time')

# not the same whether I integrate over the mean or take the mean of the integrals (the latter is not supposed to be zero, because integral [month0] = integral [month6], integral [month1] = integral[month7] and so on... but the integral of the mean is zero because qflux[0] = qflux[6] mirrored, so they compensate in each grid cell 

Cmean_int = area_integral(Cmean,area_array,landmask,'all_sfcs')# doesn't matter whether I put option all sfcs or ocean, because qflux is zero over ocean anyway
num_oceancells = np.size(landmask) - np.count_nonzero(landmask)
q_int_invweight_landzero = (Cmean_int * 1./(area_array.where(landmask==0.)))/num_oceancells

q_int_invweight_landzero = xr.DataArray(q_int_invweight_landzero, coords = [lats,lons], dims = ['lat','lon'])

correction_matrix = np.expand_dims(q_int_invweight_landzero,axis=0)
correction_matrix = np.repeat(correction_matrix,12,axis=0) 
correction_matrix = xr.DataArray(correction_matrix,coords=[time,lats,lons],dims=['time','lat','lon'])
        # the correction matrix is the same in every month, since the annual mean of a 12xlatxlon matrix of ones is a latxlon matrix of ones, so similarly if the ones are replaces by q_int_invweight_landzero

print(area_integral(Cqflux.mean('time'),area_array,landmask,'all_sfcs') - area_integral(correction_matrix.mean('time'),area_array,landmask,'all_sfcs'))

qflux_out = Cqflux - correction_matrix

print(area_integral(qflux_out.mean('time'),area_array,landmask,'all_sfcs'))


qqflux_oout = np.empty_like(qflux_out)

for i in range(0,12):
    qqflux_oout[i,:,:] = qflux_out[i,:,:] - (qflux_out.mean('time') - (Cqflux.mean('time') - correction_matrix.mean('time')))
# in theory, the stuff in the brackets should be equal to zero, but it's roughly 10^-6 in each grid cell, so I am subtracting that from the already corrected qflux')

qqflux_oout = xr.DataArray(qqflux_oout,coords=[time,lats,lons],dims=['time','lat','lon'])

print(area_integral(qqflux_oout.mean('time'),area_array,landmask,'all_sfcs'))





# # ! What is going on here ??? 
# (Cqflux.mean('time') - (correction_matrix*0.).mean('time') - (Cqflux- correction_matrix*0.).mean('time')) # this should be a matrix of zeros, but it's a matrix of smthn*10^-6 s

# #whereas
# (Cqflux.mean('time') - (Cqflux*0.).mean('time') - (Cqflux- Cqflux*0.).mean('time'))

# # is a matrix of zeros...

# and Cqflux*0. == correction_matrix*0. is a matrix of trues.... 
####################################################################
# correction_matrix[0:5,:,:].mean('time') = correction_matrix[0,:,:]
# but 
# correction_matrix[0:11,:,:].mean('time') != correction_matrix[0,:,:] -- why????
###################################################################
# real continents qflux 
RCqflux = Dataset(os.path.join(GFDL_BASE,'input/all_continents/full_continents_newbucket/ocean_qflux.nc'))

RCqflux = xr.DataArray(RCqflux.variables['ocean_qflux'][:], coords = [time,lats,lons], dims = ['time','lat','lon'])

RCmask = Dataset(os.path.join(GFDL_BASE,'input/all_continents/land.nc'))

RCmask=landfile.variables['land_mask'][:]
