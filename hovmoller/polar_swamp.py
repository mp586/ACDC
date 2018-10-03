#! /bin/python

# These are plots to show that the polar regions are a swamp
# 
# 1. Show a Hovmoller plot of humidity: 
#   -> zonal mean specific humidity seasonal cycle
#
# 2. Show a line plot of zonal mean bucket depth and the P-E there
#
# -----------------------------------------------------------------

import os
import sys
import xarray as xr
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import glob
import dask.array as da
from copy import deepcopy

# My functions
sys.path.insert(0,'../tools')
import ds_utils as dsu
import ds_plots as dsp

# --- Setup directories for saving
fig_dir = './figures'
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

# --- Load the ctrl dataset
data_dir = '/Users/tim/work/results/acdc/half_ocean_newbucket_0qflux_control'
globname = data_dir + '/*/*.nc'
file_list = glob.glob(globname)
file_list = sorted(file_list)
ds_ctrl = dsu.load_dataset(file_list)

# --- Drop first year of spinup
ds_ctrl = ds_ctrl.where(ds_ctrl['time'].dt.year.isin(range(2,11)),drop=True)

# --- Compute climatology
ds_clim = ds_ctrl.groupby('time.month').mean('time')

# --- Time vector and lat/lon
t = range(len(ds_ctrl['time']))
months = ds_clim['month']
lat = ds_ctrl['lat']
lon = ds_ctrl['lon']

# --- Specific humidity
plt.figure(figsize=(9,5))
plt.pcolormesh(months,lat,
              ds_clim['sphum'].mean('lon').isel(pfull=-1).T,
              cmap='Purples')
cbar=plt.colorbar()
cbar.set_label('kg/kg')
plt.title('Zonal mean surface specific humidity')
#plt.show()
fig_name = '%s/hov_specificHumidity' % fig_dir
plt.savefig(fig_name,transparent=True,pad_inches=0)


# --- Zonal line of pme and bucket depth
pme = ds_ctrl['precipitation'] - ds_ctrl['flux_lhe']
# Integrate pme to be in meters
pme.values = 30 * .001 * pme.values
dsp.zonalMeanVsTime(t,(ds_ctrl['bucket_depth'],pme),60,'mm/day')
plt.grid()
fig_name = '%s/polar_precip' % fig_dir
plt.savefig(fig_name,transparent=True,pad_inches=0)
#plt.show()
