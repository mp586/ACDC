import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

import sys
sys.path.insert(0, '/scratch/mp586/git_repos/ACDC/tools')
from ds_plots import *
from netCDF4 import Dataset

# load dataset 1
nc = Dataset('/scratch/mp586/git_repos/ACDC/Data/monthlyclimo_'+input('Which experiment ? ')+'.nc')
lats = nc.variables['lat'][:]
lons = nc.variables['lon'][:]

P = xr.DataArray(nc.variables['precipitation'][:])*86400

mmap(np.asarray(P.mean('dim_0')),lats,lons,show = 'Yes')

