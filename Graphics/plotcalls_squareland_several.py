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
import stats as st

testdir1=input('Enter data directory number one as string ')
testdir2=input('Enter data directory number two as string ')
testdir3=input('Enter data directory number three as string ')

runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number ')

landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
landmask=landfile.variables['land_mask'][:]
lats=landfile.variables['lat'][:]
lons=landfile.variables['lon'][:]
 

[tsurf1,tsurf1_avg,tsurf1_seasonal_avg,tsurf1_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir1,runmin,runmax,'t_surf','K')
[net_lhe1,net_lhe1_avg,net_lhe1_seasonal_avg,net_lhe1_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir1,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[precip1,precip1_avg,precip1_seasonal_avg,precip1_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir1,runmin,runmax,'precipitation','kg/m2s')

PE1_avg=precip1_avg*86400-net_lhe1_avg/28. # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC


[tsurf2,tsurf2_avg,tsurf2_seasonal_avg,tsurf2_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir2,runmin,runmax,'t_surf','K')
[net_lhe2,net_lhe2_avg,net_lhe2_seasonal_avg,net_lhe2_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir2,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[precip2,precip2_avg,precip2_seasonal_avg,precip2_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir2,runmin,runmax,'precipitation','kg/m2s')

PE2_avg=precip2_avg*86400-net_lhe2_avg/28. 


[tsurf3,tsurf3_avg,tsurf3_seasonal_avg,tsurf3_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir3,runmin,runmax,'t_surf','K')
[net_lhe3,net_lhe3_avg,net_lhe3_seasonal_avg,net_lhe3_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir3,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[precip3,precip3_avg,precip3_seasonal_avg,precip3_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir3,runmin,runmax,'precipitation','kg/m2s')

PE3_avg=precip3_avg*86400-net_lhe3_avg/28.


plotting_routines.squareland_plot_several(-90.,90.,PE1_avg,'P-E avg 1',PE2_avg,'P-E avg 2',PE3_avg,'P-E avg 3','mm/day','P-E avg','rainnorm')

