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

testdir=input('Enter data directory name as string ')
runmin=97 #input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=241 #input('Enter runmax number ')


[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir,runmin,runmax,'t_surf','K')

[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir,runmin,runmax,'precipitation','kg/m2s')

PE_avg=precipitation_avg*86400-net_lhe_avg/28. # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC

plotting_routines_kav7.aquaplanet_plot_minuszonavg(-90.,90.,PE_avg,'mm/day','P-E avg minus zonavg','rainnorm')
plotting_routines_kav7.aquaplanet_plot_minuszonavg(-90.,90.,precipitation_avg*86400,'mm/day','P avg minus zonavg','raindefault')
plotting_routines_kav7.aquaplanet_plot(-90.,90.,PE_avg,'mm/day','P-E avg','rainnorm')
plotting_routines_kav7.aquaplanet_plot(-90.,90.,precipitation_avg*86400,'mm/day','P avg','rainnorm')
plotting_routines_kav7.aquaplanet_plot(-90.,90.,tsurf_avg,'K','tsurf avg','temp')
plotting_routines_kav7.aquaplanet_plot_minuszonavg(-90.,90.,tsurf_avg,'K','tsurf avg minus zonavg','temp')

plotting_routines.globavg_tsurf_timeseries(testdir,1,runmax)
