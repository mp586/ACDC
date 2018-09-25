#!/bin/python

import xarray as xr
import cftime 


def load_dataset(filelist):
    # Load the dataset and fix the time coordinate
    # Input: list of file names *.nc to load

    with xr.set_options(enable_cftimeindex=True):
        ds = xr.open_mfdataset(filelist,autoclose='true',decode_times=False)

    time0 = ds['time']

    time2 = cftime.num2date(ds['time'][:], units=ds['time'].units,
                            calendar=ds['time'].calendar,
                            only_use_cftime_datetimes=True)

    ds['time'].values = time2

    ds = convert_units(ds)

    return ds

# ------------------------------------------------------------------------

def convert_units(ds):
    # Convert datset units to something a little more intuitive

    # Precipitation to mm/day
    ds['precipitation'].values = 86400 * ds['precipitation'].values
    ds['precipitation'].attrs['units']='mm/day'
    
    # Latent heat flux to mm/day of evaporation
    ds['flux_lhe'].values = 1/28 * ds['flux_lhe'].values
    ds['flux_lhe'].attrs['units']='mm/day'

    # Sea level pressure to kpa
    ds['slp'].values = 1 / 1000 * ds['slp'].values
    ds['slp'].attrs['units']='kPa'

    # Surface temperature to degC
    ds['t_surf'].values = -273 + ds['t_surf'].values
    ds['t_surf'].attrs['units']='degC'

    return ds


