#!/bin/bash/python

# Generate a precipitation plot 
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

def plotPrecip(ds,**kwargs):

    # Default value for months vector
    if 'month' in ds.coords.keys(): 
        months = ds['month']
    else:
        months = ds['time']

    for key in kwargs:
        if key == "months":
            months = kwargs[key]
        else:
            print('Unrecognized argument',key)

    # --- Open a nice looking figure
    plt.figure(figsize=(18,5))
    
    # --- P - E
    pme = ds['precipitation'].mean('lon').T - ds['flux_lhe'].mean('lon').T
    color_lim = np.max(np.abs(pme.values))
    plt.subplot(1,2,1)
    plt.pcolormesh(months,ds['lat'],
                   pme,
                   cmap='BrBG',vmin=-color_lim,vmax=color_lim)
    cbar = plt.colorbar()
    # Ideally, this would use ds['precipitation'].attrs['units']
    # but the groupby, compute climatology does not maintain attributes
    cbar.set_label('mm/day')
    plt.title('Zonal mean P-E')
    
    # --- Precipitation
    plt.subplot(1,2,2)
    if np.any(ds['precipitation'].values<0):
        pmax = np.max(np.abs(ds['precipitation'].mean('lon').T.values))
        plt.pcolormesh(months,ds['lat'],
                       ds['precipitation'].mean('lon').T,
                       cmap='BrBG',
                       vmin=-pmax,vmax=pmax)
    else:
        plt.pcolormesh(months,ds['lat'],
                       ds['precipitation'].mean('lon').T,
                       cmap='PuBu')

    cbar = plt.colorbar()
    cbar.set_label('mm/day')
    plt.title('Zonal mean precipitation')

    plt.show()

    return

def plotTempAndVel(ds,**kwargs):

    # Default value for months vector
    if 'month' in ds.coords.keys(): 
        months = ds['month']
    else:
        months = range(len(ds['time']))

    for key in kwargs:
        if key == "months":
            months = kwargs[key]
        else:
            print('Unrecognized argument',key)

    # --- Make a nice sized plot
    plt.figure(figsize=(18,5))

    # --- Surface temperature
    plt.subplot(1,2,1)
    plt.pcolormesh(months,ds['lat'],
                   ds['t_surf'].mean('lon').T,
                   cmap='RdBu_r',
                   vmin=-30,vmax=30)
    cbar = plt.colorbar()
    cbar.set_label('degC')
    plt.title('Zonal mean surface temperature')

    # --- Surface meridional velocity
    plt.subplot(1,2,2)
    plt.pcolormesh(months,ds['lat'],
                   ds['vcomp'].isel(pfull=-1).mean('lon').T,
                   cmap='RdBu_r',
                   vmin=-8,vmax=8)
    cbar = plt.colorbar()
    cbar.set_label('m/s')
    plt.title('Zonal mean meridional velocity')
    
    return

def plotSLP(ds,**kwargs):

    # Default value for months vector
    if 'month' in ds.coords.keys(): 
        months = ds['month']
    else:
        months = ds['time']

    for key in kwargs:
        if key == "months":
            months = kwargs[key]
        else:
            print('Unrecognized argument',key)

    # --- Nice sized plot
    plt.figure(figsize=(9,5))

    # --- Sea level pressure
    plt.pcolormesh(months,ds['lat'],
                   ds['slp'].mean('lon').T,
                   cmap='Blues')
    cbar = plt.colorbar()
    cbar.set_label('kPa')
    plt.title('Zonal mean SLP')
    
    return

def plotHumidity(ds,**kwargs):

    # Default value for months vector
    if 'month' in ds.coords.keys(): 
        months = ds['month']
    else:
        months = ds['time']

    for key in kwargs:
        if key == "months":
            months = kwargs[key]
        else:
            print('Unrecognized argument',key)

    # --- Make a nice sized plot
    plt.figure(figsize=(18,5))

    # --- Surface specific humidity
    plt.subplot(1,2,1)
    if np.any(ds['sphum']<0):
        sphum=ds['sphum'].mean('lon').isel(pfull=-1).T
        cmax = np.max(np.abs(sphum.values))
        plt.pcolormesh(months,ds['lat'],
                      sphum,
                      cmap='PRGn_r',
                      vmin=-cmax,vmax=cmax)
    else:
        plt.pcolormesh(months,ds['lat'],
                      ds['sphum'].mean('lon').isel(pfull=-1).T,
                      cmap='Purples')
    cbar=plt.colorbar()
    cbar.set_label('kg/kg')
    plt.title('Zonal mean surface specific humidity') 

    # --- Surface relative humidity
    plt.subplot(1,2,2)
    if np.any(ds['rh']<0):
        sphum=ds['rh'].mean('lon').isel(pfull=-1).T
        cmax = np.max(np.abs(sphum.values))
        plt.pcolormesh(months,ds['lat'],
                      sphum,
                      cmap='PRGn_r',
                      vmin=-cmax,vmax=cmax)
    else:
        plt.pcolormesh(months,ds['lat'],
                      ds['rh'].mean('lon').isel(pfull=-1).T,
                      cmap='Purples')
    cbar=plt.colorbar()
    cbar.set_label('%')
    plt.title('Zonal mean surface relative humidity') 

    return


