#!/bin/bash/python

# Generate a precipitation plot 
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

def zonalMeanVsTime(t,fld,sel_lat,units=''):
    # Plot the zonal mean of some field over time
    #
    # Inputs:
    #
    #   t = time vector
    #   fld = the variable, as an xarray
    #       or a tuple of fields to be plotted
    #   sel_lat = the (nearest) latitude band to grab
    #
    # Optional inputs:  
    #   
    #   units = string giving the units
    # 

    # --- Nice fig size for jupyter
    plt.figure(figsize=(9,5))

    # check if more than one field to plot
    if len(t) == len(fld):
        fld = (fld,)

    # make containers for annual mean
    fld_mean = np.zeros(len(fld))
    fld_mean_lbl = []

    # figure out number of years and make a time vector for all years
    Nyrs = int(np.ceil(len(t)/12))
    tyrs = range(Nyrs*12)

    for i in range(len(fld)):

        # Grab the data and annual averages
        line_to_plot = fld[i].mean('lon').sel(lat=sel_lat,method='nearest')

        # Repeat the annual averages over months
        annual_mean = fld[i].groupby('time.year').mean('time').mean('lon').sel(lat=sel_lat,method='nearest')
        annual_mean = np.tile(annual_mean.T.values,(12,1))
        annual_mean = annual_mean.T.reshape((1,12*Nyrs))

        # Plot line and annual average
        plt.plot(t,line_to_plot)
        plt.plot(tyrs,annual_mean.T)

        # Compute annual mean and add to the legend entries
        fld_mean[i] = line_to_plot.mean()
        fld_mean_lbl.append('Time mean = %0.1e' % fld_mean[i])



    nearest_lat = fld[0]['lat'].sel(lat=sel_lat,method='nearest')
        
    latSign = 'N'
    if nearest_lat < 0:
        latSign = 'S'


    plt.title('Zonal mean at %0.1f%s' % (nearest_lat,latSign))
    plt.ylabel(units)
    plt.legend(fld_mean_lbl)

    return
