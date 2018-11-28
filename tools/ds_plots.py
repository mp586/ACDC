#!/bin/bash/python

# Generate a precipitation plot 
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point

def pmap(mapdata,lat,pressure,cmap=None,clim=None,title=None,
        newfig=True,filepath=None,show=None,reverse_ydir=True):
    # Make a map of a zonally averaged field as a function of pressure
    # Inputs
    #   mapdata     : zonally averaged field, function of latitude & pressure
    #   lat         : latitude vector
    #   pressure    : pressure vector
    #
    # Optional Inputs
    #   newfig      : if True, make a new figure
    #   reverse_ydir: if True, reverse ydir, good for pressure coordinates

    if newfig:
        fig = plt.figure(figsize=(8,7))
    else:
        fig = []

    # Make the plot
    plt.pcolormesh(lat,pressure,mapdata)

    # set colormap
    if cmap is not None:
        plt.set_cmap(cmap)
    else:
        plt.set_cmap(plt.cm.viridis)

    # nice colorbar
    cbar = plt.colorbar(orientation='vertical',extend='both',pad=.02, shrink=0.9)
    cbar.ax.tick_params(labelsize=14) 

    if clim is not None:
        plt.clim(clim)  

    if reverse_ydir:
        plt.gca().invert_yaxis()


    if title is not None:
        plt.title(title)

    if show is not None:
        plt.show()
    
    # save figure if path given
    if filepath is not None and newfig:
        fig.savefig(filepath+'.png',dpi=300,facecolor=None,edgecolor=None,bbox_inches='tight',transparent=True,pad_inches=0.01)
        
    
    return fig,cbar

def hmap(mapdata,lat,cmap=None,clim=None,title=None,
        newfig=True,filepath=None,show=None):
    # Make a hovmoller plot, assuming climatology (months)
    # of mapdata, which should already be zonally averaged
    # Inputs
    #   mapdata     : zonally averaged field, function of latitude & time
    #   lat         : latitude vector
    #
    # Optional Inputs

    if newfig:
        fig = plt.figure(figsize=(8,7))
    else:
        fig = []

    # Assumed time vector
    months = np.arange(12)

    # Make the plot
    if np.size(mapdata,0) == 12:
        plt.pcolormesh(months,lat,mapdata.T)
    else:
        plt.pcolormesh(months,lat,mapdata.T)

    # set colormap
    if cmap is not None:
        plt.set_cmap(cmap)
    else:
        plt.set_cmap(plt.cm.viridis)

    # nice colorbar
    cbar = plt.colorbar(orientation='vertical',extend='both',pad=.02, shrink=0.9)
    cbar.ax.tick_params(labelsize=14) 

    if clim is not None:
        plt.clim(clim)  

    if title is not None:
        plt.title(title)

    if show is not None:
        plt.show()
    
    # save figure if path given
    if filepath is not None:
        fig.savefig(filepath+'.png',dpi=300,facecolor=None,edgecolor=None,bbox_inches='tight',transparent=True,pad_inches=0.01)
        
    
    return fig,cbar

def mmap(mapdata,lat,lon,title=None,cb_ttl=None,cmap=None,clim=None,
              filepath=None,show=None,sigmask=None,p=None,
              do_zonal=None,newfig=True,
              outside_color = None, outside_val=None,
              use_pcolor=False):
    # Make a map (robinson projection) of mapdata over lat/lon
    # Inputs
    #   mapdata : 2D field to be plotted
    #   lat     : latitude coordinates   
    #   lon     : longitude coordinates   
    #
    # Optional inputs
    #   title   : map title
    #   cb_ttl  : colorbar title
    #   cmap    : special colormap
    #   clim    : colorbar limits, as a 2 element array
    #   use_pcolor: default is to use contourf. Set True for pcolormesh
    
    # mp added a few lines to convert potential xarray inputs to numpy arrays
    # so that add_cyclic_point works
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    mapdata = np.asarray(mapdata)

    if newfig:
        if do_zonal==1:
            # make 2 subplots, one for map, one for zonal average next to it
            fig = plt.figure(figsize=(8,7))
        else:
            fig = plt.figure(figsize=(8,7))
    else:
        fig = []
    
    ax = plt.axes(projection=ccrs.Robinson())
    #ax.coastlines()
    ax.set_global()

    # Remove ugly center line
    cyclic_data, cyclic_lons = add_cyclic_point(mapdata,coord=lon)

    # Make the plot
    # Default: contourf
    if not use_pcolor:
        cs = plt.contourf(cyclic_lons,lat,cyclic_data,transform=ccrs.PlateCarree())
    else:
        cs = plt.pcolormesh(cyclic_lons,lat,cyclic_data,transform=ccrs.PlateCarree())

    
    # set colormap
    if cmap:
        plt.set_cmap(cmap)
    else:
        plt.set_cmap(plt.cm.viridis)
    
    # nice looking title
    if title is not None:
        plt.title(title) #,{'fontsize' : 16}) #,y=1.05)
        
    # nice colorbar
    cbar = plt.colorbar(ax=ax,orientation='horizontal',extend='both',pad=.02, shrink=0.9)
    cbar.ax.tick_params(labelsize=14) 
    
    # set colorbar limits
    if clim:
        plt.clim(clim)  
        cs.set_clim(clim[0],clim[1])
        cs.set_clim(clim)

    # colorbar title
    if cb_ttl:
        cbar.set_label(cb_ttl,fontsize=14)
    
#     # If outside_val and outside_color, check if pos (set_over) or neg (set_under) to outside_color:
#     if outside_val:
#         cm = plt.get_cmap()
#         if outside_val=='left':
#             cm.set_under(color=outside_color)
#         elif outside_val=='right':
#             cm.set_over(color=outside_color)
#         else:
#             a=9# garbarge filler...
#     else:
#         cm=0.0
    cm = plt.get_cmap()
   
    # add hatching for significance
    if p:
        cyclic_sig, cyclic_lons = add_cyclic_point(sigmask,coord=lon)
        #hatch = plt.pcolor(cyclic_lons,lat,cyclic_sig,hatch='xxx',alpha=0.,transform=ccrs.PlateCarree()); 

        CLN, CLT = np.meshgrid(cyclic_lons,lat)

        # only put the hatches where the sigmask is >0.5
        #cyclic_sig = np.ma.masked_where(cyclic_sig<p,cyclic_sig)
        lat_sig = np.ma.masked_where(cyclic_sig<p,CLT)
        lat_sig = np.ma.masked_where(np.isnan(cyclic_sig),lat_sig)
        lon_sig = np.ma.masked_where(cyclic_sig<p,CLN)
        lon_sig = np.ma.masked_where(np.isnan(cyclic_sig),lon_sig)

        hatch = ax.scatter(lon_sig,lat_sig,marker='x',s=2,c=[0.6, 0.6, 0.6],alpha=0.8,transform=ccrs.PlateCarree(),)

    ax.patch.set_alpha(0.0)
    
    if show:
        plt.show()
    
    # save figure if path given
    if filepath:
        fig.savefig(filepath+'.png',dpi=300,facecolor=None,edgecolor=None,bbox_inches='tight',transparent=True,pad_inches=0.01)
        
    if not show:
        plt.close()

    return fig,ax,cs,cbar,cm

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
    fld_mean = np.zeros(len(fld)*2)
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
        fld_mean_lbl.append('Annual average')



    nearest_lat = fld[0]['lat'].sel(lat=sel_lat,method='nearest')
        
    latSign = 'N'
    if nearest_lat < 0:
        latSign = 'S'


    plt.title('Zonal mean at %0.1f%s' % (nearest_lat,latSign))
    plt.ylabel(units)
    plt.legend(fld_mean_lbl)

    return
