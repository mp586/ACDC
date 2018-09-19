from netCDF4 import Dataset
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid, addcyclic
import xarray as xr
import pandas as pd
import os
import matplotlib.colors as colors
import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES') # personal module
import stats as st
from scipy import signal
GFDL_BASE = os.environ['GFDL_BASE']
from scipy import stats
from scipy.stats import linregress as linreg


class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
# This maps minimum, midpoint and maximum values to 0, 0.5 and 1 respectively i.e. forcing the middle of the scale bar to match your midpoint


# not using area weights yet
# def globavg_var_timeseries(testdir,varname,runmin,runmax):
# #	""" Plots & returns time series of specified variable """
# # runmin = first month for timeseries
# # runmax = last month for timeseries
# # testdir = experiment directorz name in GFDL_DATA

# # varname = variable name, e.g. 't_surf'

#     for i in range (runmin,runmax): # excludes last one!
#         runnr="{0:04}".format(i)
#         filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
#         nc = Dataset(filename,mode='r')
#         var=nc.variables[varname][:]
#         if len(np.shape(var))==4: # 4d array, sum over height
# 		var=(xr.DataArray(var)).sum(dim='dim_1')
    
#         if i==runmin:
#             timeseries=[var.mean()] # make timeseries be a list, not a float so that I can append later
#             print(type(timeseries))
# 	else:
#             timeseries.append(var.mean())

#     timeseries=np.asarray(timeseries)
#     months=np.linspace(runmin,(runmax-1),timeseries.size)
    
#     plt.plot(months,timeseries)
#     plt.title('globavg '+varname)
#     plt.xlabel('Month #')
#     plt.ylabel('Global average')
#     plt.grid(b=True)
#     plt.show()    
#     return(timeseries)


#NOT YET USING AREA WEIGHTED AVG
# def tropics_severalvars_timeseries_landonly(testdir,varname1,factor1,color1,varname2,factor2,color2,varname3,factor3,color3,height3,runmin,runmax,landmaskxr):
# #
# #	""" Plots timeseries of three variables over tropical land.
# 	# Using three different axes.  
# 	# The colors are chosen at input for each variable
# 	# Variable 3 can be at a specified height. If var3 is 3D, set height3 to 0 ""


#     for i in range (runmin,runmax):
# 		    runnr="{0:04}".format(i)
# 		    filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
# 		    nc = Dataset(filename,mode='r')

# 		    lats = nc.variables['lat'][:]
# 		    lons = nc.variables['lon'][:]
		    
# 		    var1=(xr.DataArray(nc.variables[varname1][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor1
# 		    var2=(xr.DataArray(nc.variables[varname2][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor2
# 		    if height3 != 0:
# 			    var3=(xr.DataArray(nc.variables[varname3][0,height3,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3
# 		    else:
# 			    var3=(xr.DataArray(nc.variables[varname3][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3

# 		    var1 = var1.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)
# 		    var2 = var2.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)
# 		    var3 = var3.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)


# 		    if i==runmin:
# 			    timeseries1=[var1.mean()] 
# 			    timeseries2=[var2.mean()]
# 			    timeseries3=[var3.mean()]
# 		    else:
# 			    timeseries1.append(var1.mean())
# 			    timeseries2.append(var2.mean())
# 			    timeseries3.append(var3.mean())

#     timeseries1=np.asarray(timeseries1)
#     timeseries2=np.asarray(timeseries2)
#     timeseries3=np.asarray(timeseries3)

#     months=np.linspace(runmin,(runmax-1),timeseries1.size)
    

#     fig, ax = plt.subplots()
# # Twin the x-axis twice to make independent y-axes.
#     axes = [ax, ax.twinx(), ax.twinx()] 
# # Make some space on the right side for the extra y-axis.
#     fig.subplots_adjust(right=0.75)
# ## Move the last y-axis spine over to the right by 20% of the width of the axes
#     axes[-1].spines['right'].set_position(('axes', 1.2))

# # To make the border of the right-most axis visible, we need to turn the frame
# # on. This hides the other plots, however, so we need to turn its fill off.
#     axes[-1].set_frame_on(True)
#     axes[-1].patch.set_visible(False)

#     colors = (color1, color2, color3)
#     field = (timeseries1,timeseries2,timeseries3)
#     names = (varname1,varname2,varname3)
#     i=0
#     for ax, color in zip(axes, colors):
#         data = field[i]
# 	if (i==2 and height3 !=0):
# 		ax.plot(months,data,color=color,label=names[i]+' at pfull = '+str(int(nc.variables['pfull'][height3]))+' hPa')
# 	else:	
# 		ax.plot(months,data,color=color,label=names[i])
#         ax.set_ylabel(names[i])
# 	if (i<=1):
# 		ax.set_ylim(0,timeseries1.max()+1.)
#         ax.tick_params(axis='y',colors=color)
#         i+=1
#     axes[0].set_xlabel('month #')
#     plt.legend()
#     plt.title('Tropical Land Only')
#     plt.show()


#NOT YET USING AREA WEIGHTED AVG
# def tropics_severalvars_timeseries_oceanonly(testdir,varname1,factor1,color1,varname2,factor2,color2,varname3,factor3,color3,height3,runmin,runmax,landmaskxr):

# 	# """ Plots timeseries of three variables over tropical ocean.
# 	# Using three different axes.  
# 	# The colors are chosen at input for each variable
# 	# Variable 3 can be at a specified height. If var3 is 3D, set height3 to 0 """

# # lat slice hard coded to -30. - 30. 
# # factor is needed to convert eg precip from kg/s to mm/day 

#     for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
# 	    runnr="{0:04}".format(i)
# 	    filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
# 	    nc = Dataset(filename,mode='r')
# 	    lats = nc.variables['lat'][:]
# 	    lons = nc.variables['lon'][:]
		    
# 	    var1=(xr.DataArray(nc.variables[varname1][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor1
# 	    var2=(xr.DataArray(nc.variables[varname2][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor2
# 	    if height3 != 0:
# 		    var3=(xr.DataArray(nc.variables[varname3][0,height3,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3
# 	    else:
# 		    var3=(xr.DataArray(nc.variables[varname3][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3
# 	    var1 = var1.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))!=1.)
# 	    var2 = var2.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))!=1.)
# 	    var3 = var3.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))!=1.)


# 	    if i==runmin:
# 		    timeseries1=[var1.mean()]
# 		    timeseries2=[var2.mean()]
# 		    timeseries3=[var3.mean()]
# 	    else:
# 		    timeseries1.append(var1.mean())
# 		    timeseries2.append(var2.mean())
# 		    timeseries3.append(var3.mean())

#     timeseries1=np.asarray(timeseries1)
#     timeseries2=np.asarray(timeseries2)
#     timeseries3=np.asarray(timeseries3)

#     months=np.linspace(runmin,(runmax-1),timeseries1.size)
    

#     fig, ax = plt.subplots()

# # Twin the x-axis twice to make independent y-axes.
#     axes = [ax, ax.twinx(), ax.twinx()] 

# # Make some space on the right side for the extra y-axis.
#     fig.subplots_adjust(right=0.75)

# ## Move the last y-axis spine over to the right by 20% of the width of the axes
#     axes[-1].spines['right'].set_position(('axes', 1.2))

# # To make the border of the right-most axis visible, we need to turn the frame
# # on. This hides the other plots, however, so we need to turn its fill off.
#     axes[-1].set_frame_on(True)
#     axes[-1].patch.set_visible(False)

# # And finally we get to plot things...
#     colors = (color1, color2, color3)
#     field = (timeseries1,timeseries2,timeseries3)
#     names = (varname1,varname2,varname3)
#     i=0
#     for ax, color in zip(axes, colors):
#         data = field[i]
# 	if (i==2 and height3 !=0):
# 		ax.plot(months,data,color=color,label=names[i]+' at pfull = '+str(int(nc.variables['pfull'][height3]))+' hPa')
# 	else:	
# 		ax.plot(months,data,color=color,label=names[i])
#         ax.set_ylabel(names[i])
# 	if (i<=1):
# 		ax.set_ylim(0,timeseries1.max()+1.)
#         ax.tick_params(axis='y',colors=color)
#         i+=1
#     axes[0].set_xlabel('month #')
#     plt.legend()
#     plt.title('Tropical Ocean Only')
#     plt.show()

def globavg_var_timeseries_total_and_land(outdir,testdir,model,area_array,varname,runmin,runmax,factor,landmask,minlat=-90.,maxlat=90.,select='all'):
    
#	""" Plots and returns timeseries of global average for specified variable """
# factor is needed to convert eg precip from kg/s to mm/day 

    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        if model=='isca':
		runnr="{0:04}".format(i)
	elif model=='gfdl':
		runnr="{0:03}".format(i)
        filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
        nc = Dataset(filename,mode='r')
        
	lats = nc.variables['lat'][:]
	lons = nc.variables['lon'][:]
        var=nc.variables[varname][0,:,:]*factor
	var = xr.DataArray(var,coords=[lats,lons],dims=['lat','lon'])

        if i==runmin:
		timeseries = [area_weighted_avg(var,area_array,landmask,'all_sfcs',minlat,maxlat)]
		timeseries_land = [area_weighted_avg(var,area_array,landmask,'land',minlat,maxlat)]
		timeseries_ocean = [area_weighted_avg(var,area_array,landmask,'ocean',minlat,maxlat)]

        else:
		timeseries.append(area_weighted_avg(var,area_array,landmask,'all_sfcs',minlat,maxlat))
		timeseries_land.append(area_weighted_avg(var,area_array,landmask,'land',minlat,maxlat))
		timeseries_ocean.append(area_weighted_avg(var,area_array,landmask,'ocean',minlat,maxlat))
	
    timeseries=np.asarray(timeseries)
    timeseries_land=np.asarray(timeseries_land)
    timeseries_ocean=np.asarray(timeseries_ocean)
    months=np.linspace(runmin,(runmax-1),timeseries.size)
    if select == 'all':
	    plt.plot(months,timeseries,'r',label='total')
	    plt.plot(months,timeseries_land,'g',label='land only')
	    plt.plot(months,timeseries_ocean,'b',label='ocean only')
	    plt.title('avg '+varname+' total and land/ocean only between '+str(minlat)+' - '+str(maxlat)+' N')
    
    if select == 'land':
	    plt.plot(months,timeseries_land,'g',label='land only')
	    plt.title('avg '+varname+' land only between '+str(minlat)+' - '+str(maxlat)+' N')
    elif select == 'ocean': 
	    plt.plot(months,timeseries_ocean,'b',label='ocean only')
	    plt.title('avg '+varname+' ocean only between '+str(minlat)+' - '+str(maxlat)+' N')
    elif select =='total':
	    plt.plot(months,timeseries,'r',label='total')
	    plt.title('avg '+varname+' total (land and ocean) between '+str(minlat)+' - '+str(maxlat)+' N')


    plt.xlabel('Month #')
    plt.ylabel('Area weighted average')
    plt.legend()
    plt.grid(b=True)
    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+varname+'_timeseries_'+str(runmin)+'-'+str(runmax), bbox_inches='tight')
    plt.close()
#    plt.show()    
    return(timeseries)



def globavg_var_timeseries_total_and_land_6hrly(outdir,testdir,model,area_array,varname,runmin,runmax,factor,landmask,minlat=-90.,maxlat=90.,select='all'):
    
#	""" Plots and returns timeseries of global average for specified variable """
# factor is needed to convert eg precip from kg/s to mm/day 

    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        if model=='isca':
		runnr="{0:04}".format(i)
	elif model=='gfdl':
		runnr="{0:03}".format(i)
        filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_6_hourly.nc'
        nc = Dataset(filename,mode='r')
        
	lats = nc.variables['lat'][:]
	lons = nc.variables['lon'][:]
              
        if i==runmin:
            var=xr.DataArray(nc.variables[varname][:])
        else:
            var_i=xr.DataArray(nc.variables[varname][:])
            var=xr.concat([var,var_i],'dim_0')
    
    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin)*4*30,dtype='datetime64[h]'))]
    var=xr.DataArray((var.values)*factor,coords=[time[0],lats,lons],dims=['time','lat','lon'])	


    for i in range (0,(runmax-runmin)*120):
        if i==0:
		timeseries = [area_weighted_avg(var[0,:,:],area_array,landmask,'all_sfcs',minlat,maxlat)]
		timeseries_land = [area_weighted_avg(var[0,:,:],area_array,landmask,'land',minlat,maxlat)]
		timeseries_ocean = [area_weighted_avg(var[0,:,:],area_array,landmask,'ocean',minlat,maxlat)]

        else:
		timeseries.append(area_weighted_avg(var[i,:,:],area_array,landmask,'all_sfcs',minlat,maxlat))
		timeseries_land.append(area_weighted_avg(var[i,:,:],area_array,landmask,'land',minlat,maxlat))
		timeseries_ocean.append(area_weighted_avg(var[i,:,:],area_array,landmask,'ocean',minlat,maxlat))
	
    timeseries=np.asarray(timeseries)
    timeseries_land=np.asarray(timeseries_land)
    timeseries_ocean=np.asarray(timeseries_ocean)
    months=np.linspace(runmin,(runmax-1),timeseries.size)
    if select == 'all':
	    plt.plot(months,timeseries,'r',label='total')
	    plt.plot(months,timeseries_land,'g',label='land only')
	    plt.plot(months,timeseries_ocean,'b',label='ocean only')
	    plt.title('avg '+varname+' total and land/ocean only between '+str(minlat)+' and '+str(maxlat)+' N')
    
    if select == 'land':
	    plt.plot(months,timeseries_land,'g',label='land only')
	    plt.title('avg '+varname+' land only between '+str(minlat)+' - '+str(maxlat)+' N')
    elif select == 'ocean': 
	    plt.plot(months,timeseries_ocean,'b',label='ocean only')
	    plt.title('avg '+varname+' ocean only between '+str(minlat)+' - '+str(maxlat)+' N')
    elif select =='total':
	    plt.plot(months,timeseries,'r',label='total')
	    plt.title('avg '+varname+' total (land and ocean) between '+str(minlat)+' and '+str(maxlat)+' N')


    plt.xlabel('6 hourly units # ')
    plt.ylabel('Area weighted average')
    plt.legend()
    plt.grid(b=True)
    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+varname+'_timeseries_'+str(runmin)+'-'+str(runmax)+'_'+select+'_between_'+str(minlat)+'_and_'+str(maxlat)+'.png', bbox_inches='tight')
    plt.close()
#    plt.show()    
    return(timeseries)



def globavg_var_timeseries_total_and_land_6hrly_severalvars(outdir,testdir,model,area_array,varname1,varname2,varname3,varname4,runmin,runmax,factor1,factor2,factor3,factor4,landmask,minlat=-90.,maxlat=90.,select='all'):
    
#	""" Plots and returns timeseries of global average for specified variable """
# factor is needed to convert eg precip from kg/s to mm/day 

    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        if model=='isca':
		runnr="{0:04}".format(i)
	elif model=='gfdl':
		runnr="{0:03}".format(i)
        filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_6_hourly.nc'
        nc = Dataset(filename,mode='r')
        
	lats = nc.variables['lat'][:]
	lons = nc.variables['lon'][:]
              
        if i==runmin:
            var1=xr.DataArray(nc.variables[varname1][:])
            var2=xr.DataArray(nc.variables[varname2][:])
            var3=xr.DataArray(nc.variables[varname3][:])
            var4=xr.DataArray(nc.variables[varname4][:])

        else:
            var1_i=xr.DataArray(nc.variables[varname1][:])
            var1=xr.concat([var1,var1_i],'dim_0')
            var2_i=xr.DataArray(nc.variables[varname2][:])
            var2=xr.concat([var2,var2_i],'dim_0')
            var3_i=xr.DataArray(nc.variables[varname3][:])
            var3=xr.concat([var3,var3_i],'dim_0')
            var4_i=xr.DataArray(nc.variables[varname4][:])
            var4=xr.concat([var4,var4_i],'dim_0')

    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin)*4*30,dtype='datetime64[h]'))]
    var1=xr.DataArray((var1.values)*factor1,coords=[time[0],lats,lons],dims=['time','lat','lon'])	

    var2=xr.DataArray((var2.values)*factor2,coords=[time[0],lats,lons],dims=['time','lat','lon'])

    var3=xr.DataArray((var3.values)*factor3,coords=[time[0],lats,lons],dims=['time','lat','lon'])

    var4=xr.DataArray((var4.values)*factor4,coords=[time[0],lats,lons],dims=['time','lat','lon'])


    for i in range (0,(runmax-runmin)*120):
        if i==0:
		timeseries1 = [area_weighted_avg(var1[0,:,:],area_array,landmask,'all_sfcs',minlat,maxlat)]
		timeseries1_land = [area_weighted_avg(var1[0,:,:],area_array,landmask,'land',minlat,maxlat)]
		timeseries1_ocean = [area_weighted_avg(var1[0,:,:],area_array,landmask,'ocean',minlat,maxlat)]

		timeseries2 = [area_weighted_avg(var2[0,:,:],area_array,landmask,'all_sfcs',minlat,maxlat)]
		timeseries2_land = [area_weighted_avg(var2[0,:,:],area_array,landmask,'land',minlat,maxlat)]
		timeseries2_ocean = [area_weighted_avg(var2[0,:,:],area_array,landmask,'ocean',minlat,maxlat)]


		timeseries3 = [area_weighted_avg(var3[0,:,:],area_array,landmask,'all_sfcs',minlat,maxlat)]
		timeseries3_land = [area_weighted_avg(var3[0,:,:],area_array,landmask,'land',minlat,maxlat)]
		timeseries3_ocean = [area_weighted_avg(var3[0,:,:],area_array,landmask,'ocean',minlat,maxlat)]


		timeseries4 = [area_weighted_avg(var4[0,:,:],area_array,landmask,'all_sfcs',minlat,maxlat)]
		timeseries4_land = [area_weighted_avg(var4[0,:,:],area_array,landmask,'land',minlat,maxlat)]
		timeseries4_ocean = [area_weighted_avg(var4[0,:,:],area_array,landmask,'ocean',minlat,maxlat)]

        else:
		timeseries1.append(area_weighted_avg(var1[i,:,:],area_array,landmask,'all_sfcs',minlat,maxlat))
		timeseries1_land.append(area_weighted_avg(var1[i,:,:],area_array,landmask,'land',minlat,maxlat))
		timeseries1_ocean.append(area_weighted_avg(var1[i,:,:],area_array,landmask,'ocean',minlat,maxlat))


		timeseries2.append(area_weighted_avg(var2[i,:,:],area_array,landmask,'all_sfcs',minlat,maxlat))
		timeseries2_land.append(area_weighted_avg(var2[i,:,:],area_array,landmask,'land',minlat,maxlat))
		timeseries2_ocean.append(area_weighted_avg(var2[i,:,:],area_array,landmask,'ocean',minlat,maxlat))


		timeseries3.append(area_weighted_avg(var3[i,:,:],area_array,landmask,'all_sfcs',minlat,maxlat))
		timeseries3_land.append(area_weighted_avg(var3[i,:,:],area_array,landmask,'land',minlat,maxlat))
		timeseries3_ocean.append(area_weighted_avg(var3[i,:,:],area_array,landmask,'ocean',minlat,maxlat))


		timeseries4.append(area_weighted_avg(var4[i,:,:],area_array,landmask,'all_sfcs',minlat,maxlat))
		timeseries4_land.append(area_weighted_avg(var4[i,:,:],area_array,landmask,'land',minlat,maxlat))
		timeseries4_ocean.append(area_weighted_avg(var4[i,:,:],area_array,landmask,'ocean',minlat,maxlat))

	
    timeseries1=np.asarray(timeseries1)
    timeseries1_land=np.asarray(timeseries1_land)
    timeseries1_ocean=np.asarray(timeseries1_ocean)

	
    timeseries2=np.asarray(timeseries2)
    timeseries2_land=np.asarray(timeseries2_land)
    timeseries2_ocean=np.asarray(timeseries2_ocean)

	
    timeseries3=np.asarray(timeseries3)
    timeseries3_land=np.asarray(timeseries3_land)
    timeseries3_ocean=np.asarray(timeseries3_ocean)

	
    timeseries4=np.asarray(timeseries4)
    timeseries4_land=np.asarray(timeseries4_land)
    timeseries4_ocean=np.asarray(timeseries4_ocean)

    months=np.linspace(runmin,(runmax-1),timeseries1.size)
    
    if select == 'land':

	    fig, ax = plt.subplots()
# Twin the x-axis twice to make independent y-axes.
	    axes = [ax, ax.twinx(), ax.twinx()] 
# Make some space on the right side for the extra y-axis.
	    fig.subplots_adjust(right=0.75)
## Move the last y-axis spine over to the right by 20% of the width of the axes
	    axes[-1].spines['right'].set_position(('axes', 1.2))

# To make the border of the right-most axis visible, we need to turn the frame
# on. This hides the other plots, however, so we need to turn its fill off.
	    axes[-1].set_frame_on(True)
	    axes[-1].patch.set_visible(False)

	    colors = ('g','b','r','k')
	    field = (timeseries1_land,timeseries2_land,timeseries3_land,timeseries4_land)
	    names = (varname1,varname2,varname3,varname4)
	    i=0
	    for i in range(0,4):
		    if (i<=1):
			    axes[0].plot(months,field[i],color=colors[i],label=names[i])
			    axes[0].set_ylabel(names[0]+', '+names[1])
		    else: 
			    axes[i-1].set_ylabel(names[i])
			    axes[i-1].plot(months,field[i],color=colors[i],label=names[i])
		    # if (i<1):
		    # 	    ax.set_ylim(0,timeseries1.max()+1.)
		    ax.tick_params(axis='y',colors=color)
		    axes[0].set_xlabel('month #')
		    plt.legend()
		    plt.title('avg land only between '+str(minlat)+' and '+str(maxlat)+' N')
    # elif select == 'ocean': 
    # 	    plt.plot(months,timeseries1_ocean,'g',label=varname1)
    # 	    plt.plot(months,timeseries2_ocean,'b',label=varname2)
    # 	    plt.plot(months,timeseries3_ocean,'r',label=varname3)
    # 	    plt.plot(months,timeseries4_ocean,'k',label=varname4)
    # 	    plt.title('avg ocean only between '+str(minlat)+' and '+str(maxlat)+' N')
    # elif select =='total':
    # 	    plt.plot(months,timeseries1,'g',label=varname1)
    # 	    plt.plot(months,timeseries2,'b',label=varname2)
    # 	    plt.plot(months,timeseries3,'r',label=varname3)
    # 	    plt.plot(months,timeseries4,'k',label=varname4)
    # 	    plt.title('avg total between '+str(minlat)+' and '+str(maxlat)+' N')

    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/P_E_bucket_t_timeseries_'+str(runmin)+'-'+str(runmax)+'_'+select+'_between_'+str(minlat)+'_and_'+str(maxlat)+'.png', bbox_inches='tight')
    plt.close()
#    plt.show()    




def globavg_var_timeseries_selected_points__6hrly_severalvars(outdir,testdir,model,area_array,varname1,varname2,varname3,varname4,varname5,runmin,runmax,factor1,factor2,factor3,factor4,factor5,landmask,precipitation_avg,minlat=-90.,maxlat=90.,maxormin='max',level5=39):
    
#	""" Plots and returns timeseries of global average for specified variable """
# factor is needed to convert eg precip from kg/s to mm/day 

# currently not using minlat and maxlat!

    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        if model=='isca':
		runnr="{0:04}".format(i)
	elif model=='gfdl':
		runnr="{0:03}".format(i)
        filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_6_hourly.nc'
        nc = Dataset(filename,mode='r')
        
	lats = nc.variables['lat'][:]
	lons = nc.variables['lon'][:]
              
        if i==runmin:
            var1=xr.DataArray(nc.variables[varname1][:])
            var2=xr.DataArray(nc.variables[varname2][:])
            var3=xr.DataArray(nc.variables[varname3][:])
            var4=xr.DataArray(nc.variables[varname4][:])
            var5=xr.DataArray(nc.variables[varname5][:,level5,:,:])

        else:
            var1_i=xr.DataArray(nc.variables[varname1][:])
            var1=xr.concat([var1,var1_i],'dim_0')
            var2_i=xr.DataArray(nc.variables[varname2][:])
            var2=xr.concat([var2,var2_i],'dim_0')
            var3_i=xr.DataArray(nc.variables[varname3][:])
            var3=xr.concat([var3,var3_i],'dim_0')
            var4_i=xr.DataArray(nc.variables[varname4][:])
            var4=xr.concat([var4,var4_i],'dim_0')
            var5_i=xr.DataArray(nc.variables[varname5][:,level5,:,:])
            var5=xr.concat([var5,var5_i],'dim_0')
	
    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin)*4*30,dtype='datetime64[h]'))]
    var1=xr.DataArray((var1.values)*factor1,coords=[time[0],lats,lons],dims=['time','lat','lon'])	

    var2=xr.DataArray((var2.values)*factor2,coords=[time[0],lats,lons],dims=['time','lat','lon'])

    var3=xr.DataArray((var3.values)*factor3,coords=[time[0],lats,lons],dims=['time','lat','lon'])

    var4=xr.DataArray((var4.values)*factor4,coords=[time[0],lats,lons],dims=['time','lat','lon'])
    var5=xr.DataArray((var5.values)*factor5,coords=[time[0],lats,lons],dims=['time','lat','lon'])


    if (maxormin == 'max'):
	    coords = np.where((precipitation_avg.where(landmask == 1.)) == (precipitation_avg.where(landmask == 1.)).max()) # gives a tuple
	    print (precipitation_avg.where(landmask == 1.)).where((precipitation_avg.where(landmask == 1.) == (precipitation_avg.where(landmask == 1.)).max()),drop=True)
    elif (maxormin == 'min'):
	    coords = np.where((precipitation_avg.where(landmask == 1.)) == (precipitation_avg.where(landmask == 1.)).min()) # gives a tuple# gives a tuple
	    print (precipitation_avg.where(landmask == 1.)).where((precipitation_avg.where(landmask == 1.) == (precipitation_avg.where(landmask == 1.)).min()),drop=True)
    latP = (coords[0])[0]
    lonP = (coords[1])[0]


# NB: precipitation_avg.where(precipitation_avg == precipitation_avg.max(), drop=True)
# shows the Pmax, and corresponding lat and lon 
    timeseries1=np.empty((runmax-runmin)*120)
    timeseries2=np.empty((runmax-runmin)*120)
    timeseries3=np.empty((runmax-runmin)*120)
    timeseries4=np.empty((runmax-runmin)*120)
    timeseries5=np.empty((runmax-runmin)*120)

    
    for i in range (0,(runmax-runmin)*120):
	    timeseries1[i] = var1[i,latP,lonP]
	    timeseries2[i] = var2[i,latP,lonP]
	    timeseries3[i] = var3[i,latP,lonP]
	    timeseries4[i] = var4[i,latP,lonP]
	    timeseries5[i] = var5[i,latP,lonP]

    months=np.linspace(runmin,(runmax-1),timeseries1.size)
    hours = np.linspace(0,(runmax-runmin-1)*30.*24.,timeseries1.size)
    fig, ax = plt.subplots()
# Twin the x-axis twice to make independent y-axes.
    axes = [ax, ax.twinx(), ax.twinx(), ax.twinx()] 
# Make some space on the right side for the extra y-axis.
    fig.subplots_adjust(right=0.60)
## Move the last y-axis spine over to the right by 20% of the width of the axes
    axes[-1].spines['right'].set_position(('axes', 1.2))
    axes[-2].spines['right'].set_position(('axes', 1.4))

# To make the border of the right-most axis visible, we need to turn the frame
# on. This hides the other plots, however, so we need to turn its fill off.
    # axes[-1].set_frame_on(True)
    # axes[-1].patch.set_visible(False)
    axes[-2].set_frame_on(True)
    axes[-2].patch.set_visible(False)


    colors = ('g','b','r','k','Orange')
    field = (timeseries1,timeseries2,timeseries3,timeseries4,timeseries5)
    names = (varname1,varname2,varname3,varname4,varname5)
    for i in range(0,5):
	    if (i<=1):
		    axes[0].plot(hours,field[i],color=colors[i],label=names[i])
		    axes[0].set_ylabel(names[0]+', '+names[1])
		    axes[0].legend()
	    else: 
		    axes[i-1].set_ylabel(names[i])
		    axes[i-1].plot(hours,field[i],color=colors[i],label=names[i])
		    # if (i<1):
		    # 	    ax.set_ylim(0,timeseries1.max()+1.)
		    axes[i-1].yaxis.label.set_color(colors[i])

	    ax.tick_params(axis='y',color=colors[i])
	    axes[0].set_xlabel('hour')
    plt.title('timeseries in P_'+maxormin)

    plt.show()
#    plt.close()






def globavg_var_timeseries_total_and_land_perturbed(testdir,model,area_array,varname,runmin,runmax,factor,landmask,spinup_dir,spinup_model,spinup_runmin, spinup_runmax,minlat=-90.,maxlat=90.,select='all'):
    
#	""" Plots and returns timeseries of global average for specified variable """
# factor is needed to convert eg precip from kg/s to mm/day 


	for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
		if (model=='isca') or (model == 'Isca'):
			runnr="{0:04}".format(i)
		elif model=='gfdl':
			runnr="{0:03}".format(i)
		filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_monthly.nc'


		nc = Dataset(filename,mode='r')   
		lats = nc.variables['lat'][:]
		lons = nc.variables['lon'][:]
		var=nc.variables[varname][0,:,:]*factor
	# var_land = (xr.DataArray(var)).where(landmask==1.)
	# var_ocean = (xr.DataArray(var)).where(landmask!=1.)
		var = xr.DataArray(var,coords=[lats,lons],dims=['lat','lon'])

		if i==runmin:
			timeseries = [area_weighted_avg(var,area_array,landmask,'all_sfcs',minlat,maxlat)]
			timeseries_land = [area_weighted_avg(var,area_array,landmask,'land',minlat,maxlat)]
			timeseries_ocean = [area_weighted_avg(var,area_array,landmask,'ocean',minlat,maxlat)]

		else:
			timeseries.append(area_weighted_avg(var,area_array,landmask,'all_sfcs',minlat,maxlat))
			timeseries_land.append(area_weighted_avg(var,area_array,landmask,'land',minlat,maxlat))
			timeseries_ocean.append(area_weighted_avg(var,area_array,landmask,'ocean',minlat,maxlat))


#same for spin up 
	for i in range (spinup_runmin,spinup_runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
		if (spinup_model=='isca') or (spinup_model=='Isca'):
			runnr="{0:04}".format(i)
		elif spinup_model=='gfdl':
			runnr="{0:03}".format(i)
		filename = '/scratch/mp586/'+spinup_dir+'/run'+runnr+'/atmos_monthly.nc'
		nc = Dataset(filename,mode='r')
        
		var=nc.variables[varname][0,:,:]*factor
		
		var = xr.DataArray(var,coords=[lats,lons],dims=['lat','lon'])


		if i==runmin:
			ts_spinup = [area_weighted_avg(var,area_array,landmask,'all_sfcs',minlat,maxlat)]
			ts_spinup_land = [area_weighted_avg(var,area_array,landmask,'land',minlat,maxlat)]
			ts_spinup_ocean = [area_weighted_avg(var,area_array,landmask,'ocean',minlat,maxlat)]
			
		else:
			ts_spinup.append(area_weighted_avg(var,area_array,landmask,'all_sfcs',minlat,maxlat))
			ts_spinup_land.append(area_weighted_avg(var,area_array,landmask,'land',minlat,maxlat))
			ts_spinup_ocean.append(area_weighted_avg(var,area_array,landmask,'ocean',minlat,maxlat))

	timeseries=xr.DataArray(timeseries)
	timeseries_land=xr.DataArray(timeseries_land)
	timeseries_ocean=xr.DataArray(timeseries_ocean)

	ts_spinup=xr.DataArray(ts_spinup)
	ts_spinup_land=xr.DataArray(ts_spinup_land)
	ts_spinup_ocean=xr.DataArray(ts_spinup_ocean) 

	timeseries = xr.concat([ts_spinup,timeseries],'dim_0')
	timeseries_land = xr.concat([ts_spinup_land,timeseries_land],'dim_0')
	timeseries_ocean = xr.concat([ts_spinup_ocean,timeseries_ocean],'dim_0')


	figure = plt.plot()
	months=np.linspace(runmin,((spinup_runmax+runmax)-1),timeseries.size)
	if select == 'all':
		plt.plot(months,timeseries,'r',label='total')
		plt.plot(months,timeseries_land,'g',label='land only')
		plt.plot(months,timeseries_ocean,'b',label='ocean only')
		plt.title('avg '+varname+' total and land/ocean only between '+str(minlat)+' - '+str(maxlat)+' N')
    
	if select == 'land':
		plt.plot(months,timeseries_land,'g',label='land only')
		plt.title('avg '+varname+' land only between '+str(minlat)+' - '+str(maxlat)+' N')
	elif select == 'ocean': 
		plt.plot(months,timeseries_ocean,'b',label='ocean only')
		plt.title('avg '+varname+' ocean only between '+str(minlat)+' - '+str(maxlat)+' N')
	elif select =='total':
		plt.plot(months,timeseries,'r',label='total')
		plt.title('avg '+varname+' total (land and ocean) between '+str(minlat)+' - '+str(maxlat)+' N')

	plt.xlabel('Month #')
	plt.ylabel('Global average')
	plt.legend()
	plt.grid(b=True)
	plt.show()    
	return(timeseries)




def seasonal_surface_variable(testdir,model,runmin,runmax,varname,units,factor=1.,level=None): 

    print(varname+' '+str(factor))

    plt.close()
    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        if model=='isca':
		runnr="{0:04}".format(i)
	elif model=='gfdl':
		runnr="{0:03}".format(i)
        filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
        nc = Dataset(filename,mode='r')
              
        if i==runmin:
		if (level==None) or (level=='all'):
			var=xr.DataArray(nc.variables[varname][:]) # only monthly avg for month i
		else:
			var=xr.DataArray(nc.variables[varname][:,level,:,:]) # only monthly avg for month i

        else:
		if (level==None) or (level=='all'): 
			var_i=xr.DataArray(nc.variables[varname][:])
			var=xr.concat([var,var_i],'dim_0')
		else:
			var_i=xr.DataArray(nc.variables[varname][:,level,:,:])
			var=xr.concat([var,var_i],'dim_0')

    if level == 'all':
	    var = var.sum('dim_1') # height integral
    
    lons= nc.variables['lon'][:]
    lats= nc.variables['lat'][:]
    
    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))]
    var=xr.DataArray((var.values)*factor,coords=[time[0],lats,lons],dims=['time','lat','lon'])
    var_avg=var.mean(dim='time')
    var_seasonal_avg=var.groupby('time.season').mean('time') 
    var_month_avg=var.groupby('time.month').mean('time')

    return(var,var_avg,var_seasonal_avg,var_month_avg,time)

def seasonal_surface_variable_6hrly(testdir,model,runmin,runmax,varname,units,factor=1.,level=None): 

    print(varname+' '+str(factor))

    plt.close()
    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        if model=='isca':
		runnr="{0:04}".format(i)
	elif model=='gfdl':
		runnr="{0:03}".format(i)
        filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_6_hourly.nc'
        nc = Dataset(filename,mode='r')


              
        if i==runmin:
		if (level==None) or (level=='all'):
			var=xr.DataArray(nc.variables[varname][:])
		else:
			var=xr.DataArray(nc.variables[varname][:,level,:,:])

        else:
		if (level==None) or (level=='all'): 
			var_i=xr.DataArray(nc.variables[varname][:])
			var=xr.concat([var,var_i],'dim_0')
		else:
			var_i=xr.DataArray(nc.variables[varname][:,level,:,:])
			var=xr.concat([var,var_i],'dim_0')

    if level == 'all':
	    var = var.sum('dim_1') # height integral
    
    lons= nc.variables['lon'][:]
    lats= nc.variables['lat'][:]

    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin)*4*30,dtype='datetime64[h]'))]
    var=xr.DataArray((var.values)*factor,coords=[time[0],lats,lons],dims=['time','lat','lon'])
    var_avg=var.mean(dim='time')

    return(var,var_avg,time)


def seasonal_4D_variable_6hrly(testdir,model,runmin,runmax,varname,units): 

    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        if model=='isca':
		runnr="{0:04}".format(i)
	elif model=='gfdl':
		runnr="{0:03}".format(i)
        filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_6_hourly.nc'
        nc = Dataset(filename,mode='r')
              
        if i==runmin:
            var=xr.DataArray(nc.variables[varname][:])
        else:
            var_i=xr.DataArray(nc.variables[varname][:])
            var=xr.concat([var,var_i],'dim_0')
    
    lons= nc.variables['lon'][:]
    lats= nc.variables['lat'][:]
    pres_levs= nc.variables['pfull'][:]
 
    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin)*4*30,dtype='datetime64[h]'))]
    var=xr.DataArray(var.values,coords=[time[0],pres_levs,lats,lons],dims=['time','pres_lev','lat','lon'])
    var_avg=var.mean(dim='time')

    return(var,var_avg,time)


def seasonal_4D_variable(testdir,model,runmin,runmax,varname,units): 

    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        if model=='isca':
		runnr="{0:04}".format(i)
	elif model=='gfdl':
		runnr="{0:03}".format(i)
        filename = '/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
        nc = Dataset(filename,mode='r')
              
        if i==runmin:
            var=xr.DataArray(nc.variables[varname][:]) # only monthly avg for month i
        else:
            var_i=xr.DataArray(nc.variables[varname][:])
            var=xr.concat([var,var_i],'dim_0')
    
    lons= nc.variables['lon'][:]
    lats= nc.variables['lat'][:]
    pres_levs= nc.variables['pfull'][:]

    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))]
    var=xr.DataArray(var.values,coords=[time[0],pres_levs,lats,lons],dims=['time','pres_lev','lat','lon'])
    var_avg=var.mean(dim='time')
    var_seasonal_avg=var.groupby('time.season').mean('time') 
    var_month_avg=var.groupby('time.month').mean('time')
    var_annual_avg=var.groupby('time.year').mean('time')

    return(var,var_avg,var_seasonal_avg,var_month_avg,var_annual_avg,time)

# def zonavg_plot(field,title,varname,units): #field needs to have dimensions lat|lon

#     field=xr.DataArray(field)
#     field_zonavg=field.mean(dim='lon',keep_attrs=True)
    
#     plt.close()
#     plt.plot(field['lat'],field_zonavg)
#     plt.title(title)
#     plt.ylabel(varname+' ('+units+')')
#     plt.xlabel('Latitude')
#     plt.show()

# def several_vars_zonalavg(field1,varname1,field2,varname2,field3,varname3): 

# # fields needs to have dimensions lat|lon
# # same y axis for all vars

#     field1=xr.DataArray(field1)
#     field1_zonavg=field1.mean(dim='lon',keep_attrs=True)
#     field2=xr.DataArray(field2)
#     field2_zonavg=field2.mean(dim='lon',keep_attrs=True)   
#     field3=xr.DataArray(field3)
#     field3_zonavg=field3.mean(dim='lon',keep_attrs=True)
    
#     plt.close()
#     plt.plot(field1['lat'],field1_zonavg,label=varname1)
#     plt.plot(field1['lat'],field2_zonavg,label=varname2)
#     plt.plot(field1['lat'],field3_zonavg,label=varname3)
#     plt.legend()
#     plt.xlabel('latitude')
#     plt.show()

def several_vars_zonalavg2(field1,varname1,color1,field2,varname2,color2,field3,varname3,color3,title):
# with 3 different y axes
# adapted from http://stackoverflow.com/questions/7733693/matplotlib-overlay-plots-with-different-scales


    field1=xr.DataArray(field1)
    field1_zonavg=field1.mean(dim='lon',keep_attrs=True)
    field2=xr.DataArray(field2)
    field2_zonavg=field2.mean(dim='lon',keep_attrs=True)   
    field3=xr.DataArray(field3)
    field3_zonavg=field3.mean(dim='lon',keep_attrs=True)

    fig, ax = plt.subplots()

# Twin the x-axis twice to make independent y-axes.
    axes = [ax, ax.twinx(), ax.twinx()] 

# Make some space on the right side for the extra y-axis.
    fig.subplots_adjust(right=0.75)

## Move the last y-axis spine over to the right by 20% of the width of the axes
    axes[-1].spines['right'].set_position(('axes', 1.2))

# To make the border of the right-most axis visible, we need to turn the frame
# on. This hides the other plots, however, so we need to turn its fill off.
    axes[-1].set_frame_on(True)
    axes[-1].patch.set_visible(False)

# And finally we get to plot things...
    colors = (color1, color2, color3)
    field = (field1_zonavg,field2_zonavg,field3_zonavg)
    lats=field1['lat']
    names = (varname1,varname2,varname3)
    i=0
    for ax, color in zip(axes, colors):
        data = field[i]
        ax.plot(lats,data,color=color)
        ax.set_ylabel(names[i])
	if (i>=1):
		ax.set_ylim(0,field2_zonavg.max()+1.)
        ax.tick_params(axis='y',colors=color)
        i+=1
    axes[0].set_xlabel('latitude')
    plt.title(title)
    plt.show()

def tropics_plot(lats,lons,field,units,title):
    plt.close()
    trop_minindex=np.asarray(np.where(lats>=-30.))[0,0]
    trop_maxreverseindex=np.asarray(np.where(lats[::-1]<=30.))[0,0] 
    tropical_lats=lats[trop_minindex:(lats.size-trop_maxreverseindex)]
    lon_0 = lons.mean() # trick to find central lon lat (lon from 0 to 360)
    lat_0 = lats.mean() # ~= 0 because lat from -90 to 90



    m = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=30,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')
    m.drawcoastlines()
    m.drawparallels(np.arange(-30.,31.,30.))
    m.drawmeridians(np.arange(0.,361.,60.))


# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
    lon, lat = np.meshgrid(lons, tropical_lats)
    xi, yi = m(lon, lat)

   
    cs = m.pcolor(xi,yi,field.sel(lat=tropical_lats))
    
    
# Add Colorbar
    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)

# Add Title
    plt.title(title)


    plt.show()
    return(cs)


def squareland_plot(minlat,maxlat,array,units,title,palette,contourson = False):
# kept most of the comments in this function

# plotting only the zonal average next to the map 
    plt.close()

    lats=array.lat
    lons=array.lon

    
    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] # there might be a better way of doing this!
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

    # sqlatindex=np.asarray(np.where(lats>sq_minlat))[0,0]
    # sqlatreverseindex=np.asarray(np.where(lats[::-1]<sq_maxlat))[0,0] # there might be a better way of doing this!
    # square_lats=lats[sqlatindex:(lats.size-sqlatreverseindex)]

    # sqlonindex=np.asarray(np.where(lons>sq_minlon))[0,0]
    # sqlonreverseindex=np.asarray(np.where(lons[::-1]<sq_maxlon))[0,0] # there might be a better way of doing this!
    # square_lons=lons[sqlonindex:(lons.size-sqlonreverseindex)+1]


    fig = plt.figure()
    ax1 = plt.subplot2grid((5,8), (0,1), colspan = 5, rowspan = 3)


    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# draw parallels and meridians.

    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    array, lons_cyclic = addcyclic(array, lons)
    array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
    zonavg_thin = array.mean(dim='lon')
    meravg_thin = array.mean(dim='lat')
# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.

 
    lon, lat = np.meshgrid(lons_cyclic, selected_lats)
    xi, yi = m(lon, lat)

    dlons = lons[100] - lons[99]
    dlats = lats[60] - lats[59]

#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    minval = np.absolute(array.min())
    maxval = np.absolute(array.max())
    
    if maxval >= minval:
	    minval = - maxval
    else: 
	    maxval = minval
	    minval = - minval


    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),
		      norm=MidpointNormalize(midpoint=0.),cmap='BrBG',
		      vmin=minval, vmax=maxval)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), 
		      norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r, 
		      vmin = 273.15-(maxval-273.15),
		      vmax=maxval)
    elif palette=='tempdiff': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), 
		      norm=MidpointNormalize(midpoint=0), cmap=plt.cm.RdBu_r, 
		      vmin = -maxval,
		      vmax = maxval)
    elif palette=='fromwhite': 
	    pal = plt.cm.Blues
	    pal.set_under('w',None)
	    cs = m.pcolormesh(xi,yi,array.sel(lat=selected_lats),cmap=pal,vmin=0,vmax=maxval)

    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))




#################### lon from 0 to 360 ####################
#     m = Basemap(projection='cyl',llcrnrlat=minlat,urcrnrlat=maxlat,\
#             llcrnrlon=-1.,urcrnrlon=360.,resolution='c')
# #    m.drawcoastlines()
# #m.fillcontinents(color='coral',lake_color='aqua')
# # draw parallels and meridians.
#     m.drawparallels(np.arange(minlat,maxlat+1.,30.),labels=[1,0,0,0])
#     m.drawmeridians(np.arange(0.,361.,60.),labels=[0,0,0,1])


# # Because our lon and lat variables are 1D, 
# # use meshgrid to create 2D arrays 
# # Not necessary if coordinates are already in 2D arrays.
#     lon, lat = np.meshgrid(lons, selected_lats)
#     xi, yi = m(lon, lat)

# #latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

#     if palette=='rainnorm':
#         cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG')#,latlon=True)
#     elif palette == 'raindefault':
#         cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG)#,latlon=True)
#     elif palette=='temp': 
#         cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.seismic)#,latlon=True)
#     else:
#         cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))#,latlon=True)
############################################################


    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)

    landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land_square.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    
    square_lons,square_lats=(xr.DataArray(np.meshgrid(lons, lats))).where(landmask==1.)
    square_lons = square_lons + 180. + dlons

    m.drawgreatcircle(square_lons.min(), square_lats.min(), square_lons.min(), square_lats.max()+dlats, del_s=100., color='b')
    m.drawgreatcircle(square_lons.max(), square_lats.min(), square_lons.max(), square_lats.max()+dlats, del_s=100., color='b')

    lats = [square_lats.min(),square_lats.min()] 
    lons = [square_lons.min(),square_lons.max()] 
    x, y = m(lons, lats) 
    m.plot(x, y, color='b') 

    lats = [square_lats.max()+dlats,square_lats.max()+dlats] 
    lons = [square_lons.min(),square_lons.max()] 
    x, y = m(lons, lats) 
    m.plot(x, y, color='b') 


    if contourson == True:  # add contours 
	    cont = m.contour(xi,yi,array,4,cmap='PuBu_r',
			     linewidth=5)
	    plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=12)

    plt.title(title)
   
    ax2 = plt.subplot2grid((5,8), (0,7), rowspan = 3)
    # if palette=='temp':
    # 	    plt.plot((zonavg_thin-273.15),array['lat'])
    # else: 
    plt.plot(zonavg_thin,array['lat'])

    plt.ylabel('Latitude')
    plt.xlabel(title+' ('+units+')')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.invert_xaxis()

    meravg_thin_land = (array.sel(lat=slice(square_lats.min(),square_lats.max()))).mean(dim='lat')

    ax3 = plt.subplot2grid((5,8), (4,1), colspan = 4)

    plt.plot(array['lon'],meravg_thin)
    plt.xlabel('Longitude')
    plt.ylabel(title+' ('+units+') 30S-30N only')
    
    plt.show()

    print(square_lons.min(),square_lons.max(),square_lats.min(),square_lats.max()+dlats)

def squareland_plot_forpaper(minlat,maxlat,array,units,title,palette,contourson=False):
# kept most of the comments in this function

# plotting only the zonal average next to the map 
    plt.close()

    lats=array.lat
    lons=array.lon

    
    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] # there might be a better way of doing this!
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

    # sqlatindex=np.asarray(np.where(lats>sq_minlat))[0,0]
    # sqlatreverseindex=np.asarray(np.where(lats[::-1]<sq_maxlat))[0,0] # there might be a better way of doing this!
    # square_lats=lats[sqlatindex:(lats.size-sqlatreverseindex)]

    # sqlonindex=np.asarray(np.where(lons>sq_minlon))[0,0]
    # sqlonreverseindex=np.asarray(np.where(lons[::-1]<sq_maxlon))[0,0] # there might be a better way of doing this!
    # square_lons=lons[sqlonindex:(lons.size-sqlonreverseindex)+1]


    fig = plt.figure()
    ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)


    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# draw parallels and meridians.

    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    array, lons_cyclic = addcyclic(array, lons)
    array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
    zonavg_thin = array.mean(dim='lon')
    meravg_thin = array.mean(dim='lat')
# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.

 
    lon, lat = np.meshgrid(lons_cyclic, selected_lats)
    xi, yi = m(lon, lat)

    dlons = lons[100] - lons[99]
    dlats = lats[60] - lats[59]

#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    minval = np.absolute(array.min())
    maxval = np.absolute(array.max())
    
    if maxval >= minval:
	    minval = - maxval
    else: 
	    maxval = minval
	    minval = - minval


    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',vmin=minval, vmax=maxval)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), 
		      norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r,
		      vmin = 273.15-(maxval-273.15),vmax=maxval)
        # cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.RdBu_r)
    elif palette=='fromwhite': 
	    pal = 'Blues'
	    #pal.set_under('w',None)
	    cs = m.pcolormesh(xi,yi,array.sel(lat=selected_lats),cmap=pal,vmin=0,vmax=maxval)

    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))

    if contourson == True:  # add contours 
	    cont = m.contour(xi,yi,array,4,cmap='PuBu_r',
			     linewidth=5)
	    plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=12)

    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)

    landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land_square.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    
    square_lons,square_lats=(xr.DataArray(np.meshgrid(lons, lats))).where(landmask==1.)
    square_lons = square_lons + 180. + dlons

    m.drawgreatcircle(square_lons.min(), square_lats.min(), square_lons.min(), square_lats.max()+dlats, del_s=100., color='b')
    m.drawgreatcircle(square_lons.max(), square_lats.min(), square_lons.max(), square_lats.max()+dlats, del_s=100., color='b')

    lats = [square_lats.min(),square_lats.min()] 
    lons = [square_lons.min(),square_lons.max()] 
    x, y = m(lons, lats) 
    m.plot(x, y, color='b') 

    lats = [square_lats.max()+dlats,square_lats.max()+dlats] 
    lons = [square_lons.min(),square_lons.max()] 
    x, y = m(lons, lats) 
    m.plot(x, y, color='b') 
    plt.title(title)
   
    ax2 = plt.subplot2grid((5,7), (0,6), rowspan = 3)
    # if palette=='temp':
    # 	    plt.plot((zonavg_thin-273.15),array['lat'])
    # else: 
    plt.plot(zonavg_thin,array['lat'])

    plt.ylabel('Latitude')
    plt.xlabel(title+' ('+units+')')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
#    ax2.invert_xaxis()

    plt.show()

def squareland_plot_minuszonavg(minlat,maxlat,array,units,title,palette,zonavgtitle):

    plt.close()

    lats=array.lat
    lons=array.lon
    
    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] # there might be a better way of doing this!
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

    fig = plt.figure()
    ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)

    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    
    array, lons_cyclic = addcyclic(array, lons)
    array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])


    lon, lat = np.meshgrid(lons_cyclic, selected_lats)
    xi, yi = m(lon, lat)

    zonavg = array.mean(dim='lon')
    zonavg_thin = zonavg
    zonavg = np.expand_dims(zonavg,axis=1)
    zonavg = np.repeat(zonavg,lons_cyclic.shape[0],axis=1)
    array = array - zonavg # array to plot is now minus the zonal average
    
    minval = np.absolute(array.min())
    maxval = np.absolute(array.max())
    
    if maxval >= minval:
	    minval = - maxval
    else: 
	    maxval = minval
	    minval = - minval
    
    dlons = lons[100] - lons[99]
    dlats = lats[60] - lats[59]


    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',vmin=minval, vmax=maxval)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG,vmin=minval, vmax=maxval)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.seismic,vmin=minval, vmax=maxval)
    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))

    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)

    landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land_square.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    
    square_lons,square_lats=(xr.DataArray(np.meshgrid(lons, lats))).where(landmask==1.)
    square_lons = square_lons + 180. + dlons


    m.drawgreatcircle(square_lons.min(), square_lats.min(), square_lons.min(), square_lats.max()+dlats, del_s=100., color='b')
    m.drawgreatcircle(square_lons.max(), square_lats.min(), square_lons.max(), square_lats.max()+dlats, del_s=100., color='b')

    lats = [square_lats.min(),square_lats.min()] 
    lons = [square_lons.min(),square_lons.max()] 
    x, y = m(lons, lats) 
    m.plot(x, y, color='b') 

    lats = [square_lats.max()+dlats,square_lats.max()+dlats] 
    lons = [square_lons.min(),square_lons.max()] 
    x, y = m(lons, lats) 
    m.plot(x, y, color='b') 
    plt.title(title)

    ax2 = plt.subplot2grid((5,7), (0,6), rowspan = 3)
    
    plt.plot(zonavg_thin,array['lat'])
    plt.ylabel('Latitude')
    plt.xlabel(zonavgtitle+' ('+units+')')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.invert_xaxis()
    plt.show()
    return(square_lons,square_lats)

# DONT USE UNTIL FIGURED OUT WHAT IS WRONG THERE! 
# def aquaplanet_plot_for_iPython(minlat,maxlat,array,month,units,title,palette,contourson=False):
# for some reason this plots a sort of zonal average or something (eg. for bucket depth just stripes over land as well, should be zero there and is zer o if I just use xarray .plot() )
#     lats=np.linspace(-90.,90.,len(array.dim_2))
#     lons=np.linspace(0.,360.,len(array.dim_3))

#     if month == 12:
# 	    array = array.mean('dim_0')
#     else:
# 	    array = array[month,:,:]
#     print(np.shape(lons))
    
#     minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
#     maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] 
#     selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

#     fig = plt.figure()
#     ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)


#     m = Basemap(projection='kav7',lon_0=0.,resolution='c')

#     array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
#     array,lons_cyclic = addcyclic(array, lons)
#     array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

#     lon, lat = np.meshgrid(lons_cyclic, selected_lats)
#     xi, yi = m(lon, lat)

#     zonavg_thin = array.mean(dim='lon')
#     meravg_thin = array.mean(dim='lat')
#     m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
#     m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    
#     minval = np.absolute(array.min())
#     maxval = np.absolute(array.max())
    
#     if maxval >= minval:
# 	    minval = - maxval
#     else: 
# 	    maxval = minval
# 	    minval = - minval


#     if palette=='rainnorm':
#         cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG')
#     elif palette == 'raindefault':
#         cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG)
#     elif palette=='temp': 
#         cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='RdBu_r',vmin = 273.15-(maxval-273.15),
# 			      vmax=maxval)
#     elif palette=='fromwhite':
#        	    pal = plt.cm.Blues
#             pal.set_under('w',None)
#       	    cs = m.pcolormesh(xi,yi,array,cmap=pal,vmin=0,vmax=maxval)
#     else:
#         cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))


#     cbar = m.colorbar(cs, location='right', pad="10%")
#     cbar.set_label(units)

#     if contourson == True:  # add contours 
# 	    cont = m.contour(xi,yi,array,4,cmap='PuBu_r',
# 			     linewidth=5)
# 	    plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=12)

#     plt.title(title)

#     ax2 = plt.subplot2grid((5,7), (0,6), rowspan = 3)
#     # zonal average plot
#     plt.plot(zonavg_thin,array['lat'])
#     plt.ylabel('Latitude')
#     plt.xlabel(title+' ('+units+')')
#     ax2.yaxis.tick_right()
#     ax2.yaxis.set_label_position('right')
#     ax2.invert_xaxis()

#     ax3 = plt.subplot2grid((5,7), (4,0), colspan = 4)

# #    for producing a meridional average plot
#     plt.plot(array['lon'],meravg_thin)
#     plt.xlabel('Longitude')
#     plt.ylabel(title+' ('+units+')')
#     plt.show()

def aquaplanet_plot(minlat,maxlat,array,units,title,palette,contourson = False):

    lats = array.lat
    lons = array.lon

    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] 
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]


    fig = plt.figure()
    ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)

    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
    array,lons_cyclic = addcyclic(array, lons)
    array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

    lon, lat = np.meshgrid(lons_cyclic, selected_lats)
    xi, yi = m(lon, lat)

    zonavg_thin = array.mean(dim='lon')
    meravg_thin = array.mean(dim='lat')
    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])


    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG')
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='RdBu_r',vmin = 273.15-(maxval-273.15),
			      vmax=maxval)
    elif palette=='fromwhite':
       	    pal = plt.cm.Blues
            pal.set_under('w',None)
      	    cs = m.pcolormesh(xi,yi,field,cmap=pal,vmin=0,vmax=maxval)
    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))


    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)

    if contourson == True:  # add contours 
	    cont = m.contour(xi,yi,array,4,cmap='PuBu_r',
			     linewidth=5)
	    plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=12)

    plt.title(title)

    ax2 = plt.subplot2grid((5,7), (0,6), rowspan = 3)
    # zonal average plot
    plt.plot(zonavg_thin,array['lat'])
    plt.ylabel('Latitude')
    plt.xlabel(title+' ('+units+')')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.invert_xaxis()

    ax3 = plt.subplot2grid((5,7), (4,0), colspan = 4)

#    for producing a meridional average plot
    plt.plot(meravg_thin,array['lon'])
    plt.ylabel('Longitude')
    plt.xlabel(title+' ('+units+')')
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position('right')
    ax3.invert_xaxis()
    plt.show()


def aquaplanet_plot_minuszonavg(minlat,maxlat,array,units,title,palette):

    plt.close()

    lats=array.lat
    lons=array.lon
    print(np.shape(lons))
    
    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] 
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

    array, lons_cyclic = addcyclic(array, lons)
    array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
    lon, lat = np.meshgrid(lons_cyclic, selected_lats)
    xi, yi = m(lon, lat)

    zonavg = array.mean(dim='lon')
    zonavg = np.expand_dims(zonavg,axis=1)
    zonavg = np.repeat(zonavg,lons_cyclic.shape[0],axis=1)
    array = array - zonavg # array to plot is now minus the zonal average
    
    minval = np.absolute(array.min())
    maxval = np.absolute(array.max())
    
    if maxval >= minval:
	    minval = - maxval
    else: 
	    maxval = minval
	    minval = - minval


    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',vmin=minval, vmax=maxval)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG,vmin=minval, vmax=maxval)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.seismic,vmin=minval, vmax=maxval)
    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))

    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)
    plt.title(title)

    plt.show()





# def squareland_plot_several(minlat,maxlat,array1,title1,array2,title2,array3,title3,units,title,palette):

#     plt.close()

#     lats=array1.lat
#     lons=array1.lon


#     minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
#     maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0]
#     selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]



#     fig = plt.figure()
#     ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)
    
#     m = Basemap(projection='cyl',llcrnrlat=minlat,urcrnrlat=maxlat,\
#             llcrnrlon=-181.,urcrnrlon=180.,resolution='c')

#     m.drawparallels(np.arange(minlat,maxlat+1.,30.),labels=[1,0,0,0])
#     m.drawmeridians(np.arange(-180.,181.,60.),labels=[0,0,0,1])


#     lon, lat = np.meshgrid(lons, selected_lats)
#     xi, yi = m(lon, lat)

#     landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land_square.nc'),mode='r')
#     landmask=landfile.variables['land_mask'][:]
#     landlats=landfile.variables['lat'][:]
#     landlons=landfile.variables['lon'][:]
    

#     square_lons,square_lats=(xr.DataArray(np.meshgrid(lons, lats))).where(landmask==1.)

#     xs=[square_lons.min(),square_lons.max(),square_lons.max(),square_lons.min(),square_lons.min()]
#     ys=[square_lats.min(),square_lats.min(),square_lats.max()+dlats,square_lats.max()+dlats,square_lats.min()]



#     if palette=='rainnorm':
#         cs = m.pcolor(xi,yi,array1.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',latlon=True)
#     elif palette == 'raindefault':
#         cs = m.pcolor(xi,yi,array1.sel(lat=selected_lats), cmap=plt.cm.BrBG,latlon=True)
#     elif palette=='temp':
#         cs = m.pcolor(xi,yi,array1.sel(lat=selected_lats), cmap=plt.cm.seismic,latlon=True)
#     else:
#         cs = m.pcolor(xi,yi,array1.sel(lat=selected_lats),latlon=True)

#     m.plot(xs, ys, latlon = True)

#     axes[1].set_title(title2)

#     if palette=='rainnorm':
#         cs = m.pcolor(xi,yi,array2.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',latlon=True)
#     elif palette == 'raindefault':
#         cs = m.pcolor(xi,yi,array2.sel(lat=selected_lats), cmap=plt.cm.BrBG,latlon=True)
#     elif palette=='temp': 
#         cs = m.pcolor(xi,yi,array2.sel(lat=selected_lats), cmap=plt.cm.seismic,latlon=True)
#     else:
#         cs = m.pcolor(xi,yi,array2.sel(lat=selected_lats),latlon=True)

#     m.plot(xs, ys, latlon = True)

#     axes[2].set_title(title3)

#     m = Basemap(projection='cyl',llcrnrlat=minlat,urcrnrlat=maxlat,\
#             llcrnrlon=-181.,urcrnrlon=180.,resolution='c')
#     m.drawparallels(np.arange(minlat,maxlat+1.,30.),labels=[1,0,0,0])
#     m.drawmeridians(np.arange(-180.,181.,60.),labels=[0,0,0,1])


#     lon, lat = np.meshgrid(lons, selected_lats)
#     xi, yi = m(lon, lat)

#     landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land_square.nc'),mode='r')
#     landmask=landfile.variables['land_mask'][:]
#     landlats=landfile.variables['lat'][:]
#     landlons=landfile.variables['lon'][:]

#     square_lons,square_lats=(xr.DataArray(np.meshgrid(lons, lats))).where(landmask==1.)

#     xs=[square_lons.min(),square_lons.max(),square_lons.max(),square_lons.min(),square_lons.min()]
#     ys=[square_lats.min(),square_lats.min(),square_lats.max()+dlats,square_lats.max()+dlats,square_lats.min()]


#     if palette=='rainnorm':
#         cs = m.pcolor(xi,yi,array3.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',latlon=True)
#     elif palette == 'raindefault':
#         cs = m.pcolor(xi,yi,array3.sel(lat=selected_lats), cmap=plt.cm.BrBG,latlon=True)
#     elif palette=='temp': 
#         cs = m.pcolor(xi,yi,array3.sel(lat=selected_lats), cmap=plt.cm.seismic,latlon=True)
#     else:
#         cs = m.pcolor(xi,yi,array3.sel(lat=selected_lats),latlon=True)

#     m.plot(xs, ys, latlon = True)


# # Add Title
#     plt.title(title)


# # Add Colorbar
#     cbar = plt.colorbar(cs, location='right', pad="10%")
#     cbar.set_label(units)

#     plt.show()

#     return(square_lons,square_lats)




def worldmap_inputfile(filename,varname): # assuming that the input file has the dimenstions time, lat, lon. Plots annual average

    nc=Dataset(filename,mode='r')

    array=xr.DataArray(nc.variables[varname][:])
    lats=xr.DataArray(nc.variables['lat'][:])
    lons=xr.DataArray(nc.variables['lon'][:])

    lon_0 = lons.mean() # trick to find central lon lat (lon from 0 to 360)
    lat_0 = lats.mean() # ~= 0 because lat from -90 to 90


    fig = plt.figure()
    ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)


    m = Basemap(projection='kav7',lon_0=180.,resolution='c')

    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    array = array.mean('dim_0')
    array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

    zonavg_thin = array.mean(dim='lon')
    meravg_thin = array.mean(dim='lat')
    
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    cs = m.pcolor(xi,yi,array)
    cbar = m.colorbar(cs, location='right', pad="10%")
    

    ax2 = plt.subplot2grid((5,7), (0,6), rowspan = 3)
    # zonal average plot
    plt.plot(zonavg_thin,array['lat'])
    plt.ylabel('Latitude')
    plt.xlabel(varname)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.invert_xaxis()

    ax3 = plt.subplot2grid((5,7), (4,0), colspan = 4)
    # meridional average plot
    plt.plot(meravg_thin,array['lon'])
    plt.ylabel('Longitude')
    plt.xlabel(varname)
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position('right')
    ax3.invert_xaxis()
    plt.show()



def worldmap_variable(outdir,field,units,title,palette,minval,maxval,nmb_contours=0.):
    plt.close()
    

    small = 12 #largefonts 14 # smallfonts 12 # medfonts = 14
    med = 18 #largefonts 18 # smallfonts 14 # medfonts = 16
    lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20


    lats=field.lat
    lons=field.lon
    lon_0 = lons.mean() 
    lat_0 = lats.mean() 

    fig = plt.figure(figsize = (35,15))

    ax1 = plt.subplot2grid((2,2), (0,0))


    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    field, lons = addcyclic(field, lons)
    field = np.asarray(field) 
    field,lons = shiftgrid(np.max(lons)-180.,field,lons,start=False,cyclic=np.max(lons))
    field = xr.DataArray(field,coords=[lats,lons],dims=['lat','lon'])


    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(0.,361.,60.))
    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=small)

    if (minval==0) & (maxval==0):

	    if palette=='rainnorm':
		    cs = m.pcolor(xi,yi,field,norm=MidpointNormalize(midpoint=0.),cmap='BrBG')
	    elif palette == 'raindefault':
		    cs = m.pcolor(xi,yi,field, cmap=plt.cm.BrBG)
	    elif palette=='temp': 
		    cs = m.pcolor(xi,yi,field, cmap=plt.cm.seismic)
	    else:
		    cs = m.pcolor(xi,yi,field)
    else:

	    if palette=='rainnorm':
		    cs = m.pcolor(xi,yi,field,norm=MidpointNormalize(midpoint=0.),cmap='BrBG', vmin=minval, vmax=maxval)
	    elif palette == 'raindefault':
		    cs = m.pcolor(xi,yi,field, cmap=plt.cm.BrBG, vmin=minval, vmax=maxval)
	    elif palette=='temp':
		    cs = m.pcolor(xi,yi,field,norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),vmax=maxval)




	    elif palette=='fromwhite':
		    pal = plt.cm.Blues
		    pal.set_under('w',None)
		    cs = m.pcolormesh(xi,yi,field,cmap=pal,vmin=0,vmax=maxval)

	    else:
		    cs = m.pcolor(xi,yi,field, vmin=minval, vmax=maxval)



    if nmb_contours != 0:  # add contours 
	    cont = m.contour(xi,yi,field,nmb_contours,cmap='PuBu_r', linewidth=5)
	    if cont>=1.:
		    plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=med)
	    else:
		    plt.clabel(cont, inline=2, fmt='%1.3f', fontsize=med)

    cbar = m.colorbar(cs, location='bottom', pad="10%")
    cbar.set_label(units, size=med)
    cbar.ax.tick_params(labelsize=small) 
    plt.title(title, size=lge)
    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+title+'.png', bbox_inches='tight', dpi=100)


def squareland_inputfile_4dvar(filename,varname):
	
	nc=Dataset(filename,mode='r')

	array=xr.DataArray(nc.variables[varname][:])
	array = array.sum(dim='dim=1').mean(dim='dim_0')
	lats=xr.DataArray(nc.variables['lat'][:])
	lons=xr.DataArray(nc.variables['lon'][:])

	fig = plt.figure()
	ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)


	m = Basemap(projection='kav7',lon_0=0.,resolution='c')

	array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
	array, lons_cyclic = addcyclic(array, lons)
	array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
	zonavg_thin = array.mean(dim='lon')


	lon, lat = np.meshgrid(lons_cyclic, lats)
	xi, yi = m(lon, lat)

	dlons = lons[100] - lons[99]
	dlats = lats[60] - lats[59]

      	cs = m.pcolor(xi,yi,array)
	cbar = m.colorbar(cs, location='right', pad="10%")

# Show landmask
	landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land_square.nc'),mode='r')
	landmask=landfile.variables['land_mask'][:]
	landlats=landfile.variables['lat'][:]
	landlons=landfile.variables['lon'][:]
	
	square_lons,square_lats=(xr.DataArray(np.meshgrid(landlons, landlats))).where(landmask==1.)
	square_lons = square_lons + dlons

	m.drawgreatcircle(square_lons.min(), square_lats.min(), square_lons.min(), square_lats.max()+dlats, del_s=100., color='b')
	m.drawgreatcircle(square_lons.max(), square_lats.min(), square_lons.max(), square_lats.max()+dlats, del_s=100., color='b')

	lats = [square_lats.min(),square_lats.min()] 
	lons = [square_lons.min(),square_lons.max()] 
	x, y = m(lons, lats) 
	m.plot(x, y, color='b') 

	lats = [square_lats.max()+dlats,square_lats.max()+dlats] 
	lons = [square_lons.min(),square_lons.max()] 
	x, y = m(lons, lats) 
	m.plot(x, y, color='b') 
	plt.title(varname)	

	ax2 = plt.subplot2grid((5,7), (0,6), rowspan = 3)
	
	plt.plot(zonavg_thin,array['lat'])
	plt.ylabel('Latitude')
	plt.xlabel(varname)
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position('right')
	ax2.invert_xaxis()
	plt.show()

def squareland_inputfile_3dvar(filename,varname):
	
	nc=Dataset(filename,mode='r')

	array=xr.DataArray(nc.variables[varname][:])
	array = array.mean(dim='dim_0')
	lats=xr.DataArray(nc.variables['lat'][:])
	lons=xr.DataArray(nc.variables['lon'][:])

	fig = plt.figure()
	ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)


	m = Basemap(projection='kav7',lon_0=0.,resolution='c')

	array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
	array, lons_cyclic = addcyclic(array, lons)
	array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
	zonavg_thin = array.mean(dim='lon')
	meravg_thin = array.mean(dim='lat')

	lon, lat = np.meshgrid(lons_cyclic, lats)
	xi, yi = m(lon, lat)

	dlons = lons[100] - lons[99]
	dlats = lats[60] - lats[59]

	cs = m.pcolor(xi,yi,array,norm=MidpointNormalize(midpoint=0.),cmap='seismic')
	cbar = m.colorbar(cs, location='right', pad="10%")

# Show landmask
	landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land.nc'),mode='r')
	landmask=landfile.variables['land_mask'][:]
	landlats=landfile.variables['lat'][:]
	landlons=landfile.variables['lon'][:]
	
	square_lons,square_lats=(xr.DataArray(np.meshgrid(landlons, landlats))).where(landmask==1.)
	square_lons = square_lons + dlons

	m.drawgreatcircle(square_lons.min(), square_lats.min(), square_lons.min(), square_lats.max()+dlats, del_s=100., color='b')
	m.drawgreatcircle(square_lons.max(), square_lats.min(), square_lons.max(), square_lats.max()+dlats, del_s=100., color='b')

	lats = [square_lats.min(),square_lats.min()] 
	lons = [square_lons.min(),square_lons.max()] 
	x, y = m(lons, lats) 
	m.plot(x, y, color='b') 

	lats = [square_lats.max()+dlats,square_lats.max()+dlats] 
	lons = [square_lons.min(),square_lons.max()] 
	x, y = m(lons, lats) 
	m.plot(x, y, color='b') 
	plt.title(varname)	

	ax2 = plt.subplot2grid((5,7), (0,6), rowspan = 3)
	
	plt.plot(zonavg_thin,array['lat'])
	plt.ylabel('Latitude')
	plt.xlabel(varname)
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position('right')
	ax2.invert_xaxis()

	ax3 = plt.subplot2grid((5,7), (4,0), colspan = 4)
#    for producing a meridional average plot
	plt.plot(meravg_thin,array['lon'])
	plt.ylabel('Longitude')
	plt.xlabel(varname)
	ax3.yaxis.tick_right()
	ax3.yaxis.set_label_position('right')
	ax3.invert_xaxis()
	plt.show()

def squareland_inputfile_2dvar(filename,varname):
	
	nc=Dataset(filename,mode='r')

	array=xr.DataArray(nc.variables[varname][:])
	lats=xr.DataArray(nc.variables['lat'][:])
	lons=xr.DataArray(nc.variables['lon'][:])

	fig = plt.figure()
	ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)


	m = Basemap(projection='kav7',lon_0=0.,resolution='c')

	array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
	array, lons_cyclic = addcyclic(array, lons)
	array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
	zonavg_thin = array.mean(dim='lon')
	meravg_thin = array.mean(dim='lat')

	lon, lat = np.meshgrid(lons_cyclic, lats)
	xi, yi = m(lon, lat)

	dlons = lons[100] - lons[99]
	dlats = lats[60] - lats[59]

	cs = m.pcolor(xi,yi,array,norm=MidpointNormalize(midpoint=0.),cmap='seismic')
	cbar = m.colorbar(cs, location='right', pad="10%")

# Show landmask
	landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land_square.nc'),mode='r')
	landmask=landfile.variables['land_mask'][:]
	landlats=landfile.variables['lat'][:]
	landlons=landfile.variables['lon'][:]
	
	square_lons,square_lats=(xr.DataArray(np.meshgrid(landlons, landlats))).where(landmask==1.)
	square_lons = square_lons + dlons

	m.drawgreatcircle(square_lons.min(), square_lats.min(), square_lons.min(), square_lats.max()+dlats, del_s=100., color='b')
	m.drawgreatcircle(square_lons.max(), square_lats.min(), square_lons.max(), square_lats.max()+dlats, del_s=100., color='b')

	lats = [square_lats.min(),square_lats.min()] 
	lons = [square_lons.min(),square_lons.max()] 
	x, y = m(lons, lats) 
	m.plot(x, y, color='b') 

	lats = [square_lats.max()+dlats,square_lats.max()+dlats] 
	lons = [square_lons.min(),square_lons.max()] 
	x, y = m(lons, lats) 
	m.plot(x, y, color='b') 
	plt.title(varname)	

	ax2 = plt.subplot2grid((5,7), (0,6), rowspan = 3)
	
	plt.plot(zonavg_thin,array['lat'])
	plt.ylabel('Latitude')
	plt.xlabel(varname)
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position('right')
	ax2.invert_xaxis()

	ax3 = plt.subplot2grid((5,7), (4,0), colspan = 4)
#    for producing a meridional average plot
	plt.plot(meravg_thin,array['lon'])
	plt.ylabel('Longitude')
	plt.xlabel(varname)
	ax3.yaxis.tick_right()
	ax3.yaxis.set_label_position('right')
	ax3.invert_xaxis()
	plt.show()



def aquaplanet_inputfile(filename,varname,month,palette=None):
	
	nc=Dataset(filename,mode='r')

	

 	if (month == 12): # want annual mean
		array=xr.DataArray(nc.variables[varname][:])
		array = array.mean(dim='dim_0')
	elif (month < 0): #no time dimension
		array=xr.DataArray(nc.variables[varname][:])
	else: # want specific month 0-11
		array=xr.DataArray(nc.variables[varname][month,:,:])


	lats=xr.DataArray(nc.variables['lat'][:])
	lons=xr.DataArray(nc.variables['lon'][:])

	fig = plt.figure()
	ax1 = plt.subplot2grid((5,7), (0,0), colspan = 5, rowspan = 3)


	m = Basemap(projection='kav7',lon_0=0.,resolution='c')

	array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
	array, lons_cyclic = addcyclic(array, lons)
	array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
	zonavg_thin = array.mean(dim='lon')
	meravg_thin = array.mean(dim='lat')


	lon, lat = np.meshgrid(lons_cyclic, lats)
	xi, yi = m(lon, lat)

	minval = np.absolute(array.min())
	maxval = np.absolute(array.max())
	
	if maxval >= minval:
		minval = - maxval
	else: 
		maxval = minval
		minval = - minval

	dlons = lons[100] - lons[99]
	dlats = lats[60] - lats[59]

	if palette=='rainnorm':
		cs = m.pcolor(xi,yi,array,norm=MidpointNormalize(midpoint=0.),cmap='BrBG',vmin=minval, vmax=maxval)
	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, cmap=plt.cm.BrBG)
	elif palette=='temp':
		cs = m.pcolor(xi,yi,array,norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r,vmin = 270.,vmax=maxval)

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,cmap=pal,vmin=0,vmax=maxval)
	elif palette=='tempdiff': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=0), cmap=plt.cm.RdBu_r, 
			      vmin = -maxval,
			      vmax = maxval)
	else:
		cs = m.pcolor(xi,yi,array)

	cbar = m.colorbar(cs, location='right', pad="10%")

	plt.title(varname)     

	plt.show()

def squareland_plot_correlation(minlat,maxlat,array1,array2,title):

    plt.close()

    lats=array1.lat
    lons=array1.lon
    
    array=xr.DataArray(signal.correlate2d(array1,array2,boundary='symm',mode='same'))
    array=xr.DataArray(array.values,coords=[lats,lons],dims=['lat','lon'])

    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] # there might be a better way of doing this!
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

    # sqlatindex=np.asarray(np.where(lats>sq_minlat))[0,0]
    # sqlatreverseindex=np.asarray(np.where(lats[::-1]<sq_maxlat))[0,0] # there might be a better way of doing this!
    # square_lats=lats[sqlatindex:(lats.size-sqlatreverseindex)]

    # sqlonindex=np.asarray(np.where(lons>sq_minlon))[0,0]
    # sqlonreverseindex=np.asarray(np.where(lons[::-1]<sq_maxlon))[0,0] # there might be a better way of doing this!
    # square_lons=lons[sqlonindex:(lons.size-sqlonreverseindex)+1]





    m = Basemap(projection='cyl',llcrnrlat=minlat,urcrnrlat=maxlat,\
            llcrnrlon=-181.,urcrnrlon=180.,resolution='c')
#    m.drawcoastlines()
#m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
    m.drawparallels(np.arange(minlat,maxlat+1.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,181.,60.),labels=[0,0,0,1])


# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
    lon, lat = np.meshgrid(lons, selected_lats)
    xi, yi = m(lon, lat)

#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),latlon=True,cmap='gray')

# Add Colorbar
    cbar = m.colorbar(cs, location='right', pad="10%")

    # sns.palplot(sns.color_palette("BrBG", 7))

# Show landmask
    landfile=Dataset(os.path.join(GFDL_BASE,'input/squareland/land_square.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    
#    m.contour(xi, yi, landmask, colors='k')

    square_lons,square_lats=(xr.DataArray(np.meshgrid(lons, lats))).where(landmask==1.)

# Add rectangle
    # xs = [sq_minlon,sq_maxlon,sq_maxlon,sq_minlon,sq_minlon]
    # ys = [sq_minlat,sq_minlat,sq_maxlat,sq_maxlat,sq_minlat]
    xs=[square_lons.min(),square_lons.max(),square_lons.max(),square_lons.min(),square_lons.min()]
    ys=[square_lats.min(),square_lats.min(),square_lats.max()+dlats,square_lats.max()+dlats,square_lats.min()]
    m.plot(xs, ys, latlon = True)


# Add Title
    plt.title(title)


    plt.show()

def any_configuration_plot(outdir,runmin,runmax,minlat,maxlat,array,area_array,units,title,palette,landmaskxr,nmb_contours=0,minval=None,maxval=None,month_annotate=None,save_fig=True):
# plotting only the zonal average next to the map 
# currently hard coded -30.,30. slice instead of squarelats_min, squarelats_max
    plt.close()


    small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
    med = 18 #largefonts 18 # smallfonts 14 # medfonts = 16
    lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

    lats=array.lat
    lons=array.lon
    
    # why is this not working anymore when land areas are selected ? worked in commit d110990e
    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] 
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

    landlats = np.asarray(landmaskxr.lat)
    landlons = np.asarray(landmaskxr.lon)

    landmask = np.asarray(landmaskxr)


    fig = plt.figure(figsize = (25,10))

    ax1 = plt.subplot2grid((5,8), (0,1), colspan = 5, rowspan = 3)


    m = Basemap(projection='kav7',lon_0=0.,llcrnrlon=-180.,llcrnrlat=-30.,urcrnrlon=180.,urcrnrlat=30.,resolution='c')
#    m = Basemap(projection='cyl',llcrnrlon=-180.,llcrnrlat=-30.,urcrnrlon=180.,urcrnrlat=30.,resolution='c') works, but not with kav7 projection
    array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

    zonavg_thin = area_weighted_avg(array,area_array,landmaskxr,option = 'all_sfcs',minlat=-90.,maxlat=90.,axis=1)
    meravg_thin = area_weighted_avg(array,area_array,landmaskxr,option = 'all_sfcs',minlat=-30.,maxlat=30.,axis=0)

    lons_128 = lons # non-cyclic lons, i.e. lenght = 128
    # #newline to replace shiftgrid line which is causing trouble - Doesn't work perfectly
    # lons, array = m.shiftdata(lons, datain = array, lon_0=0.)
    # array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

    array = np.asarray(array) #- This line fixes the problem!
    #the following line caused DataType error all of a sudden... Doesnt' accept xarray as input array for shiftgrid anymore.

    array, lons = addcyclic(array, lons)
    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# shiftgrid lons0 (first entry) is the longitude in lons which should correspond with the center longitude of the map. start = False --> lonsin is the end latitude, not the beginning.
    # this doesn't work for some reason
    #array, lons = shiftgrid(np.max(lons)-100.,array,lons,start=True,cyclic=np.max(lons))

    array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=small)

    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    dlons = lons[100] - lons[99]
    dlats = lats[60] - lats[59]


    if minval==None and maxval==None:
	    minval = array.min()
	    maxval = array.max()
    
    minval = np.absolute(minval)
    maxval = np.absolute(maxval)

    if maxval >= minval:
	    minval = - maxval
    else: 
	    maxval = minval
	    minval = - minval


    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',vmin=minval, vmax=maxval)
    elif palette == 'PE_scale':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='bwr_r',vmin=minval, vmax=maxval)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG)
    elif palette=='temp':
	    cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),vmax=maxval)
#	    cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r,vmin = 273.15-30.,vmax=273.15+30.) # forpaper

    elif palette=='fromwhite': 
	    pal = plt.cm.Blues
	    pal.set_under('w',None)
	    cs = m.pcolormesh(xi,yi,array.sel(lat=selected_lats),cmap=pal,vmin=0,vmax=maxval)

    elif palette=='bucket': 
	    pal = plt.cm.Greens
	    pal.set_under('w',None)
	    cs = m.pcolormesh(xi,yi,array.sel(lat=selected_lats),cmap=pal,vmin=0,vmax=maxval)

    elif palette=='tempdiff': 
	    cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), 
		      norm=MidpointNormalize(midpoint=0), cmap=plt.cm.RdBu_r, 
		      vmin = -maxval,
		      vmax = maxval)
    elif palette=='slp':
	    cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.coolwarm,vmin=minval,vmax=maxval)



    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))


    if nmb_contours != 0:  # add contours 
	    cont = m.contour(xi,yi,array,nmb_contours,cmap='PuBu_r', linewidth=5)
	    if cont>=1.:
		    plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=med)
	    else:
		    plt.clabel(cont, inline=2, fmt='%1.3f', fontsize=med)



# Add Colorbar
    cbar = m.colorbar(cs, location='bottom', pad="10%") # usually on right 
    cbar.set_label(units, size=med)
    cbar.ax.tick_params(labelsize=small) 

    # sns.palplot(sns.color_palette("BrBG", 7))

# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
    landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

    landmask, lons_cyclic = addcyclic(landmask, landlons)

    if np.any(landmask != 0.):
	    m.contour(xi,yi,landmask, 1)

    plt.title(title, size=lge)

    
    if month_annotate >= 1:
	    plt.annotate('Month #'+str(month_annotate), xy=(0.15,0.8), xycoords='figure fraction')
	    return fig

    else:
   
	    ax2 = plt.subplot2grid((5,8), (0,6), rowspan = 3)

	    plt.plot(zonavg_thin,lats)
	    plt.ylabel('Latitude', size=med)
	    plt.xlabel(title+' ('+units+')', size=med)
	    ax2.yaxis.tick_right()
	    ax2.yaxis.set_label_position('right')
	    ax2.tick_params(axis='both', which='major', labelsize=small)
	    ax2.invert_xaxis()

# 	    ax3 = plt.subplot2grid((5,8), (4,1), colspan = 4)
# 	    plt.plot(lons_128,meravg_thin)
# 	    plt.xlabel('Longitude')
# 	    plt.ylabel(title+' ('+units+') 30S-30N')
# #	    plt.tight_layout()

	    manager = plt.get_current_fig_manager()
	    manager.window.showMaximized()
#	    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+title+'_'+str(runmin)+'-'+str(runmax)+'_highres.png', format = 'png', dpi = 400, bbox_inches='tight')

	    if save_fig == True:
		    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+title+'_'+str(runmin)+'-'+str(runmax)+'.png', format = 'png', bbox_inches='tight')
#	    plt.show()
    return fig


def any_configuration_plot_allmonths(outdir,runmin,runmax,minlat,maxlat,array,area_array,units,title,palette,landmaskxr,nmb_contours=0,minval=None,maxval=None,month_annotate=None,save_fig=True):
# plotting only the zonal average next to the map 
# currently hard coded -30.,30. slice instead of squarelats_min, squarelats_max
    plt.close()


    small = 10 #largefonts 14 # smallfonts 10 # medfonts = 14
    med = 14 #largefonts 18 # smallfonts 14 # medfonts = 16
    lge = 18 #largefonts 22 # smallfonts 18 # medfonts = 20

    time = array.month
    lats=array.lat
    lons=array.lon
    
    # why is this not working anymore when land areas are selected ? worked in commit d110990e
    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] 
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

    landlats = np.asarray(landmaskxr.lat)
    landlons = np.asarray(landmaskxr.lon)

    landmask = np.asarray(landmaskxr)


    fig = plt.figure(figsize = (25,10))

    array = xr.DataArray(array,coords=[time,lats,lons],dims=['time','lat','lon'])
    lons_128 = lons # non-cyclic lons, i.e. lenght = 128
    array = np.asarray(array) #- This line fixes the problem!
    array, lons = addcyclic(array, lons)
    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
    array = xr.DataArray(array,coords=[time,lats,lons],dims=['time','lat','lon'])


    m = Basemap(projection='kav7',lon_0=0.,llcrnrlon=-180.,llcrnrlat=-30.,urcrnrlon=180.,urcrnrlat=30.,resolution='c')

    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=small)
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
    landmask, lons_cyclic = addcyclic(landmask, landlons)

    months = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

    for i in range (0,12):
	    if (i < 4):
		    ax = plt.subplot2grid((3,4), (0,i))
	    elif (i >= 4) and (i < 8):
		    ax = plt.subplot2grid((3,4), (1,i-4))
	    elif (i >= 8):
		    ax = plt.subplot2grid((3,4), (2,i-8))

	    if minval==None and maxval==None:
		    minval = array[i,:,:].min()
		    maxval = array[i,:,:].max()

	    minval = np.absolute(minval)
	    maxval = np.absolute(maxval)

	    if maxval >= minval:
		    minval = - maxval
	    else: 
		    maxval = minval
		    minval = - minval


	    if palette=='rainnorm':
		cs = m.pcolor(xi,yi,array[i,:,:].sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',vmin=minval, vmax=maxval)
	    elif palette == 'PE_scale':
		cs = m.pcolor(xi,yi,array[i,:,:].sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='bwr_r',vmin=minval, vmax=maxval)
	    elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array[i,:,:].sel(lat=selected_lats), cmap=plt.cm.BrBG)
	    elif palette=='temp':
		   # cs = m.pcolor(xi,yi,array[i,:,:].sel(lat=selected_lats),norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),vmax=maxval)
		    cs = m.pcolor(xi,yi,array[i,:,:].sel(lat=selected_lats),norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r,vmin = 273.15-30.,vmax=273.15+30.) # forpaper

	    elif palette=='fromwhite': 
		    pal = plt.cm.Blues
		    pal.set_under('w',None)
		    cs = m.pcolormesh(xi,yi,array[i,:,:].sel(lat=selected_lats),cmap=pal,vmin=0,vmax=maxval)

	    elif palette=='bucket': 
		    pal = plt.cm.Greens
		    pal.set_under('w',None)
		    cs = m.pcolormesh(xi,yi,array[i,:,:].sel(lat=selected_lats),cmap=pal,vmin=0,vmax=maxval)

	    elif palette=='tempdiff': 
		    cs = m.pcolor(xi,yi,array[i,:,:].sel(lat=selected_lats), 
			      norm=MidpointNormalize(midpoint=0), cmap=plt.cm.RdBu_r, 
			      vmin = -maxval,
			      vmax = maxval)
	    elif palette=='slp':
		    cs = m.pcolor(xi,yi,array[i,:,:].sel(lat=selected_lats), cmap=plt.cm.coolwarm,vmin=minval,vmax=maxval)



	    else:
		cs = m.pcolor(xi,yi,array[i,:,:].sel(lat=selected_lats))


	    if nmb_contours != 0:  # add contours 
		    cont = m.contour(xi,yi,array[i,:,:],nmb_contours,cmap='PuBu_r', linewidth=5)
		    if cont>=1.:
			    plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=med)
		    else:
			    plt.clabel(cont, inline=2, fmt='%1.3f', fontsize=med)

	# Add rectangles

	    if np.any(landmask != 0.):
		    m.contour(xi,yi,landmask, 1)
	    
	    ax.annotate(months[i], xy=(0.,1.), xycoords='axes fraction')


	# Add Colorbar
    cbar = m.colorbar(cs, location='bottom', pad="10%") # usually on right 
    cbar.set_label(units, size=med)
    cbar.ax.tick_params(labelsize=small) 

	    # sns.palplot(sns.color_palette("BrBG", 7))

	# Read landmask

    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()

    if save_fig == True:
	    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+title+'_'+str(runmin)+'-'+str(runmax)+'.png', format = 'png', bbox_inches='tight')
    #	    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+title+'_'+str(runmin)+'-'+str(runmax)+'_highres.png', format = 'png', dpi = 400, bbox_inches='tight')

    plt.show()
    return fig











def animated_map(testdir,outdir,array,units,title,plot_title,palette,imin,imax,minval=0,maxval=1000,landlats = None, landlons = None, landmask = 0.):

# Can be used for a surface variable, e.g. to animate the climatology of evaporation

    plt.close()
    lats=array.lat
    lons=array.lon

    fig = plt.figure(figsize=(25,10))
    m = Basemap(projection='kav7',lon_0=0.,resolution='c')
    array = np.asarray(array)
    array, lons = addcyclic(array, lons)
    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# draw parallels and meridians.

    #m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    #m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    array = xr.DataArray(array)


# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.

    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)



    # print minval, maxval

    for idx in range (imin,imax): # exclues last value
	    if palette == 'rainnorm':

		    # minval = np.absolute(array.min())
		    # maxval = np.absolute(array.max())
		    
    
		    # if maxval >= minval:
		    # 	    minval = - maxval
		    # else: 
		    # 	    maxval = minval
		    # 	    minval = - minval
		    cs = m.pcolormesh(xi,yi,array[idx,:,:],cmap='BrBG',
				      norm=MidpointNormalize(midpoint=0.),
				      vmin=minval,vmax=maxval)
	    elif palette == 'raindefault':
		    cs = m.pcolor(xi,yi,array[idx,:,:],cmap=plt.cm.BrBG)
	    
	    elif palette=='temp':
		    cs = m.pcolor(xi,yi,array[idx,:,:],norm=MidpointNormalize(midpoint=273.15), cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),vmax=maxval)
	    elif palette=='fromwhite': 
		    pal = plt.cm.Blues
		    pal.set_under('w',None)
		    cs = m.pcolormesh(xi,yi,array[idx,:,:],cmap=pal,vmin=0,vmax=maxval)
	    elif palette=='slp':
		    minval = np.min(array)
		    maxval = np.max(array)
		    cs = m.pcolor(xi,yi,array[idx,:,:], cmap=plt.cm.coolwarm,vmin=minval,vmax=maxval)

	    cbar = m.colorbar(cs, location='bottom', pad="10%")
	    cbar.set_label(units)





# plot continent contours, not working 
	    # if np.any(landmask != 0.):

	    # 	    landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
	    # 	    landmask, lons_cyclic = addcyclic(landmask, landlons)
	    # 	    m.contour(xi,yi,landmask, 1)



	    plt.title(title)
	    plt.annotate('Month #'+str(idx+1), xy=(0.15,0.8), xycoords='figure fraction')
	    plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+plot_title+'_month'+str(1000+idx)+'.png',bbox_inches='tight')
	    fig.clf()
   

#+testdir+'/'
    os.system('convert -delay 50 /scratch/mp586/Code/Graphics/'+outdir+'/'+plot_title+'_month*.png /scratch/mp586/Code/Graphics/'+testdir+'/'+plot_title+'.gif')
    # Use ffmeg to convert animated gif to mp4
    # os.system('ffmpeg -f gif -i /scratch/mp586/Code/Graphics/'+testdir+'/'+plot_title+'.gif /scratch/mp586/Code/Graphics/'+testdir+'/'+plot_title+'.mp4')

	
def winds_at_heightlevel(uwind,vwind,level,array,palette,units,minval,maxval,landmaskxr,landlats,landlons,veclen=10):

# Plots every third wind vector at specified height level
# onto a world map of the 'array' which could be e.g. precip

# uwind and vwind are 4D arrays, 
# level should be between 0 and 39 for MiMA
# array is the underlying plot (e.g. Tsurf, precip, ...)
# palette is for the underlying plot
# units are for the underlying plot
	landmask = np.asarray(landmaskxr)

	fig = plt.figure(figsize = (25,10))
	m = Basemap(projection='kav7',lon_0=0.,resolution='c')
	lons = uwind.lon
	lats = uwind.lat
	pres = uwind.pres_lev[level]	
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	plt.title('Wind at '+str(pres)+' hPa')

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
	landmask, landlons = addcyclic(landmask, landlons)

	lon, lat = np.meshgrid(lons_shift, lats)
	xi, yi = m(lon, lat)

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)

	cbar = m.colorbar(cs, location='right', pad="10%")
	cbar.set_label(units)

	Q = plt.quiver(xi[::3,::3], yi[::3,::3], uwind[level,::3,::3], vwind[level,::3,::3], units='width')
	qk = plt.quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{m}{s}$', 
			   labelpos='E', coordinates='figure')

	# plt.show()
	return fig


# def winds_at_4_levels(): 
# fig, axes = plt.subplots(2,2,sharex = True, sharey = True, figsize = (25,10))
# axes[0,0].set_title("Level 39")



	
def winds_one_level(outdir,runmin,runmax,plt_title,uwind,vwind,array,palette,units,minval,maxval,landmaskxr,veclen=10,level=39, units_numerator = 'm', units_denom = 's',save = False):

# Plots every third wind vector at specified height level
# onto a world map of the 'array' which could be e.g. precip

# uwind and vwind are 4D arrays, 
# level should be between 0 and 39 for MiMA
# array is the underlying plot (e.g. Tsurf, precip, ...)
# palette is for the underlying plot
# units are for the underlying plot

	small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
	med = 18 #largefonts 18 # smallfonts 14 # medfonts = 16
	lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

	landlats = np.asarray(landmaskxr.lat)
	landlons = np.asarray(landmaskxr.lon)
	landmask = np.asarray(landmaskxr)

	fig, ax = plt.subplots(figsize = (25,10))

# if fig, ax = plt.subplots(0,0, figsize = (25,10)) then ax is a numpy array --> can't do quiver plot on it

	m = Basemap(projection='kav7',lon_0=0.,resolution='c')
	lons = uwind.lon
	lats = uwind.lat
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=med)
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=med)


	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
	landmask, landlons = addcyclic(landmask, landlons)

	lon, lat = np.meshgrid(lons_shift, lats)
	xi, yi = m(lon, lat)

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)
	elif palette == 'PE_scale':
		cs = m.pcolor(xi,yi,array,norm=MidpointNormalize(midpoint=0.),cmap='bwr_r',vmin=minval, vmax=maxval)
	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)

	cbar = m.colorbar(cs, location='right', pad="10%")
	cbar.set_label(units)

	Q = ax.quiver(xi[::3,::3], yi[::3,::3], uwind[::3,::3], vwind[::3,::3], units='width')
	ax.quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure')

	if save == True:
		fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/'+plt_title+'_level'+str(level)+'_'+str(runmin)+'-'+str(runmax)+'.png', dpi=100) # if add bbox_inches = 'tight' the quiverkey is not saved!

	return fig

def winds_seasons_one_level(uwind_in,vwind_in,level,array_in,palette,units,minval,maxval,landmaskxr,outdir,runmin,runmax,units_numerator='m',units_denom='s',quivkey=4,veclen=1.): 

	landlats = np.asarray(landmaskxr.lat)
	landlons = np.asarray(landmaskxr.lon)
	
	fig, axes = plt.subplots(2,2, figsize = (25, 10))

	axes[0,0].set_title('MAM')
	uwind = uwind_in.sel(season='MAM')
	vwind = vwind_in.sel(season='MAM')
	array = array_in.sel(season='MAM')

	landmask = np.asarray(landmaskxr)

#	fig = plt.figure()
	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,0])
	lons = uwind.lon
	lats = uwind.lat
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
	landmask, landlons = addcyclic(landmask, landlons)

	lon, lat = np.meshgrid(lons_shift, lats)
	xi, yi = m(lon, lat)

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)


	Q = axes[0,0].quiver(xi[::quivkey,::quivkey], yi[::quivkey,::quivkey], uwind[::quivkey,::quivkey], vwind[::quivkey,::quivkey], units='width')
	qk = axes[0,0].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', 
			   labelpos='E', coordinates='figure')


	axes[0,1].set_title('JJA')
	uwind = uwind_in.sel(season='JJA')
	vwind = vwind_in.sel(season='JJA')
	array = array_in.sel(season='JJA')

	fig = plt.figure()
	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,1])
	lons = uwind.lon
	lats = uwind.lat
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)


	Q = axes[0,1].quiver(xi[::quivkey,::quivkey], yi[::quivkey,::quivkey], uwind[::quivkey,::quivkey], vwind[::quivkey,::quivkey], units='width')

	qk = axes[0,1].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', 
			   labelpos='E', coordinates='figure')


	axes[1,0].set_title('SON')
	uwind = uwind_in.sel(season='SON')
	vwind = vwind_in.sel(season='SON')
	array = array_in.sel(season='SON')

	fig = plt.figure()
	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[1,0])
	lons = uwind.lon
	lats = uwind.lat
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)


	Q = axes[1,0].quiver(xi[::quivkey,::quivkey], yi[::quivkey,::quivkey], uwind[::quivkey,::quivkey], vwind[::quivkey,::quivkey], units='width')

	qk = axes[1,0].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', 
			   labelpos='E', coordinates='figure')

	axes[1,1].set_title('DJF')
	uwind = uwind_in.sel(season='DJF')
	vwind = vwind_in.sel(season='DJF')
	array = array_in.sel(season='DJF')

	fig = plt.figure()
	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[1,1])
	lons = uwind.lon
	lats = uwind.lat
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)


	Q = axes[1,1].quiver(xi[::quivkey,::quivkey], yi[::quivkey,::quivkey], uwind[::quivkey,::quivkey], vwind[::quivkey,::quivkey], units='width')
	qk = axes[1,1].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', 
			   labelpos='E', coordinates='figure')	
	cbar = fig.colorbar(cs,ax=axes)
	cbar.set_label(units)

# not working for some reason.... only saves white space
#	fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/winds_4_seasons_'+str(runmin)+'-'+str(runmax)+'.png', dpi=100)

	plt.show()



def winds_seasons(uwind_in,vwind_in,level,array_in,palette,units,minval,maxval,landmaskxr,outdir,runmin,runmax,units_numerator='m',units_denom='s',quivkey=4,veclen=1.): 

	landlats = np.asarray(landmaskxr.lat)
	landlons = np.asarray(landmaskxr.lon)
	
	fig, axes = plt.subplots(2,2, figsize = (25, 10))

	axes[0,0].set_title('MAM')
	uwind = uwind_in.sel(season='MAM')
	vwind = vwind_in.sel(season='MAM')
	array = array_in.sel(season='MAM')

	landmask = np.asarray(landmaskxr)

#	fig = plt.figure()
	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,0])
	lons = uwind.lon
	lats = uwind.lat
	pres = uwind.pres_lev[level]	
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
	landmask, landlons = addcyclic(landmask, landlons)

	lon, lat = np.meshgrid(lons_shift, lats)
	xi, yi = m(lon, lat)

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)


	Q = axes[0,0].quiver(xi[::quivkey,::quivkey], yi[::quivkey,::quivkey], uwind[level,::quivkey,::quivkey], vwind[level,::quivkey,::quivkey], units='width')
	qk = axes[0,0].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', 
			   labelpos='E', coordinates='figure')


	axes[0,1].set_title('JJA')
	uwind = uwind_in.sel(season='JJA')
	vwind = vwind_in.sel(season='JJA')
	array = array_in.sel(season='JJA')

	fig = plt.figure()
	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,1])
	lons = uwind.lon
	lats = uwind.lat
	pres = uwind.pres_lev[level]	
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)


	Q = axes[0,1].quiver(xi[::quivkey,::quivkey], yi[::quivkey,::quivkey], uwind[level,::quivkey,::quivkey], vwind[level,::quivkey,::quivkey], units='width')

	qk = axes[0,1].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', 
			   labelpos='E', coordinates='figure')


	axes[1,0].set_title('SON')
	uwind = uwind_in.sel(season='SON')
	vwind = vwind_in.sel(season='SON')
	array = array_in.sel(season='SON')

	fig = plt.figure()
	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[1,0])
	lons = uwind.lon
	lats = uwind.lat
	pres = uwind.pres_lev[level]	
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)


	Q = axes[1,0].quiver(xi[::quivkey,::quivkey], yi[::quivkey,::quivkey], uwind[level,::quivkey,::quivkey], vwind[level,::quivkey,::quivkey], units='width')

	qk = axes[1,0].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', 
			   labelpos='E', coordinates='figure')

	axes[1,1].set_title('DJF')
	uwind = uwind_in.sel(season='DJF')
	vwind = vwind_in.sel(season='DJF')
	array = array_in.sel(season='DJF')

	fig = plt.figure()
	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[1,1])
	lons = uwind.lon
	lats = uwind.lat
	pres = uwind.pres_lev[level]	
	uwind, lons_cyclic = addcyclic(uwind, lons)
	vwind, lons_cyclic = addcyclic(vwind, lons)
	uwind = np.asarray(uwind)
	vwind = np.asarray(vwind)
	uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))
	vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
			       cyclic=np.max(lons_cyclic))	

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,1,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	array, lons_cyclic = addcyclic(array, lons)
	array = np.asarray(array)
	array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				     start=False,cyclic=np.max(lons_cyclic))
	array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)

	if palette=='rainnorm':

		if maxval >= minval:
			minval = - maxval
		else: 
			maxval = minval
			minval = - minval

		cs = m.pcolor(xi,yi,array,
			      norm=MidpointNormalize(midpoint=0.),
			      cmap='BrBG',vmin=minval, vmax=maxval)

	elif palette == 'raindefault':
		cs = m.pcolor(xi,yi,array, 
			      cmap=plt.cm.BrBG)

	elif palette=='temp': 
		cs = m.pcolor(xi,yi,array, 
			      norm=MidpointNormalize(midpoint=273.15), 
			      cmap=plt.cm.RdBu_r,vmin = 273.15-(maxval-273.15),
			      vmax=maxval) 

	elif palette=='fromwhite': 
		pal = plt.cm.Blues
		pal.set_under('w',None)
		cs = m.pcolormesh(xi,yi,array,
				  cmap=pal,vmin=0,vmax=maxval)

	else:
		cs = m.pcolor(xi,yi,array)


	Q = axes[1,1].quiver(xi[::quivkey,::quivkey], yi[::quivkey,::quivkey], uwind[level,::quivkey,::quivkey], vwind[level,::quivkey,::quivkey], units='width')
	qk = axes[1,1].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', 
			   labelpos='E', coordinates='figure')	
	cbar = fig.colorbar(cs,ax=axes)
	cbar.set_label(units)

# not working for some reason.... only saves white space
#	plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/winds_4_seasons_'+str(runmin)+'-'+str(runmax)+'.png', bbox_inches='tight', dpi=100)

	plt.show()








def winds_anomaly(uwind,vwind,landmaskxr,landlats,landlons,level=39, minval = -8, maxval = 8):
# colors = u wind minus zonal avg wind
# vectors = absolute wind field 
	landmask = np.asarray(landmaskxr)

	fig = plt.figure(figsize = (25,10))
	m = Basemap(projection='kav7',lon_0=0.,resolution='c')
	lons = uwind.lon
	lats = uwind.lat
	pres = uwind.pres_lev[level]	
	uwind = uwind[level,:,:]
	vwind = vwind[level,:,:]
	uwind,lons_shift = shiftgrid(np.max(lons)-180.,uwind,lons,start=False,
			       cyclic=np.max(lons))
	vwind,lons_shift = shiftgrid(np.max(lons)-180.,vwind,lons,start=False,
			       cyclic=np.max(lons))	
	uwind, lons_cyclic = addcyclic(uwind, lons_shift)
	vwind, lons_cyclic = addcyclic(vwind, lons_shift)

	lon, lat = np.meshgrid(lons_cyclic, lats)
	xi, yi = m(lon, lat)

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	plt.title('Wind at '+str(pres)+' hPa')
	uwind = xr.DataArray(uwind)
	zonavg = uwind.mean(dim='dim_1')
	zonavg = np.expand_dims(zonavg,axis=1)
	zonavg = np.repeat(zonavg,lons_cyclic.shape[0],axis=1)
	array = uwind - zonavg

	# array, array_lons = shiftgrid(np.max(lons)-180.,array,array.lon,
	# 			     start=False,cyclic=np.max(lons))
	# array, lons_cyclic = addcyclic(array, lons)
	# array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

	landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
	landmask, lons_cyclic = addcyclic(landmask, landlons)

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)
	
	cs = m.pcolor(xi,yi,array, 
		      norm=MidpointNormalize(midpoint=0), cmap=plt.cm.RdBu_r, 
		      vmin = minval,
		      vmax = maxval)

	cbar = m.colorbar(cs, location='right', pad="10%")
	cbar.set_label('m/s anomaly')

	Q = plt.quiver(xi[::5,::5], yi[::5,::5], uwind[::5,::5], vwind[::5,::5], units='width')
	qk = plt.quiverkey(Q, 0.9, 0.9, 10, r'$10 \frac{m}{s}$', 
			   labelpos='E', coordinates='figure')

	plt.show()
	return fig

def winds_anomaly_uv_vectors(uwind,vwind,landmaskxr,landlats,landlons,level=39):
# colours = u wind - zonal avg u wind 
# vectors = (u - u_zonavg, v - v_zonavg)

	landmask = np.asarray(landmaskxr)

	fig = plt.figure(figsize = (25,10))
	m = Basemap(projection='kav7',lon_0=0.,resolution='c')
	lons = uwind.lon
	lats = uwind.lat
	pres = uwind.pres_lev[level]	
	uwind = uwind[level,:,:]
	vwind = vwind[level,:,:]
	uwind,lons_shift = shiftgrid(np.max(lons)-180.,uwind,lons,start=False,
			       cyclic=np.max(lons))
	vwind,lons_shift = shiftgrid(np.max(lons)-180.,vwind,lons,start=False,
			       cyclic=np.max(lons))	
	uwind, lons_cyclic = addcyclic(uwind, lons_shift)
	vwind, lons_cyclic = addcyclic(vwind, lons_shift)

	lon, lat = np.meshgrid(lons_cyclic, lats)
	xi, yi = m(lon, lat)

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

	plt.title('Wind at '+str(pres)+' hPa')
	uwind = xr.DataArray(uwind)
	zonavg_u = uwind.mean(dim='dim_1')
	zonavg_u = np.expand_dims(zonavg_u,axis=1)
	zonavg_u = np.repeat(zonavg_u,lons_cyclic.shape[0],axis=1)
	array = uwind - zonavg_u

	vwind = xr.DataArray(vwind)
	zonavg_v = vwind.mean(dim='dim_1')
	zonavg_v = np.expand_dims(zonavg_v,axis=1)
	zonavg_v = np.repeat(zonavg_v,lons_cyclic.shape[0],axis=1)


	landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
	landmask, lons_cyclic = addcyclic(landmask, landlons)

	if np.any(landmask != 0.):
		m.contour(xi,yi,landmask, 1)
# just to get white background
	pal = plt.cm.Blues
	pal.set_under('w',None)
	cs = m.pcolor(xi,yi,array*0., cmap=pal,vmin=0.,vmax=0.)

	Q = plt.quiver(xi[::3,::3], yi[::3,::3], (uwind-zonavg_u)[::3,::3], (vwind-zonavg_v)[::3,::3], units='width')
	qk = plt.quiverkey(Q, 0.9, 0.9, 10, r'$10 \frac{m}{s}$', 
			   labelpos='E', coordinates='figure')

	plt.show()
	return fig








def plot_a_climatology(clim_field,area_array,landmaskxr):
# Plots each month of the climatology as separate curve
# Only works for climatology input field (i.e. 'dim_0' = 12)	
	lats = landmaskxr.lat
	lons = landmaskxr.lon
	time = np.linspace(1,12,12)
	clim_field = xr.DataArray(clim_field, coords = [time, lats, lons], dims = ['time','lat','lon'])
	zonavg = np.empty([12,len(lats)])
	for i in range (0,12):
		zonavg[i,:,] = area_weighted_avg(clim_field[i,:,:],area_array,landmaskxr,option = 'all_sfcs',minlat=-90.,maxlat=90.,axis=1)
	
	colors = ['Yellow','Orange','r','Purple','Blue','Gray',
		  'Brown','Black','Green','Cyan','Teal','Navy']
	
	for i in range (0,12):
		plt.plot(lats,zonavg[i,:],colors[i],label='Month '+str(i+1))
	plt.legend()
	plt.xlabel('Latitude')
	plt.ylabel('Ocean Heat Transport (W/m^2)')
	plt.show()


# def global_areaweighted_avg_2d(array):
	
# 	lats = array.lat
# 	lons = array.lon
# 	latr = np.deg2rad(lats)
# 	# zonalmean = array.mean(dim = 'lon')
# 	# globavg_check = np.average(zonalmean,  weights = (np.cos(latr)))
	
# 	# print(globavg_check)

# 	weights = (np.cos(latr))/(np.sum(np.cos(latr))*len(lons)) # grid cell area/total sfc area --> only works if taking the average over the entire globe, doesn't work to select landmask areas!

# 	weights_2d = np.expand_dims(weights, axis = 1) # for lons
# 	weights_2d = np.repeat(weights_2d, len(lons), axis = 1)
	
# 	weights_2d = xr.DataArray(weights_2d, coords=[lats,lons],dims=['lat','lon'])

# 	# weights_2d.plot()
# 	# plt.show()
# 	# print(np.sum(weights_2d))

# 	array.plot()
# 	plt.show()

	
# 	weighted_array = xr.DataArray((array*weights_2d), coords=[lats,lons],dims=['lat','lon'])
# 	globavg = weighted_array.sum()

# 	# cos_1d = np.cos(latr)
# 	# cos_2d = np.expand_dims(cos_1d, axis = 1) # for lons
# 	# cos_2d = np.repeat(cos_2d, len(lons), axis = 1)

# 	# check2 = np.average(array, cos_2d)

# 	return globavg, weighted_array #, globint 


def area_weighted_avg(array,area_array,landmaskxr,option,minlat=-90.,maxlat=90.,axis=None):

	lats = array.lat
	lons = array.lon
	landmask = np.asarray(landmaskxr)
	array = xr.DataArray(array, coords=[lats,lons], dims = ['lat','lon'])
	area_array = xr.DataArray(area_array, coords=[lats,lons], dims = ['lat','lon'])

	if (minlat!=-90. and maxlat!=90.):
		array = array.sel(lat=slice(minlat,maxlat))
		landmask = np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))
		area_array = area_array.sel(lat=slice(minlat,maxlat))

	if option=='all_sfcs': # meaning both land and ocean
		w_avg = np.average(array, axis=axis, weights=area_array)
#		w_avg_check = (np.sum(array*area_array))/np.sum(area_array) # is the saame
#		aquaplanet_plot(-90.,90.,array,'1','all_sfc',0)

	elif option=='ocean':
		ma = np.ma.array(array, mask=np.isnan(array.where(landmask!=1.)))
		w_avg = np.ma.average(ma,axis=axis,weights=area_array)


		# ma = xr.DataArray(ma, coords=[lats,lons], dims = ['lat','lon'])
		# aquaplanet_plot(-90.,90.,ma,'1','ocean masked array',0)

		# xr.DataArray(mo).plot()
		# plt.show()	
		# plt.close()
	elif option=='land': 
		ma = np.ma.array(array, mask=np.isnan(array.where(landmask==1.)))
		w_avg = np.ma.average(ma,axis=axis,weights=area_array)

		# ma = xr.DataArray(ma, coords=[lats,lons], dims = ['lat','lon'])

		# aquaplanet_plot(-90.,90.,ma,'1','land masked array',0)

		# why can't I plot land mask and then ocean mask or vice  
		# versa if calling ocean option first... ?
		# is it a plotting issue or is the mask then wrong for the 2nd option?
		# xr.DataArray(ml).plot() 
		# plt.show()
		# plt.close()
	return w_avg



def area_weighted_avg_4D(array,area_array,landmaskxr,option,minlat=-90.,maxlat=90):

	num_dims = np.size(np.shape(array))

	lats = array.lat
	lons = array.lon
	dim0_name = array.dims[0]
	dim_0 = array[dim0_name]
	landmask = np.asarray(landmaskxr)
	axis = (num_dims - 2, num_dims - 1)

	if num_dims == 3:
		area_array = xr.DataArray(area_array, coords=[dim_0,lats,lons], dims = [dim0_name,'lat','lon'])

	elif num_dims == 4:
		dim1_name = array.dims[1]
		dim_1 = array[dim1_name]
		area_array = xr.DataArray(area_array, coords=[dim_0,dim_1,lats,lons], dims = [dim0_name,dim1_name,'lat','lon'])

	if (minlat!=-90. and maxlat!=90.):
		array = array.sel(lat=slice(minlat,maxlat))
		landmask = np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))
		area_array = area_array.sel(lat=slice(minlat,maxlat))

	if option=='all_sfcs': # meaning both land and ocean
		w_avg = np.average(array, axis=axis, weights=area_array)
#		w_avg_check = (np.sum(array*area_array))/np.sum(area_array) # is the saame
#		aquaplanet_plot(-90.,90.,array,'1','all_sfc',0)

	elif option=='ocean':
		ma = np.ma.array(array, mask=np.isnan(array.where(landmask!=1.)))
		w_avg = np.ma.average(ma,axis=axis,weights=area_array)


		# ma = xr.DataArray(ma, coords=[lats,lons], dims = ['lat','lon'])
		# aquaplanet_plot(-90.,90.,ma,'1','ocean masked array',0)

		# xr.DataArray(mo).plot()
		# plt.show()	
		# plt.close()
	elif option=='land': 
		ma = np.ma.array(array, mask=np.isnan(array.where(landmask==1.)))
		w_avg = np.ma.average(ma,axis=axis,weights=area_array)

		# ma = xr.DataArray(ma, coords=[lats,lons], dims = ['lat','lon'])

		# aquaplanet_plot(-90.,90.,ma,'1','land masked array',0)

		# why can't I plot land mask and then ocean mask or vice  
		# versa if calling ocean option first... ?
		# is it a plotting issue or is the mask then wrong for the 2nd option?
		# xr.DataArray(ml).plot() 
		# plt.show()
		# plt.close()

	if num_dims == 3:
		return xr.DataArray(w_avg, coords = [dim_0], dims = [dim0_name])
	elif num_dims == 4: 
		return xr.DataArray(w_avg, coords = [dim_0,dim_1], dims = [dim0_name,dim1_name])



def area_integral(array,area_array,landmaskxr,option,minlat=-90.,maxlat=90.,factor=1.):

 	lats = array.lat
	lons = array.lon
	landmask = np.asarray(landmaskxr)
	area_array = xr.DataArray(area_array, coords=[lats,lons], dims = ['lat','lon'])


	array = area_array*array*factor
	

	if (minlat!=-90. and maxlat!=90.):
		array = array.sel(lat=slice(minlat,maxlat))
		landmask = np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))
		area_array = area_array.sel(lat=slice(minlat,maxlat))

	if option=='all_sfcs': # meaning both land and ocean 
		integ = array.sum()
	elif option=='ocean':
		integ = (array.where(landmask!=1.)).sum()
	elif option=='land': 
		integ = (array.where(landmask==1.)).sum()
	return integ






	# latr = np.deg2rad(lats)
	# cos_1d = np.cos(latr)
	# cos_2d = np.expand_dims(cos_1d, axis = 1) # for lons
	# cos_2d = np.repeat(cos_2d, len(lons), axis = 1)	
	# cos_2d = xr.DataArray(cos_2d, coords=[lats,lons], dims = ['lat','lon'])
	
	# landmask = np.asarray(landmaskxr)

	# if (minlat!=-90. and maxlat!=90.):
	# 	array = array.sel(lat=slice(minlat,maxlat))
	# 	landmask = np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))
	# 	cos_2d = cos_2d.sel(lat=slice(minlat,maxlat))

	# if option=='all_sfcs': # meaning both land and ocean NOT entier globe 
	# 	ma =  np.ma.MaskedArray(array, mask=np.isnan(array.where(landmask!=1.)))
	# 	w_avg = np.average(ma, axis=axis,weights=cos_2d)
	# elif option=='ocean':
	# 	ma =  np.ma.MaskedArray(array, mask=np.isnan(array.where(landmask!=1.)))
	# 	w_avg = np.ma.average(ma,axis=axis,weights=cos_2d)
	# elif option=='land': 
	# 	ma = np.ma.MaskedArray(array, mask=np.isnan(array.where(landmask==1.)))
	# 	w_avg = np.ma.average(ma,axis=axis,weights=cos_2d)
	# return w_avg


def mass_streamfunction(testdir,model,runmin,runmax,v='vcomp', a=6376.0e3, g=9.8):
    """Calculate the mass streamfunction for the atmosphere.
    Based on a vertical integral of the meridional wind.
    Ref: Physics of Climate, Peixoto & Oort, 1992.  p158.
    Parameters
    ----------
    data :  xarray.DataSet
        GCM output data
    v : str, optional
        The name of the meridional flow field in `data`.  Default: 'vcomp'
    a : float, optional
        The radius of the planet. Default: Earth 6317km
    g : float, optional
        Surface gravity. Default: Earth 9.8m/s^2
    Returns
    -------
    streamfunction : xarray.DataArray
        The meridional mass streamfunction.

    Author = James Penn
    """

    for i in range(runmin,runmax):
	    if model=='isca':
		    runnr="{0:04}".format(i)
	    elif model=='gfdl':
		    runnr="{0:03}".format(i)
	    data = xr.open_dataset('/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_monthly.nc')
	    vbar = data[v].mean('lon')
	    if i == runmin:
		    msf = np.empty_like(vbar)
		    msf = msf.repeat((runmax-runmin),axis=0) 
	    c = 2*np.pi*a*np.cos(vbar.lat*np.pi/180) / g
# take a diff of half levels, and assign to pfull coordinates
	    dp = xr.DataArray(data.phalf.diff('phalf').values*100, coords=[('pfull', data.pfull)])
	    msf[i-runmin] = (np.cumsum(vbar*dp, axis=vbar.dims.index('pfull')))*c
	    # msf[i-runmin] = (np.cumsum(vbar*dp, axis='pfull')*c)
	    # why cumsum and not # (np.sum(vbar*dp, axis = vbar.dims.index('pfull')))*c
    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))]
    msf=xr.DataArray(msf*1e-10,coords=[time[0],data.pfull,data.lat],dims=['time','pfull','lat'])
    msf_avg=msf.mean(dim='time')
    msf_seasonal_avg=msf.groupby('time.season').mean('time') 
    msf_month_avg=msf.groupby('time.month').mean('time')
    return (msf,msf_avg,msf_seasonal_avg,msf_month_avg)


# In [57]: for i in range(runmin,runmax):
#     ...:         runnr="{0:04}".format(i)
#     ...:         data = xr.open_dataset('/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_monthly.nc')
#     ...:         if i == runmin:
#     ...:             msf = np.empty_like(data[v].mean('lon'))
#     ...:         vbar = data[v].mean('lon')
#     ...:         msf = msf.repeat((runmax-runmin),axis=0)
#     ...:         c = 2*np.pi*a*np.cos(vbar.lat*np.pi/180) / g
#     ...: # take a diff of half levels, and assign to pfull coordinates
#     ...:         dp = xr.DataArray(data.phalf.diff('phalf').values*100, coords=[('pfull', data.pfull)])
#     ...:         print(i-runmin)
#     ...:         msf[i-runmin] = (np.cumsum(vbar*dp, axis='pfull')*c)[0]
#     ...:      

def plot_streamfunction(msf_array,title,units='10^10 kg/s'):

	matplotlib.rcParams['contour.negative_linestyle']= 'dashed'

	lats = msf_array.lat
	pfull = msf_array.pfull
	
	y, p = np.meshgrid(lats, pfull)

	cset1 = plt.contourf(y, p, msf_array, norm=MidpointNormalize(midpoint=0.),
                     cmap='RdBu_r', vmin=-20, vmax=20)

	cont = plt.contour(y,p,msf_array, 20, colors = 'k', linewidth=5)
	plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=14)
	cbar = plt.colorbar(cset1)
	cbar.set_label(units)
	plt.title(title)
	plt.xlabel('Latitude N')
	plt.ylabel('Pressure (hPa)')
	plt.gca().invert_yaxis()
	plt.show()

def plot_streamfunction_seasonal(msf_array,units='10^10 kg/s'):
#colorbar not properly working (limits are weird)
	matplotlib.rcParams['contour.negative_linestyle']= 'dashed'

	lats = msf_array.lat
	pfull = msf_array.pfull
	
	y, p = np.meshgrid(lats, pfull)


	fig, axes = plt.subplots(2,2, sharey = True, figsize = (25, 10))
	fig.gca().invert_yaxis()

	cset = axes[0,0].contourf(y, p, msf_array.sel(season='MAM'), 
				   norm=MidpointNormalize(midpoint=0.),cmap='RdBu_r',vmin = -30.,vmax = 30.)
	cont = axes[0,0].contour(y, p, msf_array.sel(season='MAM'), 20, colors = 'k', linewidth=5)
	axes[0,0].clabel(cont, inline=2, fmt='%1.1f',fontsize=14)
	axes[0,0].set_title('MAM')
	axes[0,0].set_xlabel('Latitude N')
	axes[0,0].set_ylabel('Pressure (hPa)')

	cset = axes[0,1].contourf(y, p, msf_array.sel(season='JJA'), 
				   norm=MidpointNormalize(midpoint=0.),cmap='RdBu_r',vmin = -30.,vmax = 30.)
	cont = axes[0,1].contour(y, p, msf_array.sel(season='JJA'), 20, colors = 'k', linewidth=5)
	axes[0,1].clabel(cont, inline=2, fmt='%1.1f',fontsize=14)
	axes[0,1].set_title('JJA')
	axes[0,1].set_xlabel('Latitude N')
	axes[0,1].set_ylabel('Pressure (hPa)')

	cset3 = axes[1,1].contourf(y, p, msf_array.sel(season='DJF'), 
				   norm=MidpointNormalize(midpoint=0.),cmap='RdBu_r',vmin = -30.,vmax = 30.)
	cont = axes[1,1].contour(y, p, msf_array.sel(season='DJF'), 20, colors = 'k', linewidth=5)
	axes[1,1].clabel(cont, inline=2, fmt='%1.1f',fontsize=14)
	axes[1,1].set_title('DJF')
	axes[1,1].set_xlabel('Latitude N')
	axes[1,1].set_ylabel('Pressure (hPa)')

	cset = axes[1,0].contourf(y, p, msf_array.sel(season='SON'), 
				   norm=MidpointNormalize(midpoint=0.),cmap='RdBu_r',vmin = -30.,vmax = 30.)
	cont = axes[1,0].contour(y, p, msf_array.sel(season='SON'), 20, colors = 'k', linewidth=5)
	axes[1,0].clabel(cont, inline=2, fmt='%1.1f',fontsize=14)
	axes[1,0].set_title('SON')
	axes[1,0].set_xlabel('Latitude N')
	axes[1,0].set_ylabel('Pressure (hPa)')


	cbar = fig.colorbar(cset3,ax=axes) # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
	cbar.set_label(units)

	plt.show()




def rh_P_E_T(outdir,runmin,runmax,rh_avg,precipitation_avg,net_lhe_avg,tsurf_avg,landmask):



	small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
	med = 18 #largefonts 18 # smallfonts 14 # medfonts = 16
	lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

	rh_avg_tropical_land = rh_avg.where(landmask==1.).sel(lat=slice(-30.,30.))
	precip_avg_tropical_land = precipitation_avg.where(landmask==1.).sel(lat=slice(-30.,30.))
	evap_avg_tropical_land = net_lhe_avg.where(landmask==1.).sel(lat=slice(-30.,30.))
	tsurf_avg_tropical_land = tsurf_avg.where(landmask==1.).sel(lat=slice(-30.,30.))

	rh_land_1d = np.asarray(rh_avg_tropical_land).flatten()
	P_land_1d = np.asarray(precip_avg_tropical_land).flatten()
	E_land_1d = np.asarray(evap_avg_tropical_land).flatten()
	tsurf_land_1d = np.asarray(tsurf_avg_tropical_land).flatten()

	rh_avg_tropical_ocean = rh_avg.where(landmask==0.).sel(lat=slice(-30.,30.))
	precip_avg_tropical_ocean = precipitation_avg.where(landmask==0.).sel(lat=slice(-30.,30.))
	evap_avg_tropical_ocean = net_lhe_avg.where(landmask==0.).sel(lat=slice(-30.,30.))
	tsurf_avg_tropical_ocean = tsurf_avg.where(landmask==0.).sel(lat=slice(-30.,30.))

	rh_ocean_1d = np.asarray(rh_avg_tropical_ocean).flatten()
	P_ocean_1d = np.asarray(precip_avg_tropical_ocean).flatten()
	E_ocean_1d = np.asarray(evap_avg_tropical_ocean).flatten()
	tsurf_ocean_1d = np.asarray(tsurf_avg_tropical_ocean).flatten()


	

# 	rh_P_E = np.stack((rh_1d,P_1d,E_1d),axis=1)

	mask = ~np.isnan(rh_land_1d)
	[slope, intercept, r_value, p_value, std_err] = stats.linregress(rh_land_1d[mask],P_land_1d[mask])
	line_P = slope*rh_land_1d + intercept
	[slope, intercept, r_value, p_value, std_err] = stats.linregress(rh_land_1d[mask],E_land_1d[mask])
	line_E = slope*rh_land_1d + intercept

# 	plt.plot(rh_land_1d,P_land_1d,'b.',label = 'P land')
# 	plt.plot(rh_land_1d,E_land_1d,'g.',label = 'E land')
# #	plt.plot(rh_1d, line_P, 'b', label = 'P_regr')
# #	plt.plot(rh_1d, line_E, 'g', label = 'E_regr')

# 	plt.legend()
	# plt.xlabel('RH %')
	# plt.ylabel('P and E (mm/d)')
	# plt.title('P and E versus RH, annual mean (land only)')
	# manager = plt.get_current_fig_manager()
	# manager.window.showMaximized()	
	# plt.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/RH_P_E_land_'+str(runmin)+'-'+str(runmax), bbox_inches='tight', dpi=100)
	# plt.show()

	fig, ax = plt.subplots(1,2,sharey=True,figsize=(25,10))
	ax[0].plot(rh_land_1d,P_land_1d,'b.',label = 'P land (tropics)')
	ax[0].plot(rh_land_1d,E_land_1d,'g.',label = 'E land (tropics)')
	ax[1].plot(rh_ocean_1d,P_ocean_1d,'b.',label = 'P ocean (tropics)')
	ax[1].plot(rh_ocean_1d,E_ocean_1d,'g.',label = 'E ocean (tropics)')
	ax[0].set_xlabel("RH (%)",fontsize = lge)
	ax[1].set_xlabel("RH (%)",fontsize = lge)
	ax[0].tick_params(labelsize = med)
	ax[1].tick_params(labelsize = med)
	ax[0].legend(fontsize = lge)
	ax[1].legend(fontsize = lge)
	ax[0].set_ylabel('P and E (mm/day)',fontsize = lge)

#	fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/RH_P_E_land_oc_'+str(runmin)+'-'+str(runmax)+'_highres.pdf', bbox_inches='tight', dpi=400)
	fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/RH_P_E_land_oc_'+str(runmin)+'-'+str(runmax)+'.png', bbox_inches='tight', dpi=100)

	fig, ax = plt.subplots(1,2, figsize=(25,10))
	ax[0].plot(rh_land_1d,P_land_1d,'b.',label = 'P land (tropics)')
	ax[0].plot(rh_land_1d,E_land_1d,'g.',label = 'E land (tropics)')
	ax[1].plot(rh_land_1d,tsurf_land_1d,'r.',label = 'T land (tropics)')
	ax[0].set_xlabel("RH (%)",fontsize = lge)
	ax[1].set_xlabel("RH (%)",fontsize = lge)
	ax[0].tick_params(labelsize = med)
	ax[1].tick_params(labelsize = med)
	ax[0].legend(fontsize = lge)
	ax[1].legend(fontsize = lge)
	ax[0].set_ylabel('P and E (mm/day)',fontsize = lge)
	ax[1].set_ylabel('tsurf (K)',fontsize = lge)

	fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/RH_P_E_T_land_'+str(runmin)+'-'+str(runmax)+'.png', bbox_inches='tight', dpi=100)

	fig, ax = plt.subplots(1,2, figsize=(25,10))
	ax[0].plot(rh_ocean_1d,P_ocean_1d,'b.',label = 'P ocean (tropics)')
	ax[0].plot(rh_ocean_1d,E_ocean_1d,'g.',label = 'E ocean (tropics)')
	ax[1].plot(rh_ocean_1d,tsurf_ocean_1d,'r.',label = 'T ocean (tropics)')
	ax[0].set_xlabel("RH (%)",fontsize = lge)
	ax[1].set_xlabel("RH (%)",fontsize = lge)
	ax[0].tick_params(labelsize = med)
	ax[1].tick_params(labelsize = med)
	ax[0].legend(fontsize = lge)
	ax[1].legend(fontsize = lge)
	ax[0].set_ylabel('P and E (mm/day)',fontsize = lge)
	ax[1].set_ylabel('tsurf (K)',fontsize = lge)

	fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/RH_P_E_T_ocean_'+str(runmin)+'-'+str(runmax)+'.png', bbox_inches='tight', dpi=100)


	fig, ax = plt.subplots(2,1,sharex = True,figsize = (15,10))
	rh_oc_1d = np.asarray(rh_avg.where(landmask==0.).sel(lat=slice(-30.,30.))).flatten()
	PE_avg = precipitation_avg - net_lhe_avg
	PE_oc_1d = np.asarray(PE_avg.where(landmask==0.).sel(lat=slice(-30.,30.))).flatten()
	PE_land_1d = np.asarray(PE_avg.where(landmask==1.).sel(lat=slice(-30.,30.))).flatten()
	PE_all_1d = np.asarray(PE_avg.sel(lat=slice(-30.,30.))).flatten()
	rh_all_1d = np.asarray(rh_avg.sel(lat=slice(-30.,30.))).flatten()

	ax[0].plot(rh_land_1d,PE_land_1d,'g.', label='P-E land (tropics)')
	ax[1].plot(rh_oc_1d,PE_oc_1d,'b.', label='P-E ocean (tropics)')
#	ax[2].plot(rh_all_1d,PE_all_1d,'k.', label='P-E allsfcs')
#	ax[0].legend()
#	ax[1].legend()
#	ax[2].legend()
#	ax[2].set_xlabel("RH %")
	ax[1].set_xlabel("RH (%)",fontsize = lge)
	ax[0].tick_params(labelsize = med)
	ax[1].tick_params(labelsize = med)
	ax[0].legend(fontsize = lge)
	ax[1].legend(fontsize = lge)
	ax[0].set_ylabel('P-E (mm/day)',fontsize = lge)
	ax[1].set_ylabel('P-E (mm/day)',fontsize = lge)

#	fig.suptitle('P-E versus RH')
#	fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/RH_PE_land_oc_all_'+str(runmin)+'-'+str(runmax)+'_highres.pdf', bbox_inches='tight', dpi=400)
	fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/RH_PE_land_oc_all_'+str(runmin)+'-'+str(runmax)+'.png', bbox_inches='tight', dpi=100)
#	fig.show()

	fig2, ax2 = plt.subplots(3,1,sharex = True,figsize = (25,10))
	P_all_1d = np.asarray(precipitation_avg.sel(lat=slice(-30.,30.))).flatten()
	E_all_1d = np.asarray(net_lhe_avg.sel(lat=slice(-30.,30.))).flatten()

	ax2[0].plot(rh_all_1d,P_all_1d,'b.', label='P tropics')
	ax2[1].plot(rh_all_1d,E_all_1d,'g.', label='E tropics')
	ax2[2].plot(rh_all_1d,PE_all_1d,'k.', label='P-E tropics')
	ax2[0].legend()
	ax2[1].legend()
	ax2[2].legend()
	ax2[2].set_xlabel('RH %')
	ax2[0].set_ylabel('P (mm/d)')
	ax2[1].set_ylabel('E (mm/d)')
	ax2[2].set_ylabel('P-E (mm/d)')
	fig2.suptitle('P, E and P-E versus RH (tropics)')
	fig2.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/RH_P_E_all_'+str(runmin)+'-'+str(runmax)+'.png', bbox_inches='tight', dpi=100)
#	fig2.show()

# NB if I need a (2,2) subplot, need to use ax[0, 1] -- there has to be a space between the comma and the 1 otherwise it doesn't work! i.e. .plot(...) doesn't work and also can't access the subplot that I want! OR fig, axes = plt.subplots(2,2....) and then can call axes[0,0]...

def rh_P_E_change(outdir,runmin,runmax,rh_avg,rh_avg_ctl,precipitation_avg,precipitation_avg_ctl,net_lhe_avg,net_lhe_avg_ctl,tsurf_avg,tsurf_avg_ctl,landmask,sfc='all',minlat=-30.,maxlat=30.):
	

	small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
	med = 18 #largefonts 18 # smallfonts 14 # medfonts = 16
	lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

	if sfc == 'all':
		rh_avg = rh_avg.sel(lat=slice(minlat,maxlat))
		P_avg = precipitation_avg.sel(lat=slice(minlat,maxlat))
		T_avg = tsurf_avg.sel(lat=slice(minlat,maxlat))
		E_avg = net_lhe_avg.sel(lat=slice(minlat,maxlat))

		rh_avg_ctl = rh_avg_ctl.sel(lat=slice(minlat,maxlat))
		P_avg_ctl = precipitation_avg_ctl.sel(lat=slice(minlat,maxlat))
		T_avg_ctl = tsurf_avg_ctl.sel(lat=slice(minlat,maxlat))
		E_avg_ctl = net_lhe_avg_ctl.sel(lat=slice(minlat,maxlat))

		rh_1d = np.asarray(rh_avg).flatten()
		P_1d = np.asarray(P_avg).flatten()
		E_1d = np.asarray(E_avg).flatten()
		T_1d = np.asarray(T_avg).flatten()

		rh_ctl_1d = np.asarray(rh_avg_ctl).flatten()
		P_ctl_1d = np.asarray(P_avg_ctl).flatten()
		E_ctl_1d = np.asarray(E_avg_ctl).flatten()
		T_ctl_1d = np.asarray(T_avg_ctl).flatten()

	elif (sfc == 'land') or (sfc == 'ocean'):
		if sfc == 'land':
			condition = 1.
		elif sfc == 'ocean':
			condition = 0.
		
		rh_avg = rh_avg.where(landmask==condition).sel(lat=slice(minlat,maxlat))
		P_avg = precipitation_avg.where(landmask==condition).sel(lat=slice(minlat,maxlat))
		T_avg = tsurf_avg.where(landmask==condition).sel(lat=slice(minlat,maxlat))
		E_avg = net_lhe_avg.where(landmask==condition).sel(lat=slice(minlat,maxlat))

		rh_avg_ctl = rh_avg_ctl.where(landmask==condition).sel(lat=slice(minlat,maxlat))
		P_avg_ctl = precipitation_avg_ctl.where(landmask==condition).sel(lat=slice(minlat,maxlat))
		T_avg_ctl = tsurf_avg_ctl.where(landmask==condition).sel(lat=slice(minlat,maxlat))
		E_avg_ctl = net_lhe_avg_ctl.where(landmask==condition).sel(lat=slice(minlat,maxlat))

		rh_1d = np.asarray(rh_avg).flatten()
		P_1d = np.asarray(P_avg).flatten()
		E_1d = np.asarray(E_avg).flatten()
		T_1d = np.asarray(T_avg).flatten()

		rh_ctl_1d = np.asarray(rh_avg_ctl).flatten()
		P_ctl_1d = np.asarray(P_avg_ctl).flatten()
		E_ctl_1d = np.asarray(E_avg_ctl).flatten()
		T_ctl_1d = np.asarray(T_avg_ctl).flatten()

	mask = ~np.isnan(rh_ctl_1d)

	fig, ax = plt.subplots(1,2, sharey = False, figsize=(25,10))
	ax[0].plot(P_ctl_1d,(P_1d - P_ctl_1d),'g.', label = sfc)
	ax[0].set_xlabel("P control (mm/d)",fontsize = lge)
	ax[1].set_xlabel("T change (K)",fontsize = lge)
	ax[1].plot((T_1d - T_ctl_1d),(P_1d - P_ctl_1d),'r.', label = sfc)
	ax[0].tick_params(labelsize = med)
	ax[1].tick_params(labelsize = med)
	ax[0].legend(fontsize = lge)
	ax[1].legend(fontsize = lge)
	ax[0].set_ylabel('P change (mm/d)', fontsize = lge)
	ax[1].set_ylabel('P change (mm/d)', fontsize = lge)

	ax[0].set_xlim(precipitation_avg.min(),precipitation_avg.max())
	ax[1].set_xlim((tsurf_avg-tsurf_avg_ctl).min(),(tsurf_avg-tsurf_avg_ctl).max())

	[k,dy,r,p,stderr] = linreg(P_ctl_1d[mask],(P_1d - P_ctl_1d)[mask]) # aa = 8.4, dq = -32
	x1 = np.linspace(np.min((P_ctl_1d)[mask]),np.max((P_ctl_1d)[mask]),500)
	y = k*x1 + dy
	ax[0].plot(x1,y,'g-')
	ax[0].annotate('k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r = '+str("%.2f" % r+', p = '+str(p)), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)

	[k,dy,r,p,stderr] = linreg((T_1d - T_ctl_1d)[mask],(P_1d - P_ctl_1d)[mask]) # aa = 8.4, dq = -32
	x1 = np.linspace(np.min((T_1d - T_ctl_1d)[mask]),np.max((T_1d - T_ctl_1d)[mask]),500)
	y = k*x1 + dy
	ax[1].plot(x1,y,'r-')
	ax[1].annotate('k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r = '+str("%.2f" % r+', p = '+str(p)), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)


	
	fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/Delta_P_vs_Delta_T_and_P_control_'+str(runmin)+'-'+str(runmax)+'_'+sfc+'_between_'+str(minlat)+'N_and_'+str(maxlat)+'N.png', bbox_inches='tight', dpi=100)

	fig2, ax2 = plt.subplots(1,2, sharey = False, figsize=(25,10))

	ax2[0].plot((P_ctl_1d-E_ctl_1d),(P_1d - E_1d - (P_ctl_1d-E_ctl_1d)),'k.', label = sfc)
	ax2[0].set_xlabel("P-E control (mm/d)",fontsize = lge)
	ax2[1].set_xlabel("T change (K)",fontsize = lge)
	ax2[1].plot((T_1d - T_ctl_1d),(P_1d - E_1d - (P_ctl_1d-E_ctl_1d)),'r.', label = sfc)
	ax2[0].tick_params(labelsize = med)
	ax2[1].tick_params(labelsize = med)
	ax2[0].legend(fontsize = lge)
	ax2[1].legend(fontsize = lge)
	ax2[0].set_ylabel('P-E change (mm/d)', fontsize = lge)
	ax2[1].set_ylabel('P-E change (mm/d)', fontsize = lge)

	ax2[0].set_xlim((precipitation_avg - net_lhe_avg).min(),(precipitation_avg - net_lhe_avg).max())
	ax2[1].set_xlim((tsurf_avg-tsurf_avg_ctl).min(),(tsurf_avg-tsurf_avg_ctl).max())

	[k,dy,r,p,stderr] = linreg((P_ctl_1d-E_ctl_1d)[mask],(P_1d - E_1d - (P_ctl_1d-E_ctl_1d))[mask]) # aa = 8.4, dq = -32
	x1 = np.linspace(np.min((P_ctl_1d-E_ctl_1d)[mask]),np.max((P_ctl_1d-E_ctl_1d)[mask]),500)
	y = k*x1 + dy
	ax2[0].plot(x1,y,'k-')
	ax2[0].annotate('k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r = '+str("%.2f" % r+', p = '+str(p)), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)

	[k,dy,r,p,stderr] = linreg((T_1d - T_ctl_1d)[mask],(P_1d - E_1d - (P_ctl_1d-E_ctl_1d))[mask]) # aa = 8.4, dq = -32
	x1 = np.linspace(np.min((T_1d - T_ctl_1d)[mask]),np.max((T_1d - T_ctl_1d)[mask]),500)
	y = k*x1 + dy
	ax2[1].plot(x1,y,'r-')
	ax2[1].annotate('k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r = '+str("%.2f" % r+', p = '+str(p)), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)


	fig2.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/Delta_P-E_vs_Delta_T_and_P-E_control_'+str(runmin)+'-'+str(runmax)+'_'+sfc+'_between_'+str(minlat)+'N_and_'+str(maxlat)+'N.png', bbox_inches='tight', dpi=100)


def land_sea_contrast(array_avg, array_avg_ctl, area_array, landmaskxr, outdir, runmin, runmax, units, plot_title, colorbar, minval = None, maxval = None): 

	print(maxval)

	landmask = np.asarray(landmaskxr)
	zonavg_delta_oc = ((array_avg - array_avg_ctl).where(landmask == 0.)).mean('lon')

	lats = array_avg.lat
	lons = array_avg.lon

	zonavg_delta_oc_2D = np.expand_dims(zonavg_delta_oc, axis = 1)
	zonavg_delta_oc_2D = np.repeat(zonavg_delta_oc_2D, len(lons), axis = 1)
	zonavg_delta_oc_2D = xr.DataArray(zonavg_delta_oc_2D, coords = [lats, lons], dims = ['lat','lon'])

	plot = any_configuration_plot(outdir,runmin,runmax,-90.,90.,(array_avg-array_avg_ctl)/zonavg_delta_oc_2D,area_array,units,plot_title,colorbar,landmaskxr,minval = minval, maxval=maxval)

	plot.show()

def make_var_seasonal(var):

	lats = var.lat
	lons = var.lon
	time = var.time
#	time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))]
	var=xr.DataArray((var.values),coords=[time,lats,lons],dims=['time','lat','lon'])
	var_avg=var.mean(dim='time')
	var_seasonal_avg=var.groupby('time.season').mean('time') 
	var_month_avg=var.groupby('time.month').mean('time')

	return(var,var_avg,var_seasonal_avg,var_month_avg,time)


