from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid, addcyclic
import xarray as xr
import pandas as pd
import os
import matplotlib.colors as colors
import sys
sys.path.insert(0, '/scratch/mp586/PYCODES') # personal module
import stats as st
from scipy import signal


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



def globavg_var_timeseries(testdir,varname,runmin,runmax):
    
    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        runnr="{0:03}".format(i)
        filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
        nc = Dataset(filename,mode='r')
        
        var=nc.variables[varname][:]
        
        if i==runmin:
            timeseries=[var.mean()] # make timeseries be a list, not a float so that I can append later
            print(type(timeseries))
        else:
            timeseries.append(var.mean())

    timeseries=np.asarray(timeseries)
    #print(timeseries)
    months=np.linspace(runmin,(runmax-1),timeseries.size)
    
    plt.plot(months,timeseries)
    plt.title('globavg '+varname)
    plt.xlabel('Month #')
    plt.ylabel('Global average')
    plt.grid(b=True)
    plt.show()    
    return(timeseries)

def tropics_severalvars_timeseries_landonly(testdir,varname1,factor1,varname2,factor2,varname3,factor3,height3,runmin,runmax,squareland):

# for squareland this is naturally only for tropics 
# for continents -- selected latitude slice

    if squareland == 'true': 
	    print('Squareland mode')
	    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
	    landmask=landfile.variables['land_mask'][:]
	    landlats=landfile.variables['lat'][:]
	    landlons=landfile.variables['lon'][:]

	    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
		    runnr="{0:03}".format(i)
		    filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
		    nc = Dataset(filename,mode='r')
	    
		    var1=(xr.DataArray(nc.variables[varname1][:])).where(landmask==1.)*factor1
		    var2=(xr.DataArray(nc.variables[varname2][:])).where(landmask==1.)*factor2
		    if height3 != 0:
			    var3=(xr.DataArray(nc.variables[varname3][:,height3,:,:])).where(landmask==1.)*factor3
		    else:
			    var3=(xr.DataArray(nc.variables[varname3][:])).where(landmask==1.)*factor3


		    if i==runmin:
			    timeseries1=[var1.mean()] # make timeseries be a list, not a float so that I can append later
			    timeseries2=[var2.mean()]
			    timeseries3=[var3.mean()]
		    else:
			    timeseries1.append(var1.mean())
			    timeseries2.append(var2.mean())
			    timeseries3.append(var3.mean())



    else:
	    print('Continental landmask')
	    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land.nc',mode='r')
	    landmask=xr.DataArray(landfile.variables['land_mask'][:])
	    landlats=landfile.variables['lat'][:]
	    landlons=landfile.variables['lon'][:]
	    landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon'])

	    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
		    runnr="{0:03}".format(i)
		    filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
		    nc = Dataset(filename,mode='r')

		    lats = nc.variables['lat'][:]
		    lons = nc.variables['lon'][:]
		    
		    var1=(xr.DataArray(nc.variables[varname1][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor1
		    var2=(xr.DataArray(nc.variables[varname2][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor2
		    if height3 != 0:
			    var3=(xr.DataArray(nc.variables[varname3][0,height3,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3
		    else:
			    var3=(xr.DataArray(nc.variables[varname3][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3

		    var1 = var1.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)
		    print(np.shape(var1))
		    var2 = var2.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)
		    var3 = var3.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))==1.)


		    if i==runmin:
			    timeseries1=[var1.mean()] # make timeseries be a list, not a float so that I can append later
			    timeseries2=[var2.mean()]
			    timeseries3=[var3.mean()]
		    else:
			    timeseries1.append(var1.mean())
			    timeseries2.append(var2.mean())
			    timeseries3.append(var3.mean())

		    print(runnr)



    timeseries1=np.asarray(timeseries1)
    timeseries2=np.asarray(timeseries2)
    timeseries3=np.asarray(timeseries3)

    #print(timeseries)
    months=np.linspace(runmin,(runmax-1),timeseries1.size)
    

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
    colors = ('Green', 'Red', 'Blue')
    field = (timeseries1,timeseries2,timeseries3)
    names = (varname1,varname2,varname3)
    i=0
    for ax, color in zip(axes, colors):
        data = field[i]
	if (i==2 and height3 !=0):
		ax.plot(months,data,color=color,label=names[i]+' at pfull = '+str(int(nc.variables['pfull'][height3]))+' hPa')
	else:	
		ax.plot(months,data,color=color,label=names[i])
        ax.set_ylabel(names[i])
	if (i<=1):
		ax.set_ylim(0,timeseries1.max()+1.)
        ax.tick_params(axis='y',colors=color)
        i+=1
    axes[0].set_xlabel('month #')
    plt.legend()
    plt.title('Tropical Land Only')
    plt.show()

def tropics_severalvars_timeseries_oceanonly(testdir,varname1,factor1,varname2,factor2,varname3,factor3,height3,runmin,runmax,squareland):


    if squareland == 'true': 
	    print('Squareland mode')
	    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
	    landmask=landfile.variables['land_mask'][:]
	    landlats=landfile.variables['lat'][:]
	    landlons=landfile.variables['lon'][:]

	    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
		    runnr="{0:03}".format(i)
		    filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
		    nc = Dataset(filename,mode='r')
	    
		    lats = nc.variables['lat'][:]
		    lons = nc.variables['lon'][:]
		    
		    var1=(xr.DataArray(nc.variables[varname1][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor1
		    var2=(xr.DataArray(nc.variables[varname2][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor2
		    if height3 != 0:
			    var3=(xr.DataArray(nc.variables[varname3][0,height3,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3
		    else:
			    var3=(xr.DataArray(nc.variables[varname3][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3

		    
		    var1 = var1.sel(lat=slice(-30.,30.))
		    var2 = var2.sel(lat=slice(-30.,30.))
		    var3 = var3.sel(lat=slice(-30.,30.))

		    if i==runmin:
			    timeseries1=[var1.mean()] # make timeseries be a list, not a float so that I can append later
			    timeseries2=[var2.mean()]
			    timeseries3=[var3.mean()]
		    else:
			    timeseries1.append(var1.mean())
			    timeseries2.append(var2.mean())
			    timeseries3.append(var3.mean())

    else:
	    print('Continental landmask')
	    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land.nc',mode='r')
	    landmask=xr.DataArray(landfile.variables['land_mask'][:])
	    landlats=landfile.variables['lat'][:]
	    landlons=landfile.variables['lon'][:]
	    landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon'])

	    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
		    runnr="{0:03}".format(i)
		    filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
		    nc = Dataset(filename,mode='r')

		    lats = nc.variables['lat'][:]
		    lons = nc.variables['lon'][:]
		    
		    var1=(xr.DataArray(nc.variables[varname1][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor1
		    var2=(xr.DataArray(nc.variables[varname2][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor2
		    if height3 != 0:
			    var3=(xr.DataArray(nc.variables[varname3][0,height3,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3
		    else:
			    var3=(xr.DataArray(nc.variables[varname3][0,:,:],coords=[lats,lons],dims=['lat','lon']))*factor3

		    var1 = var1.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))!=1.)
		    print(np.shape(var1))
		    var2 = var2.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))!=1.)
		    var3 = var3.sel(lat=slice(-30.,30.)).where(np.asarray(landmaskxr.sel(lat=slice(-30.,30.)))!=1.)


		    if i==runmin:
			    timeseries1=[var1.mean()] # make timeseries be a list, not a float so that I can append later
			    timeseries2=[var2.mean()]
			    timeseries3=[var3.mean()]
		    else:
			    timeseries1.append(var1.mean())
			    timeseries2.append(var2.mean())
			    timeseries3.append(var3.mean())

		    print(runnr)



    timeseries1=np.asarray(timeseries1)
    timeseries2=np.asarray(timeseries2)
    timeseries3=np.asarray(timeseries3)

    #print(timeseries)
    months=np.linspace(runmin,(runmax-1),timeseries1.size)
    

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
    colors = ('Green', 'Red', 'Blue')
    field = (timeseries1,timeseries2,timeseries3)
    names = (varname1,varname2,varname3)
    i=0
    for ax, color in zip(axes, colors):
        data = field[i]
	if (i==2 and height3 !=0):
		ax.plot(months,data,color=color,label=names[i]+' at pfull = '+str(int(nc.variables['pfull'][height3]))+' hPa')
	else:	
		ax.plot(months,data,color=color,label=names[i])
        ax.set_ylabel(names[i])
	if (i<=1):
		ax.set_ylim(0,timeseries1.max()+1.)
        ax.tick_params(axis='y',colors=color)
        i+=1
    axes[0].set_xlabel('month #')
    plt.legend()
    plt.title('Tropical Ocean Only')
    plt.show()

def globavg_var_timeseries_total_and_land(testdir,varname,runmin,runmax,factor,squareland):
    
    if squareland == 'true': 
	    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
	    landmask=landfile.variables['land_mask'][:]
	    landlats=landfile.variables['lat'][:]
	    landlons=landfile.variables['lon'][:]
    else:
	    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land.nc',mode='r')
	    landmask=landfile.variables['land_mask'][:]
	    landlats=landfile.variables['lat'][:]
	    landlons=landfile.variables['lon'][:]

    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        runnr="{0:03}".format(i)
        filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
        nc = Dataset(filename,mode='r')
        
        var=nc.variables[varname][:]*factor
	var_land = (xr.DataArray(var)).where(landmask==1.)
	var_ocean = (xr.DataArray(var)).where(landmask!=1.)
        if i==runmin:
            timeseries=[var.mean()] # make timeseries be a list, not a float so that I can append later
	    timeseries_land=[var_land.mean()]
	    timeseries_ocean=[var_ocean.mean()]
          
        else:
            timeseries.append(var.mean())
	    timeseries_land.append(var_land.mean())
	    timeseries_ocean.append(var_ocean.mean())


    timeseries=np.asarray(timeseries)
    timeseries_land=np.asarray(timeseries_land)
    timeseries_ocean=np.asarray(timeseries_ocean)
    months=np.linspace(runmin,(runmax-1),timeseries.size)
    plt.plot(months,timeseries,'r',label='total')
    plt.plot(months,timeseries_land,'g',label='land only')
    plt.plot(months,timeseries_ocean,'b',label='ocean only')
    plt.title('globavg '+varname+' total and land/ocean only')
    plt.xlabel('Month #')
    plt.ylabel('Global average')
    plt.legend()
    plt.grid(b=True)
    plt.show()    
    return(timeseries)


def seasonal_surface_variable(testdir,runmin,runmax,varname,units): # only works for surface variables (dims time, lat, lon) at the moment
    
    plt.close()
    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        runnr="{0:03}".format(i)
        filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
        nc = Dataset(filename,mode='r')
              
        if i==runmin:
            var=xr.DataArray(nc.variables[varname][:]) # only monthly avg for month i
        else:
            var_i=xr.DataArray(nc.variables[varname][:])
            var=xr.concat([var,var_i],'dim_0')
    
    
    lons= nc.variables['lon'][:]
    lats= nc.variables['lat'][:]
    
    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))]
    var=xr.DataArray(var.values,coords=[time[0],lats,lons],dims=['time','lat','lon'])
    var_avg=var.mean(dim='time')
    var_seasonal_avg=var.groupby('time.season').mean('time') 
    var_month_avg=var.groupby('time.month').mean('time')
# for several dimension mean: XY.mean(dim=('lat','lon'))

#     plot=tropics_plot(lats,lons,var_avg,units,'Average '+varname)
#     print(type(plot))
    JJA='JJA'
    DJF='DJF'
    MAM='MAM'
    SON='SON'
#     #does not work if I write .....sel(season='JJA') because the resulting array has the dimensions (1,nrlats,nrlons), but if I select indirectly using the above definitions, it works!
#     plot=tropics_plot(lats,lons,var_seasonal_avg.sel(season=JJA),units,'JJA '+varname)
   

# # this is not producing one panel plot but opens each separately
#     fig=plt.figure()
#     fig.add_subplot(2,2,1)
#     tropics_plot(lats,lons,var_seasonal_avg.sel(season=DJF),units,'DJF '+varname)
#     fig.add_subplot(2,2,2)
#     tropics_plot(lats,lons,var_seasonal_avg.sel(season=MAM),units,'MAM '+varname)
#     fig.add_subplot(2,2,3)
#    # tropics_plot(lats,lons,var_seasonal_avg.sel(season=JJA),units,'JJA '+varname)
#     fig.add_subplot(2,2,4)
#     tropics_plot(lats,lons,var_seasonal_avg.sel(season=SON),units,'SON '+varname)
# # Fine-tune figure; make subplots farther from each other.


#    zonavg_plot(var_seasonal_avg.sel(season=JJA),'zonal mean JJA '+varname,varname,units)

    zonavg_plot(var_avg,'zonal mean avg '+varname,varname,units)   
    
    return(var,var_avg,var_seasonal_avg,var_month_avg,time)


def zonavg_plot(field,title,varname,units): #field needs to have dimensions lat|lon

    field=xr.DataArray(field)
    field_zonavg=field.mean(dim='lon',keep_attrs=True)
    
    plt.close()
    plt.plot(field['lat'],field_zonavg)
    plt.title(title)
    plt.ylabel(varname+' ('+units+')')
    plt.xlabel('Latitude')
    plt.show()

def several_vars_zonalavg(field1,varname1,field2,varname2,field3,varname3): #fields needs to have dimensions lat|lon
    
    field1=xr.DataArray(field1)
    field1_zonavg=field1.mean(dim='lon',keep_attrs=True)
    field2=xr.DataArray(field2)
    field2_zonavg=field2.mean(dim='lon',keep_attrs=True)   
    field3=xr.DataArray(field3)
    field3_zonavg=field3.mean(dim='lon',keep_attrs=True)
    
    plt.close()
    plt.plot(field1['lat'],field1_zonavg,label=varname1)
    plt.plot(field1['lat'],field2_zonavg,label=varname2)
    plt.plot(field1['lat'],field3_zonavg,label=varname3)
    plt.legend()
    plt.xlabel('latitude')
    plt.show()

# adapted from http://stackoverflow.com/questions/7733693/matplotlib-overlay-plots-with-different-scales
# with 3 different y axes
def several_vars_zonalavg2(field1,varname1,field2,varname2,field3,varname3,title):

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
    colors = ('Green', 'Red', 'Blue')
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
    trop_maxreverseindex=np.asarray(np.where(lats[::-1]<=30.))[0,0] # there might be a better way of doing this!
    tropical_lats=lats[trop_minindex:(lats.size-trop_maxreverseindex)]
    lon_0 = lons.mean() # trick to find central lon lat (lon from 0 to 360)
    lat_0 = lats.mean() # ~= 0 because lat from -90 to 90



    m = Basemap(projection='cyl',llcrnrlat=-30,urcrnrlat=30,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')
    m.drawcoastlines()
#m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
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

def squareland_plot(minlat,maxlat,array,units,title,palette):

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
    ax1 = plt.subplot2grid((3,7), (0,0), colspan = 5, rowspan = 3)


    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# draw parallels and meridians.

    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    array, lons_cyclic = addcyclic(array, lons)
    array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
    zonavg_thin = array.mean(dim='lon')
# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.


    lon, lat = np.meshgrid(lons_cyclic, selected_lats)
    xi, yi = m(lon, lat)

    dlons = lons[100] - lons[99]
    dlats = lats[60] - lats[59]

#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG')
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.seismic)
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


# Add Colorbar
    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)
    # sns.palplot(sns.color_palette("BrBG", 7))

# Show landmask
    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    
#    m.contour(xi, yi, landmask, colors='k')

    square_lons,square_lats=(xr.DataArray(np.meshgrid(lons, lats))).where(landmask==1.)
    square_lons = square_lons + 180. + dlons

# adjust square so that it optically fits eg evaporation    
#    square_lons = square_lons + dlons
#    square_lats = square_lats + dlats



# Add rectangle
    # xs = [sq_minlon,sq_maxlon,sq_maxlon,sq_minlon,sq_minlon]
    # ys = [sq_minlat,sq_minlat,sq_maxlat,sq_maxlat,sq_minlat]
    #xs=[square_lons.min(),square_lons.max(),square_lons.max(),square_lons.min(),square_lons.min()]
    #ys=[square_lats.min(),square_lats.min(),square_lats.max()+dlats,square_lats.max()+dlats,square_lats.min()]
    #m.plot(xs, ys)

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
   
    ax2 = plt.subplot2grid((3,7), (0,6), rowspan = 3)
    
    plt.plot(zonavg_thin,array['lat'])
    plt.ylabel('Latitude')
    plt.xlabel(title+' ('+units+')')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.invert_xaxis()
    plt.show()

    plt.show()
    print(square_lons.min(),square_lons.max(),square_lats.min(),square_lats.max()+dlats)



    return(square_lons,square_lats)


def squareland_plot_minuszonavg(minlat,maxlat,array,units,title,palette,zonavgtitle):

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
    ax1 = plt.subplot2grid((3,7), (0,0), colspan = 5, rowspan = 3)

    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# draw parallels and meridians.
    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    
    array, lons_cyclic = addcyclic(array, lons)
    array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.


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
    
# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.


    dlons = lons[100] - lons[99]
    dlats = lats[60] - lats[59]

#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',vmin=minval, vmax=maxval)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG,vmin=minval, vmax=maxval)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.seismic,vmin=minval, vmax=maxval)
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


# Add Colorbar
    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)
    # sns.palplot(sns.color_palette("BrBG", 7))

# Show landmask
    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    
#    m.contour(xi, yi, landmask, colors='k')

    square_lons,square_lats=(xr.DataArray(np.meshgrid(lons, lats))).where(landmask==1.)
    square_lons = square_lons + 180. + dlons

# optically adjust square posiiton os that it fits shape in eg evaporation plot
#    square_lons = square_lons + dlons
#    square_lats = square_lats + dlats


# Add rectangle
    # xs = [sq_minlon,sq_maxlon,sq_maxlon,sq_minlon,sq_minlon]
    # ys = [sq_minlat,sq_minlat,sq_maxlat,sq_maxlat,sq_minlat]
    #xs=[square_lons.min(),square_lons.max(),square_lons.max(),square_lons.min(),square_lons.min()]
    #ys=[square_lats.min(),square_lats.min(),square_lats.max()+dlats,square_lats.max()+dlats,square_lats.min()]
    #m.plot(xs, ys)

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

    ax2 = plt.subplot2grid((3,7), (0,6), rowspan = 3)
    
    plt.plot(zonavg_thin,array['lat'])
    plt.ylabel('Latitude')
    plt.xlabel(zonavgtitle+' ('+units+')')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.invert_xaxis()
    plt.show()
    print(square_lons.min(),square_lons.max(),square_lats.min(),square_lats.max()+dlats)



    return(square_lons,square_lats)

def aquaplanet_plot(minlat,maxlat,array,units,title,palette):

    plt.close()

    lats=array.lat
    lons=array.lon
    print(np.shape(lons))
    
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
    ax1 = plt.subplot2grid((3,7), (0,0), colspan = 5, rowspan = 3)


    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# draw parallels and meridians.
    array, lons_cyclic = addcyclic(array, lons)
    array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.


    lon, lat = np.meshgrid(lons_cyclic, selected_lats)
    xi, yi = m(lon, lat)

    zonavg_thin = array.mean(dim='lon')
    
    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])

#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG')
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.seismic)
    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))


# Add Colorbar
    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)
    # sns.palplot(sns.color_palette("BrBG", 7))



# Add Title
    plt.title(title)


    ax2 = plt.subplot2grid((3,7), (0,6), rowspan = 3)
    
    plt.plot(zonavg_thin,array['lat'])
    plt.ylabel('Latitude')
    plt.xlabel(title+' ('+units+')')
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.invert_xaxis()
    plt.show()


    plt.show()


def aquaplanet_plot_minuszonavg(minlat,maxlat,array,units,title,palette):

    plt.close()

    lats=array.lat
    lons=array.lon
    print(np.shape(lons))
    
    minlatindex=np.asarray(np.where(lats>=minlat))[0,0]
    maxlatreverseindex=np.asarray(np.where(lats[::-1]<=maxlat))[0,0] # there might be a better way of doing this!
    selected_lats=lats[minlatindex:(lats.size-maxlatreverseindex)+1]

    # sqlatindex=np.asarray(np.where(lats>sq_minlat))[0,0]
    # sqlatreverseindex=np.asarray(np.where(lats[::-1]<sq_maxlat))[0,0] # there might be a better way of doing this!
    # square_lats=lats[sqlatindex:(lats.size-sqlatreverseindex)]

    # sqlonindex=np.asarray(np.where(lons>sq_minlon))[0,0]
    # sqlonreverseindex=np.asarray(np.where(lons[::-1]<sq_maxlon))[0,0] # there might be a better way of doing this!
    # square_lons=lons[sqlonindex:(lons.size-sqlonreverseindex)+1]





    m = Basemap(projection='kav7',lon_0=0.,resolution='c')

    array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# draw parallels and meridians.
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


#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',vmin=minval, vmax=maxval)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG,vmin=minval, vmax=maxval)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.seismic,vmin=minval, vmax=maxval)
    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats))


# Add Colorbar
    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)
    # sns.palplot(sns.color_palette("BrBG", 7))



# Add Title
    plt.title(title)


    plt.show()





def squareland_plot_several(minlat,maxlat,array1,title1,array2,title2,array3,title3,units,title,palette):

    plt.close()

    lats=array1.lat
    lons=array1.lon


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
    ax1 = plt.subplot2grid((3,7), (0,0), colspan = 5, rowspan = 3)
    
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

# Show landmask
    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
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


#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array1.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',latlon=True)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array1.sel(lat=selected_lats), cmap=plt.cm.BrBG,latlon=True)
    elif palette=='temp':
        cs = m.pcolor(xi,yi,array1.sel(lat=selected_lats), cmap=plt.cm.seismic,latlon=True)
    else:
        cs = m.pcolor(xi,yi,array1.sel(lat=selected_lats),latlon=True)

    m.plot(xs, ys, latlon = True)

    axes[1].set_title(title2)

    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array2.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',latlon=True)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array2.sel(lat=selected_lats), cmap=plt.cm.BrBG,latlon=True)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array2.sel(lat=selected_lats), cmap=plt.cm.seismic,latlon=True)
    else:
        cs = m.pcolor(xi,yi,array2.sel(lat=selected_lats),latlon=True)

    m.plot(xs, ys, latlon = True)

    axes[2].set_title(title3)

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

# Show landmask
    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
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


    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array3.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',latlon=True)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array3.sel(lat=selected_lats), cmap=plt.cm.BrBG,latlon=True)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array3.sel(lat=selected_lats), cmap=plt.cm.seismic,latlon=True)
    else:
        cs = m.pcolor(xi,yi,array3.sel(lat=selected_lats),latlon=True)

    m.plot(xs, ys, latlon = True)


# Add Title
    plt.title(title)


# Add Colorbar
    cbar = plt.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)

    plt.show()

    return(square_lons,square_lats)




def worldmap_inputfile(filename,varname): # assuming that the input file has the dimenstions time, lat, lon

    nc=Dataset(filename,mode='r')

    array=xr.DataArray(nc.variables[varname][:])

    lats=xr.DataArray(nc.variables['lat'][:])
    lons=xr.DataArray(nc.variables['lon'][:])

    lon_0 = lons.mean() # trick to find central lon lat (lon from 0 to 360)
    lat_0 = lats.mean() # ~= 0 because lat from -90 to 90


    fig = plt.figure()
    ax1 = plt.subplot2grid((3,7), (0,0), colspan = 5, rowspan = 3)


    m = Basemap(projection='kav7',lon_0=180.,resolution='c')

# draw parallels and meridians.
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
    array = array.mean('dim_0')
    array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

    zonavg_thin = array.mean(dim='lon')
# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.

    
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    cs = m.pcolor(xi,yi,array)
    cbar = m.colorbar(cs, location='right', pad="10%")
    

    ax2 = plt.subplot2grid((3,7), (0,6), rowspan = 3)
    
    plt.plot(zonavg_thin,array['lat'])
    plt.ylabel('Latitude')
    plt.xlabel(varname)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position('right')
    ax2.invert_xaxis()
    plt.show()






def worldmap_variable(field,units,title,palette,minval,maxval):
    plt.close()
    
    lats=field.lat
    lons=field.lon
    lon_0 = lons.mean() # trick to find central lon lat (lon from 0 to 360)
    lat_0 = lats.mean() # ~= 0 because lat from -90 to 90



    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
            llcrnrlon=0,urcrnrlon=360,resolution='c')
    m.drawcoastlines()
#m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(0.,361.,60.))


# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)




    if (minval==0) & (maxval==0):

	    if palette=='rainnorm':
		    cs = m.pcolor(xi,yi,field,norm=MidpointNormalize(midpoint=0.),cmap='BrBG', latlon=True)
	    elif palette == 'raindefault':
		    cs = m.pcolor(xi,yi,field, cmap=plt.cm.BrBG,latlon=True)
	    elif palette=='temp': 
		    cs = m.pcolor(xi,yi,field, cmap=plt.cm.seismic,latlon=True)
	    else:
		    cs = m.pcolor(xi,yi,field,latlon=True)
    else:

	    if palette=='rainnorm':
		    cs = m.pcolor(xi,yi,field,norm=MidpointNormalize(midpoint=0.),cmap='BrBG', vmin=minval, vmax=maxval, latlon=True)
	    elif palette == 'raindefault':
		    cs = m.pcolor(xi,yi,field, cmap=plt.cm.BrBG, vmin=minval, vmax=maxval, latlon=True)
	    elif palette=='temp': 
		    cs = m.pcolor(xi,yi,field, cmap=plt.cm.seismic, vmin=minval, vmax=maxval, latlon=True)
	    else:
		    cs = m.pcolor(xi,yi,field, vmin=minval, vmax=maxval, latlon=True)

    


# Add Colorbar
    cbar = m.colorbar(cs, location='right', pad="10%")
    cbar.set_label(units)
    # sns.palplot(sns.color_palette("BrBG", 7))

# Add Title
    plt.title(title)


    plt.show()
    return(cs)


def squareland_inputfile(filename,varname):
	
	nc=Dataset(filename,mode='r')

	array=xr.DataArray(nc.variables[varname][:])
	array = array.mean(dim='dim_0')
	lats=xr.DataArray(nc.variables['lat'][:])
	lons=xr.DataArray(nc.variables['lon'][:])

	fig = plt.figure()
	ax1 = plt.subplot2grid((3,7), (0,0), colspan = 5, rowspan = 3)


	m = Basemap(projection='kav7',lon_0=0.,resolution='c')

	array,lons = shiftgrid(np.max(lons)-180.,array,lons,start=False,cyclic=np.max(lons))
# draw parallels and meridians.

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0])
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1])
	array, lons_cyclic = addcyclic(array, lons)
	array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
	zonavg_thin = array.mean(dim='lon')
# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.


	lon, lat = np.meshgrid(lons_cyclic, lats)
	xi, yi = m(lon, lat)

	dlons = lons[100] - lons[99]
	dlats = lats[60] - lats[59]

      	cs = m.pcolor(xi,yi,array)
# Add Colorbar
	cbar = m.colorbar(cs, location='right', pad="10%")
# sns.palplot(sns.color_palette("BrBG", 7))

# Show landmask
	landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
	landmask=landfile.variables['land_mask'][:]
	landlats=landfile.variables['lat'][:]
	landlons=landfile.variables['lon'][:]
	
#    m.contour(xi, yi, landmask, colors='k')

	square_lons,square_lats=(xr.DataArray(np.meshgrid(landlons, landlats))).where(landmask==1.)
	square_lons = square_lons + dlons

# adjust square so that it optically fits eg evaporation    
#    square_lons = square_lons + dlons
#    square_lats = square_lats + dlats



# Add rectangle
# xs = [sq_minlon,sq_maxlon,sq_maxlon,sq_minlon,sq_minlon]
# ys = [sq_minlat,sq_minlat,sq_maxlat,sq_maxlat,sq_minlat]
    #xs=[square_lons.min(),square_lons.max(),square_lons.max(),square_lons.min(),square_lons.min()]
    #ys=[square_lats.min(),square_lats.min(),square_lats.max()+dlats,square_lats.max()+dlats,square_lats.min()]
    #m.plot(xs, ys)

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

	ax2 = plt.subplot2grid((3,7), (0,6), rowspan = 3)
	
	plt.plot(zonavg_thin,array['lat'])
	plt.ylabel('Latitude')
	plt.xlabel(varname)
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position('right')
	ax2.invert_xaxis()
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
    landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land_square.nc',mode='r')
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

