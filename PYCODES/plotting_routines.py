from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
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




def globavg_tsurf_timeseries(testdir,runmin,runmax):
    
    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        runnr="{0:03}".format(i)
        filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
        nc = Dataset(filename,mode='r')
        
        tsurf=nc.variables['t_surf'][:]
        
        if i==runmin:
            timeseries=[tsurf.mean()] # make timeseries be a list, not a float so that I can append later
            print(type(timeseries))
        else:
            timeseries.append(tsurf.mean())

    timeseries=np.asarray(timeseries)
    #print(timeseries)
    months=np.linspace(runmin,(runmax-1),timeseries.size)
    plot1=plt.plot(months,timeseries)
    plt.title('Tsurf (globavg)')
    plt.xlabel('Month #')
    plt.ylabel('Global average tsurf (K)')
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
    cbar = m.colorbar(cs, location='bottom', pad="10%")
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





    m = Basemap(projection='cyl',llcrnrlat=minlat,urcrnrlat=maxlat,\
            llcrnrlon=-181.,urcrnrlon=180.,resolution='c')
# draw parallels and meridians.
    m.drawparallels(np.arange(minlat,maxlat+1.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180.,181.,60.),labels=[0,0,0,1])


# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
    lon, lat = np.meshgrid(lons, selected_lats)
    xi, yi = m(lon, lat)

#latlon=True --> wrap lons so that map can go from -180 to 180, otherwise there is nothing from -180 to 0

    if palette=='rainnorm':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),norm=MidpointNormalize(midpoint=0.),cmap='BrBG',latlon=True)
    elif palette == 'raindefault':
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.BrBG,latlon=True)
    elif palette=='temp': 
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats), cmap=plt.cm.seismic,latlon=True)
    else:
        cs = m.pcolor(xi,yi,array.sel(lat=selected_lats),latlon=True)




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

# Add rectangle
    # xs = [sq_minlon,sq_maxlon,sq_maxlon,sq_minlon,sq_minlon]
    # ys = [sq_minlat,sq_minlat,sq_maxlat,sq_maxlat,sq_minlat]
    xs=[square_lons.min(),square_lons.max(),square_lons.max(),square_lons.min(),square_lons.min()]
    ys=[square_lats.min(),square_lats.min(),square_lats.max(),square_lats.max(),square_lats.min()]
    m.plot(xs, ys, latlon = True)


# Add Title
    plt.title(title)


    plt.show()

    return(square_lons,square_lats)


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



    fig, axes = plt.subplots(3,1)
    axes[0].set_title(title1)
    
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
    ys=[square_lats.min(),square_lats.min(),square_lats.max(),square_lats.max(),square_lats.min()]


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
    ys=[square_lats.min(),square_lats.min(),square_lats.max(),square_lats.max(),square_lats.min()]


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
       
    lats=np.linspace(-90.,90.,array.shape[1]+1)
    lons=np.linspace(0.,360.,array.shape[2]+1)


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

   
    cs = m.pcolor(xi,yi,array.mean(dim='dim_0'))
 
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
    ys=[square_lats.min(),square_lats.min(),square_lats.max(),square_lats.max(),square_lats.min()]
    m.plot(xs, ys, latlon = True)


# Add Title
    plt.title(title)


    plt.show()

