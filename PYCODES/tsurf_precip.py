from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
import stats as st # personal module 

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
    plt.show()    
    return(timeseries)


def timeavg_precip(testdir,runmin,runmax):
    
    plt.close()
    for i in range (runmin,runmax): # excludes last one! i.e. not from 1 - 12 but from 1 - 11!
        runnr="{0:03}".format(i)
        filename = '/scratch/mp586/GFDL_DATA/'+testdir+'/run'+runnr+'/atmos_monthly.nc'
        nc = Dataset(filename,mode='r')
              
        if i==runmin:
            ccrain=xr.DataArray(nc.variables['convection_rain'][:])+xr.DataArray(nc.variables['condensation_rain'][:]) # only monthly avg for month i
        else:
            ccrain_i=xr.DataArray(nc.variables['convection_rain'][:])+xr.DataArray(nc.variables['condensation_rain'][:])
            ccrain=xr.concat([ccrain,ccrain_i],'dim_0')
    
    
    lons= nc.variables['lon'][:]
    lats= nc.variables['lat'][:]
    
    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))]
    ccrain=xr.DataArray(ccrain.values,coords=[time[0],lats,lons],dims=['time','lat','lon'])
    ccrain_avg=ccrain.mean(dim='time')
    ccrain_seasonal_avg=ccrain.groupby('time.season').mean('time') 

    plot=tropics_plot(lats,lons,ccrain_avg*86400,'mm/day','Average Rainfall')
    print(type(plot))
    JJA='JJA'
    DJF='DJF'
    MAM='MAM'
    SON='SON'
    #does not work if I write .....sel(season='JJA') because the resulting array has the dimensions (1,nrlats,nrlons), but if I select indirectly using the above definitions, it works!
    plot=tropics_plot(lats,lons,ccrain_seasonal_avg.sel(season=JJA)*86400,'mm/day','JJA Rainfall')
   

# this is not producing one panel plot but opens each separately
    fig=plt.figure()
    fig.add_subplot(2,2,1)
    tropics_plot(lats,lons,ccrain_seasonal_avg.sel(season=DJF)*86400,'mm/day','DJF Rainfall')
    fig.add_subplot(2,2,2)
    tropics_plot(lats,lons,ccrain_seasonal_avg.sel(season=MAM)*86400,'mm/day','MAM Rainfall')
    fig.add_subplot(2,2,3)
  #  tropics_plot(lats,lons,ccrain_seasonal_avg.sel(season=JJA)*86400,'mm/day','JJA Rainfall')
    fig.add_subplot(2,2,4)
    tropics_plot(lats,lons,ccrain_seasonal_avg.sel(season=SON)*86400,'mm/day','SON Rainfall')
# Fine-tune figure; make subplots farther from each other.
   

    zonavg_plot(ccrain_seasonal_avg.sel(season=JJA)*86400,'zonal mean JJA Rainfall','Rainfall','mm/day')



    return(ccrain,ccrain_avg,ccrain_seasonal_avg,time)


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
    var_seasonal_avg=var.groupby('time.season').mean('time') # for several dimension mean: XY.mean(dim=('lat','lon'))

    plot=tropics_plot(lats,lons,var_avg,units,'Average '+varname)
    print(type(plot))
    JJA='JJA'
    DJF='DJF'
    MAM='MAM'
    SON='SON'
   

# this is not producing one panel plot but opens each separately
    fig=plt.figure()
    fig.add_subplot(2,2,1)
    tropics_plot(lats,lons,var_seasonal_avg.sel(season=DJF),units,'DJF '+varname)
    fig.add_subplot(2,2,2)
    tropics_plot(lats,lons,var_seasonal_avg.sel(season=MAM),units,'MAM '+varname)
    fig.add_subplot(2,2,3)
    tropics_plot(lats,lons,var_seasonal_avg.sel(season=JJA),units,'JJA '+varname)
    fig.add_subplot(2,2,4)
    tropics_plot(lats,lons,var_seasonal_avg.sel(season=SON),units,'SON '+varname)
# Fine-tune figure; make subplots farther from each other.


    zonavg_plot(var_seasonal_avg.sel(season=JJA),'zonal mean JJA '+varname,varname,units)
    
    
    return(var,var_avg,var_seasonal_avg,time)


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


  #  plt.show()
    return(cs)


# runmin=1
# runmax=360
# testdir='co2_test_86'
# temp=globavg_tsurf_timeseries(testdir,runmin,runmax)



# globavg t after 360 months of spinup : 270K --> yes, depends on model setup, this one is a bit cold.



runmin=85 # careful, first month needs to be january!
runmax=240
testdir='full_qflux'

k=globavg_tsurf_timeseries(testdir,runmin,runmax)

[ccrain,ccrain_avg,ccrain_seasonal_avg,time]=timeavg_precip(testdir,runmin,runmax) #do this one extra because ccrain is not a field in the files, but is the sum of two fields.
[tsurf,tsurf_avg,tsurf_seasonal_avg,time]=seasonal_surface_variable(testdir,runmin,runmax,'t_surf','K')
[net_sw,net_sw_avg,net_sw_seasonal_avg,time]=seasonal_surface_variable(testdir,runmin,runmax,'flux_sw','W/m^2') # shortwave net flux at surface (UP)
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,time]=seasonal_surface_variable(testdir,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[convection_rain,convection_rain_avg,convection_rain_seasonal_avg,time]=seasonal_surface_variable(testdir,runmin,runmax,'convection_rain','kg/m^2s') # latent heat flux at surface (UP)
[condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,time]=seasonal_surface_variable(testdir,runmin,runmax,'condensation_rain','kg/m^2s') # latent heat flux at surface (UP)


JJA='JJA'
DJF='DJF'
MAM='MAM'
SON='SON'

#this is not the normalized cross correlation (values not between -1 and 1)
#from scipy import signal 

[r1,pval1]=st.pattern_corr_2d(ccrain_seasonal_avg.sel(season=DJF),tsurf_seasonal_avg.sel(season=DJF))
print('Spatial correlation of ccrain and tsurf for DJF ='+str(r1))

[r2,pval2]=st.pattern_corr_2d(ccrain_seasonal_avg.sel(season=JJA),tsurf_seasonal_avg.sel(season=JJA))
print('Spatial correlation of ccrain and tsurf for JJA ='+str(r2))

[r3,pval3]=st.pattern_corr_2d(ccrain_seasonal_avg.sel(season=JJA),net_sw_seasonal_avg.sel(season=JJA)) 
print('Spatial correlation of ccrain and net sw for JJA ='+str(r3))

[r4,pval4]=st.pattern_corr_2d(tsurf_seasonal_avg.sel(season=JJA),net_sw_seasonal_avg.sel(season=JJA))
print('Spatial correlation of tsurf and net sw for JJA ='+str(r4))

several_vars_zonalavg2(ccrain_seasonal_avg.sel(season=JJA)*86400.,'ccrain (mm/day)',net_sw_seasonal_avg.sel(season=JJA),'net_sw (W/m^2)',net_lhe_seasonal_avg.sel(season=JJA),'net_lhe (W/m^2)','JJA')

several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=JJA),'tsurf (K)',net_sw_seasonal_avg.sel(season=JJA),'net_sw (W/m^2)',net_lhe_seasonal_avg.sel(season=JJA),'net_lhe (W/m^2)','JJA')

several_vars_zonalavg2(ccrain_seasonal_avg.sel(season=DJF)*86400.,'ccrain (mm/day)',net_sw_seasonal_avg.sel(season=DJF),'net_sw (W/m^2)',net_lhe_seasonal_avg.sel(season=DJF),'net_lhe (W/m^2)','DJF')

several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=DJF),'tsurf (K)',net_sw_seasonal_avg.sel(season=DJF),'net_sw (W/m^2)',net_lhe_seasonal_avg.sel(season=DJF),'net_lhe (W/m^2)','DJF')

several_vars_zonalavg(ccrain_seasonal_avg.sel(season=DJF)*86400.,'ccrain (mm/day)',condensation_rain_seasonal_avg.sel(season=DJF)*86400.,'cond_rain (mm/day)',convection_rain_seasonal_avg.sel(season=DJF)*86400.,'conv_rain (mm/day)')

several_vars_zonalavg(ccrain_seasonal_avg.sel(season=JJA)*86400.,'ccrain (mm/day)',condensation_rain_seasonal_avg.sel(season=JJA)*86400.,'cond_rain (mm/day)',convection_rain_seasonal_avg.sel(season=JJA)*86400.,'conv_rain (mm/day)')


#http://stackoverflow.com/questions/6991471/computing-cross-correlation-function

#import cv2 as cv

# resultNp=np.zeros((tsurf_avg.shape[0],tsurf_avg.shape[1]))

# tsurf_cv=cv.fromarray(np.float32(tsurf_seasonal_avg.sel(season=JJA)))
# ccrain_cv=cv.fromarray(np.float32(ccrain_seasonal_avg.sel(season=JJA)))
# resultNp_cv=cv.fromarray(np.float32(resultNp))


# cv.MatchTemplate(tsurf_cv,ccrain_cv,resultNp_cv,cv.CV_TM_CCORR_NORMED)

# resultNp = np.asarray(resultNp_cv)

# print('Spatial correlation of ccrain and tsurf for JJA ='+str(resultNp.mean()))


# this is not working yet does not know fromarray
