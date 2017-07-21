from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
from scipy.stats import *


def pattern_corr_2d(arr1,arr2):

    arr1=np.asarray(arr1)
    arr2=np.asarray(arr2)
    arr1_1d=arr1.flatten() #turn MxN array into 1X(M*N) vector
    arr2_1d=arr2.flatten()

    r,pval=pearsonr(arr1_1d,arr2_1d)

    return(r,pval)
