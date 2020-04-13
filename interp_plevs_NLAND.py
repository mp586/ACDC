from isca.util import interpolate_output


infiles = ['LandworldLakes_acdc_atmos_mean_5_to_50',    'Northland_dark_acdc_atmos_mean_5_to_50',
'Landworld_acdc_atmos_mean_5_to_50',         'Northland_dry_acdc_atmos_mean_5_to_50',
'Northland_bright_acdc_atmos_mean_5_to_50',  'Northland_empty_acdc_atmos_mean_5_to_50']

from isca.util import interpolate_output
for infile_in in infiles:
    infile = '/scratch/mp586/Isca_DATA/Northland/'+infile_in+'.nc'
    outfile = '/scratch/mp586/Isca_DATA/Northland/'+infile_in+'_interp.nc'
    interpolate_output(infile, outfile, p_levs='INPUT', var_names=['slp', 'height', 'omega', 'ucomp', 'vcomp', 'temp','rh','sphum'])






