# Purpose: this script will calculate the difference in surface air temperature

# By: Ty Janoski
# Edited: 02.02.22

# import statements

import numpy as np
import xarray as xr

# create function for reading in cesm-LE output
# note: each ensemble member starts on a different year
# please keep this in mind when combining datasets
def read_in(exp,mon,ens,var):
    """
    Use xarray to read in a netCDF file.

    Keyword arguments:
    exp -- CO2 scenario
    mon -- starting month in which CO2 is altered
    ens -- ensemble number
    var -- model output variable
    """
    filein = '/dx01/janoski/cesm/output/b40.1850.cam5-lens.'+exp+'.'+str(
         f"{mon:02d}")+'.'+str(f"{ens:02d}")+'.h1_'+var+'.nc'
    return(xr.open_dataset(filein,chunks=None))

# get starting month
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()
    
if (m==7):
    end=51
else:
    end=101

# iterate through each ensemble member
for e in range(1,end,1):
    
    print(e)

    # get the difference in surface skin temperature as well
    ctrl = read_in('ctrl',m,e,'TREFHT')
    exp = read_in('4xCO2',m,e,'TREFHT')
    
    dSAT = exp.TREFHT - ctrl.TREFHT
    del exp, ctrl
    
    # combine all data arrays to make dataset
    # save it out as netCDF file
    out = dSAT.rename('dSAT')
    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_dSAT.nc'
    
    out.to_netcdf(pathout)
    
# iterate through each ensemble member
for e in range(1,end,1):
    
    print(e)

    # get the difference in surface skin temperature as well
    ctrl = read_in('ctrl',m,e,'TS')
    exp = read_in('4xCO2',m,e,'TS')
    
    dTS = exp.TS - ctrl.TS
    del exp, ctrl
    
    # combine all data arrays to make dataset
    # save it out as netCDF file
    out = dTS.rename('TS')
    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_dTS.nc'
    
    out.to_netcdf(pathout)
