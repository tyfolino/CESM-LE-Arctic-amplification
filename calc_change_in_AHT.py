# Purpose: Calculate the moist static energy (MSE) transport for output from the CESM model - both to verify control climate
# and to find difference after 4xCO2

# By: Ty Janoski
# Updated: 12.30.21

# import statments
import numpy as np
import xarray as xr

def read_in(exp,mon,ens,var):
    """
    Use xarray to read in a netCDF file.

    Keyword arguments:
    exp -- CO2 scenario
    mon -- starting month in which CO2 is altered
    ens -- ensemble number
    var -- model output variable
    """
    filein = '/dx02/janoski/cesm/ctrl_4xCO2_transports/b40.1850.cam5-lens.'+exp+'.'+str(
        f"{mon:02d}")+'.'+str(f"{ens:02d}")+'.h1_'+var+'.nc'
    return(xr.open_dataset(filein,chunks=None))

# Before getting into the 4xCO$_2$ simulations, we will see if the meridional MSE transport is reasonable in the control climate simulations.

# Note: because of the computational intensity of calculating MSE, they are calculated separately in the script AHT_explicit.py and AHT_as_residual.py

# Here we will be using these premade files

# get starting month
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()
    
if(m==1):
    end=101
elif(m==7):
     end=51
        
for e in range(1,end,1):
    if(e%5==0):
        print(e)
    ctrl = read_in('ctrl',m,e,'Fwall_70N_expl')
    exp = read_in('4xCO2',m,e,'Fwall_70N_expl')
    diff = exp - ctrl
    
    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_AHT_70N_expl.nc'
    diff.to_netcdf(pathout)
    
for e in range(1,end,1):
    if(e%5==0):
        print(e)
    ctrl = read_in('ctrl',m,e,'Fconv_resi')
    exp = read_in('4xCO2',m,e,'Fconv_resi')
    diff = exp - ctrl
    
    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_AHT_resi.nc'
    diff.to_netcdf(pathout)
    
for e in range(1,end,1):
    if(e%5==0):
        print(e)
    ctrl = read_in('ctrl',m,e,'Wconv')
    exp = read_in('4xCO2',m,e,'Wconv')
    diff = exp - ctrl
    
    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_Wconv.nc'
    diff.to_netcdf(pathout)
    
for e in range(1,end,1):
    if(e%5==0):
        print(e)
    ctrl = read_in('ctrl',m,e,'Sconv')
    exp = read_in('4xCO2',m,e,'Sconv')
    diff = exp - ctrl
    
    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_Sconv.nc'
    diff.to_netcdf(pathout)