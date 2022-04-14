# In this script, we will use the temporally regridded CAM5 kernels
# to compute the radiataive perturbations at the surface and TOA
# from changes in surface albedo.

# By: Ty Janoski
# Updated: 12.28.21


# import statements
import numpy as np
import xarray as xr
import sys

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

# get kernels and repeat them for additional year
if (m==1):
    kern = xr.open_dataset('/dx02/janoski/cesm/daily_kernels/alb.kernel.nc',use_cftime=True).transpose('time','lat','lon')
    end=101
elif(m==7):
    # if initialization month is July, make sure we move July 1st to the first spot on the time axis
    kern = xr.open_dataset('/dx02/janoski/cesm/daily_kernels/alb.kernel.nc',use_cftime=True).transpose('time','lat','lon').roll(
        time=-181,roll_coords=True)
    end=51
FSNS = np.tile(kern.FSNS,(2,1,1))
FSNSC = np.tile(kern.FSNSC,(2,1,1))
FSNT = np.tile(kern.FSNT,(2,1,1))
FSNTC = np.tile(kern.FSNTC,(2,1,1))

for e in range(1,end,1):
    print(e)
    # get control albedo
    net = read_in('ctrl',m,e,'FSNS').FSNS
    down = read_in('ctrl',m,e,'FSDS').FSDS

    α_ctrl = (down - net)/(down.where(down>0))

    # same for 4xCO2
    net = read_in('4xCO2',m,e,'FSNS').FSNS
    down = read_in('4xCO2',m,e,'FSDS').FSDS

    α_exp = (down - net)/(down.where(down>0))

    # find difference, multiply by 100 to make percent
    diff = (α_exp - α_ctrl)*100
    
    # we want 0's, not NaN's, in polar night
    diff = diff.fillna(0)

    # calculate radiative perturbation
    ΔFSNS = (diff * FSNS).rename('FSNS')
    ΔFSNSC = (diff * FSNSC).rename('FSNSC')
    ΔFSNT = (diff * FSNT).rename('FSNT')
    ΔFSNTC = (diff * FSNTC).rename('FSNTC')

    # make a dataset to save out
    out = xr.merge([ΔFSNS,ΔFSNSC,ΔFSNT,ΔFSNTC])

    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
            f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_albedo.nc'
    
    out.to_netcdf(pathout)