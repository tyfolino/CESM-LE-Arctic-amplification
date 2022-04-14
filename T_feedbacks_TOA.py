# Purpose: this script will calculate the difference in air temperature, surface skin temperature, then
# use the CAM5 kernels to produce the radiative perturbations at the TOA from the
# lapse rate and Planck feedbacks.

# Kernel work flow: kernels are produced in regrid_kernels_z&time.ipynb

# By: Ty Janoski
# Edited: 12.28.21

# import statements

import numpy as np
import xarray as xr
import ngl
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
    return(xr.open_dataset(filein,chunks=None,use_cftime=True))

# get starting month
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()

# because CESM output uses hybrid sigma vertical coordinates, the upper levels are isobaric surfaces
upper = np.array([  3.64346569,   7.59481965,  14.35663225,  24.61222   ,
        38.26829977,  54.59547974,  72.01245055,  87.82123029,
       103.31712663, 121.54724076, 142.99403876, 168.22507977])
# for lower level, we will go by 25 hPa
plevs = np.concatenate([upper,np.arange(200,1001,25)])

# get midpoints (= interface levels) for layer thickness later
mids = (plevs[1:] + plevs[:-1])/2

# make a list of kernel variables
# note: TOA and surface fluxes are partitioned slightly differently
# The following is for the TOA fluxes
kern_names = ['FLNT','FLNTC']

print('Reading in kernels...')

# The following step is to shift the kernels so that they start on Jul 1 if these are July init. simulations
if (m==7):
    tskerns = [np.roll(np.load('/dx02/janoski/cesm/daily_kernels/ts.kernel.'+k+'.npy'),-181,axis=0) for k in kern_names]
    tkerns = [np.roll(np.load('/dx02/janoski/cesm/daily_kernels/t.kernel.'+k+'.npy'),-181,axis=0) for k in kern_names]
    end=51
else:
    tskerns = [np.load('/dx02/janoski/cesm/daily_kernels/ts.kernel.'+k+'.npy') for k in kern_names]
    tkerns = [np.load('/dx02/janoski/cesm/daily_kernels/t.kernel.'+k+'.npy') for k in kern_names]
    end=101

# iterate through each ensemble member
for e in range(1,end,1):
    print(e)
    
    # read in kernels, control, and 4xCO2 output
    # grab a few variable attributes for later
    ctrl = read_in('ctrl',m,e,'T')
    trop_c = read_in('ctrl',m,e,'TROP_P')
    exp = read_in('4xCO2',m,e,'T')
    trop_e = read_in('4xCO2',m,e,'TROP_P')
    ps = exp.PS
    time = ps.time
    lat = ps.lat
    lon = ps.lon

    # for ctrl/exp, regrid T (air temp) to p levels using PyNGL vinth2p function
    # pressure levels below surface become nan
    ctrl = ngl.vinth2p(ctrl.T,ctrl.hyam,ctrl.hybm,plevs,ctrl.PS,1,ctrl.P0/100,1,kxtrp=False)
    ctrl[(ctrl==1e30) | (np.expand_dims(plevs,axis=[0,2,3]) < np.expand_dims(trop_c.TROP_P/100,axis=1))] = np.nan
    exp = ngl.vinth2p(exp.T,exp.hyam,exp.hybm,plevs,exp.PS,1,exp.P0/100,1,kxtrp=False)
    exp[(exp==1e30) | (np.expand_dims(plevs,axis=[0,2,3]) < np.expand_dims(trop_e.TROP_P/100,axis=1))] = np.nan
    
    # take difference of 4xCO2 and ctrl air temp
    dT = exp - ctrl
    del exp, ctrl, trop_c, trop_e

    # get surface air temperature
    ctrl = read_in('ctrl',m,e,'TREFHT')
    exp = read_in('4xCO2',m,e,'TREFHT')
    dTS = np.array(exp.TREFHT - ctrl.TREFHT)
    
    del exp, ctrl
    
    # Make a 4D version of dTS for separating lapse rate and Planck feedbacks
    dTS_uniform = np.expand_dims(dTS,axis=1)
    
    # change in lapse is defined as dT - dTS_uniform
    dT_lapse = dT - dTS_uniform
    
    del dT
    
    # make array of layer thicknesses for vertical integration - ps at bottom, 0 at top, difference in interface levels in between
    dp = np.diff(np.concatenate([np.zeros((730,1,192,288)),np.broadcast_to(np.expand_dims(mids,axis=[0,2,3]),(730,44,192,288)),
                np.expand_dims(np.array(ps),axis=1)/100],axis=1),axis=1)
    del ps
    
    # make blank list to store output
    rad_pert_p = []
    rad_pert_l = []
    
    for k in range(len(kern_names)):
    
        tkern = tkerns[k]
        # put NaNs where kernel is below the surface
        tkern = np.where(np.isfinite(dT_lapse),tkern,np.nan)
        tskern = tskerns[k]
        
        planck = (np.nansum(tkern * dp * dTS_uniform,axis=1) + (dTS * tskern))  * -1
        lapse = np.nansum(tkern * dp * (dT_lapse),axis=1) * -1
        
        # make our output into a data array to save out as netcdf file
        da = xr.DataArray(data=planck,dims=['time','lat','lon'],coords=dict(
        lon=(["lon"], lon.data),
        lat=(["lat"], lat.data),
        time=time.data,
            ),
                attrs=dict(
                units="W/m^2",
            ),
        ).rename(kern_names[k])
        
        rad_pert_p.append(da)
        
        # make our output into a data array to save out as netcdf file
        da = xr.DataArray(data=lapse,dims=['time','lat','lon'],coords=dict(
        lon=(["lon"], lon.data),
        lat=(["lat"], lat.data),
        time=time.data,
            ),
                attrs=dict(
                units="W/m^2",
            ),
        ).rename(kern_names[k])
        
        rad_pert_l.append(da)
        del planck,lapse
    
    del dp, dTS, dTS_uniform, dT_lapse
    
    # combine all data arrays to make dataset
    # save it out as netCDF file
    out = xr.merge(rad_pert_p)
    out['lat'].attrs = lat.attrs
    out['lon'].attrs = lon.attrs
    out['time'].attrs = time.attrs

    
    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_Planck_tropo_TOA.nc'
    
    out.to_netcdf(pathout)
    
    out = xr.merge(rad_pert_l)
    out['lat'].attrs = lat.attrs
    out['lon'].attrs = lon.attrs
    out['time'].attrs = time.attrs

    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_lapse_tropo_TOA.nc'
    
    out.to_netcdf(pathout)