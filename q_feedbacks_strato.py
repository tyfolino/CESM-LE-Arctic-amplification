# Purpose: this script will calculate the difference in the natural log of specific humidity,
# multiply by the CAM5 water vapor kernels, and perform a vertical integration from the surface to the tropopause.
# note: kernels were already normalized by log(q) factor

# By: Ty Janoski
# Edited: 01.18.22

# import statements

import numpy as np
import xarray as xr
import ngl

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

# because CESM output uses hybrid sigma vertical coordinates, the upper levels are isobaric surfaces
upper = np.array([  3.64346569,   7.59481965,  14.35663225,  24.61222   ,
        38.26829977,  54.59547974,  72.01245055,  87.82123029,
       103.31712663, 121.54724076, 142.99403876, 168.22507977])
# for lower level, we will go by 25 hPa
plevs = np.concatenate([upper,np.arange(200,1001,25)])

# get midpoints (= interface levels) for layer thickness later
mids = (plevs[1:] + plevs[:-1])/2

# make a list of kernel variables
kern_names = ['FLNS','FLNSC','FLNT','FLNTC','FSNS','FSNSC','FSNT','FSNTC']

# get starting month
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()
    
if(m==1):
    end=101
elif(m==7):
     end=51

# iterate through each ensemble member
for e in range(1,end,1):
    
    print(e)
    
    # read in kernels, control, and 4xCO2 output
    # grab a few variable attributes for later
    ctrl = read_in('ctrl',m,e,'Q')
    trop_c = read_in('ctrl',m,e,'TROP_P')
    exp = read_in('4xCO2',m,e,'Q')
    trop_e = read_in('4xCO2',m,e,'TROP_P')
    ps = exp.PS
    time = ps.time
    lat = ps.lat
    lon = ps.lon

    # for ctrl/exp, regrid the log of Q (spec hum) to p levels
    # pressure levels below surface become nan
    ctrl = ngl.vinth2p(np.log(ctrl.Q),ctrl.hyam,ctrl.hybm,plevs,ctrl.PS,1,ctrl.P0/100,1,kxtrp=False)
    ctrl[(ctrl==1e30) | (np.expand_dims(plevs,axis=[0,2,3]) > np.expand_dims(trop_c.TROP_P/100,axis=1))] = np.nan
    exp = ngl.vinth2p(np.log(exp.Q),exp.hyam,exp.hybm,plevs,exp.PS,1,exp.P0/100,1,kxtrp=False)
    exp[(exp==1e30) | (np.expand_dims(plevs,axis=[0,2,3]) > np.expand_dims(trop_e.TROP_P/100,axis=1))] = np.nan
    
    # take difference
    diff = exp - ctrl
    del exp, ctrl, trop_c, trop_e

    # make array of layer thicknesses for vertical integration - ps at bottom, 0 at top, difference in interface levels in between
    dp = np.diff(np.concatenate([np.zeros((730,1,192,288)),np.broadcast_to(np.expand_dims(mids,axis=[0,2,3]),(730,44,192,288)),
                np.expand_dims(np.array(ps),axis=1)/100],axis=1),axis=1)
    del ps
    
    # make blank list to store output
    rad_pert = []
    
    for k in kern_names[:]:
        
        # if starting in Jul, shift forward so that day 1 = jul 1st
        if(m==1):
            kern = np.load('/dx02/janoski/cesm/daily_kernels/q.kernel.'+k+'.npy')
        if(m==7):
            kern = np.roll(np.load('/dx02/janoski/cesm/daily_kernels/q.kernel.'+k+'.npy'),
                           -181,axis=0)
            
        # because of way longwave is defined, we need longwave kernels to be opposite sign
        if(k in ['FLNS','FLNSC','FLNT','FLNTC']):
            product = np.nansum(kern * dp * diff,axis=1) * -1
        else:
            product = np.nansum(kern * dp * diff,axis=1)
        
        
        # make our output into a data array to save out as netcdf file
        da = xr.DataArray(data=product,dims=['time','lat','lon'],coords=dict(
        lon=(["lon"], lon.data),
        lat=(["lat"], lat.data),
        time=time.data,
            ),
                attrs=dict(
                units="W/m^2",
            ),
        ).rename(k)
        
        rad_pert.append(da)
        del kern, product
    
    del dp, diff
    
    # combine all data arrays to make dataset
    # save it out as netCDF file
    out = xr.merge(rad_pert)
    out['lat'].attrs = lat.attrs
    out['lon'].attrs = lon.attrs
    out['time'].attrs = time.attrs
    pathout = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_q_strato.nc'
    
    out.to_netcdf(pathout)
