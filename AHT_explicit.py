#!/usr/bin/env python
# coding: utf-8

# Purpose: Calculate the moist static energy (MSE) transport for output from the CESM model - both to verify control climate
# and to find difference after 4xCO2

# By: Ty Janoski
# Updated: 03.30.2023

# for Arctic/Antarctic, change latitude, sign of terms, and output file name


# import statments
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

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

# few constants, a is radius of Earth, Cp is specific heat of air, L is enthalpy of vaporization, g is gravity
a = 6.371e6
Cp = 1003
L = 2.5e6
g = 9.81

# for some reason, hyai and hybi are not in the July initalization netCDFs,
# so let's grab it from one of the control simulations.
hyai = xr.open_dataset('/dx01/janoski/cesm/output/b40.1850.cam5-lens.ctrl.01.01.h1_covar.nc',
                     drop_variables=['hyam','hybm','P0','time_bnds','PS','VQ','VZ','VT']).hyai
hybi = xr.open_dataset('/dx01/janoski/cesm/output/b40.1850.cam5-lens.ctrl.01.01.h1_covar.nc',
                     drop_variables=['hyam','hybm','P0','time_bnds','PS','VQ','VZ','VT']).hybi

# get starting month
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()
    
if(m==1):
    end=101
elif(m==7):
     end=51

# iterate through control and 4xCO2 simulations
for s in ['ctrl','4xCO2']:
    # read in simulations (one ensemble member at a time)
    for e in range(1,end,1):
        print(e)
        if(m==1):
            ds = read_in(s,m,e,'covar')
        elif(m==7):
            VQ = read_in(s,m,e,'VQ')
            VT = read_in(s,m,e,'VT')
            VZ = read_in(s,m,e,'VZ')
            ds = xr.merge([VQ,VT,VZ,hyai,hybi])

        # we are going to need the closest latitude to our 70 degree N band
        lat = float(ds.lat.sel(lat=-70,method='nearest'))
        
        # to get units of W/m^2, we divide by the area of the Arctic (north of our lat)
        area = 2*np.pi*a**2 * (np.sin(np.deg2rad(-70)) - np.sin(np.deg2rad(-90)))
        if(e==1):
            print(area)
        
        # let's slice ds to save space
        ds = ds.sel(lat=lat)

        # calculate the distance between points on this latitude circle (288 lon points in CESM)
        # dx = circumference of earth at this latitude / number of longitude points
        dx = 2 * np.pi * a * np.cos(np.deg2rad(lat))/len(ds.lon)

        # calculate difference in pressure levels at interfaces
        # first part converts from hybrid-sigma to pressure, second part does the vertical gradient
        # third part switch time and vertical coordinates order for compatibility
        dp = np.array((ds.hyai * ds.P0 + ds.hybi * ds.PS).diff(dim='ilev').transpose('time','ilev','lon'))

        # get components of MSE transport
        VQ = -1 * (((ds.VQ * L * dp).sum(dim='lev') * dx).sum(dim='lon') / g).rename('VQ')
        VT = -1 * (((ds.VT * Cp * dp).sum(dim='lev') * dx).sum(dim='lon') / g).rename('VT')
        VZ = -1 * (((ds.VZ * dp).sum(dim='lev') * dx).sum(dim='lon') / g).rename('VZ')

        # get components of MSE transport
        # the below method is only if you want to avoid integrating over lon
    #     VQ = (((ds.VQ * L * dp).sum(dim='lev')*dx) / g).rename('VQ')
    #     VT = (((ds.VT * Cp * dp).sum(dim='lev')*dx) / g).rename('VT')
    #     VZ = (((ds.VZ * dp).sum(dim='lev')*dx) / g).rename('VZ')


        # MSE will be the sum of the three (excluding kinetic energy which is orders of magnitude smaller than the others)
        MSE = (VQ + VT + VZ).rename('MSE')

        ds_out = xr.merge([VQ,VT,VZ,MSE])/area

        filename = 'b40.1850.cam5-lens.'+s+'.'+str(f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_Fwall_70S_expl.nc'
        ds_out.to_netcdf('/dx02/janoski/cesm/ctrl_4xCO2_transports/'+filename)