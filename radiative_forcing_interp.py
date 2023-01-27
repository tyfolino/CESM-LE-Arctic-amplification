# Note: this script replaces the time_interp.ipynb

# By: Ty Janoski
# Updated: 08.18.22

# import statments
import numpy as np
import xarray as xr
import nc_time_axis
from cftime import DatetimeNoLeap
import xesmf as xe

# read in the CAM4 instantaneous 4xCO2 RF file
orig = xr.open_dataset('/dx02/janoski/cesm/inst_rf/forcing.TOA.4xCO2.inst.nc',use_cftime=True)

# make dates we will interpolate to
dates = xr.cftime_range(start="0001-01-01 12:00:00",end="0002-12-31 12:00:00",freq='D',calendar='noleap')

# interpolate to noon of each day
new = orig.interp(time=dates,method='linear')

# now, let's regrid it to match CESM
fin = xr.open_dataset('/dx01/janoski/cesm/output/b40.1850.cam5-lens.4xCO2.01.01.h1_PS.nc')
lat = fin.lat
lon = fin.lon
del fin

# make dataarray to regrid too
ds_out = xr.Dataset({'lat': (['lat'], lat.data),
                     'lon': (['lon'], lon.data),
                    }
                   )

# make regridder
regridder = xe.Regridder(new.FLNS, ds_out, 'bilinear',periodic=True)

# regrid!
new_regrid = regridder(new)

# now just save it out
new_regrid.to_netcdf('/dx02/janoski/cesm/inst_rf/forcing.TOA.4xCO2.inst.regridded.nc')