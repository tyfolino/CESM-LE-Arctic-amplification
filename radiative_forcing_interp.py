# Note: this script replaces the time_interp.ipynb

# By: Ty Janoski
# Updated: 11.09.21

# import statments
import numpy as np
import xarray as xr
import nc_time_axis
from cftime import DatetimeNoLeap
import xesmf as xe

# read in the CAM4 instantaneous 4xCO2 RF file
RFin = xr.open_dataset('/dx02/janoski/cesm/inst_rf/forcing.TOA.4xCO2.nc',use_cftime=True)

# make sure to set boundary conditions on either side so that the interpolation is periodic
pre = RFin.isel(time=-1)
pre['time'] = DatetimeNoLeap(0,12,30,11,30,0,0,has_year_zero=True)

post=RFin.isel(time=0)
post['time'] = DatetimeNoLeap(2,1,1,0,0,0,0,has_year_zero=True)

# concatenate before interpolation
CAM4_b4_interp = xr.concat([pre,RFin,post],dim='time')

# make dates we will interpolate to
dates = xr.cftime_range(start="0001-01-01 12:00:00",end="0001-12-31 12:00:00",freq='D',calendar='noleap')

# interpolate to noon of each day
CAM4 = CAM4_b4_interp.interp(time=dates,method='linear')

# we could stop here, but it would be nice to make it two years long to match our CESM output
CAM4_y1 = CAM4.copy()

CAM4_y2 = CAM4.copy()
CAM4_y2['time'] = xr.cftime_range(start="0002-01-01 12:00:00",end="0002-12-31 12:00:00",freq='D',calendar='noleap')

# concat one last time
CAM4 = xr.concat([CAM4_y1,CAM4_y2],dim='time')

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
regridder = xe.Regridder(CAM4.FLNS, ds_out, 'bilinear',periodic=True)

# regrid!
CAM4_regrid = regridder(CAM4)

# now just save it out
CAM4_regrid.to_netcdf('/dx02/janoski/cesm/inst_rf/CAM4_inst_4xCO2_RF_2yr_regridded.nc')