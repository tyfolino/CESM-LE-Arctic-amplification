# import statements

import xarray as xr
import numpy as np

def preprocess(ds):
    dsnew = ds.copy()
    dsnew['time'] = np.arange(1,731,1)
    return dsnew

# get starting month
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()
    
areacella = xr.open_dataarray('/dx01/janoski/cesm/output/gridarea.nc')

# read in ctrl ICEFRAC
sic_ctrl = xr.open_mfdataset('/dx01/janoski/cesm/output/b40.1850.cam5-lens.ctrl.'+
                             str(f"{m:02d}")+'.*.h1_ICEFRAC.nc',preprocess=preprocess,
                             use_cftime=True,combine='nested',concat_dim='ens').ICEFRAC
# 4xCO2
sic_exp = xr.open_mfdataset('/dx01/janoski/cesm/output/b40.1850.cam5-lens.4xCO2.'+
                             str(f"{m:02d}")+'.*.h1_ICEFRAC.nc',preprocess=preprocess,
                             use_cftime=True,combine='nested',concat_dim='ens').ICEFRAC
# take difference
sic_diff = sic_exp - sic_ctrl

# multiply ice area fraction by cell area to get sea ice area (limit lat>0 to make it Arctic)
sia_diff_abs = (sic_diff * areacella).where(sic_diff.lat>0).sum(dim=['lat','lon'])

# same for control
sia_ctrl = (sic_ctrl * areacella).where(sic_ctrl.lat>0).sum(dim=['lat','lon'])

# also get difference in sea ice area as a percent
sia_diff_perc = (sia_diff_abs/sia_ctrl.where(sia_ctrl>0)) * 100 

# save out
pathout = '/dx02/janoski/cesm/sea_ice/ctrl/b40.1850.cam5-lens.'+str(
    f"{m:02d}")+'.SIC.nc'
sic_ctrl.to_netcdf(pathout)

pathout = '/dx02/janoski/cesm/sea_ice/diff/b40.1850.cam5-lens.'+str(
    f"{m:02d}")+'.dSIC.nc'
sic_diff.to_netcdf(pathout)

pathout = '/dx02/janoski/cesm/sea_ice/diff/b40.1850.cam5-lens.'+str(
    f"{m:02d}")+'.dSIA_abs.nc'
sia_diff_abs.to_netcdf(pathout)

pathout = '/dx02/janoski/cesm/sea_ice/ctrl/b40.1850.cam5-lens.'+str(
    f"{m:02d}")+'.SIA.nc'
sia_ctrl.to_netcdf(pathout)

pathout = '/dx02/janoski/cesm/sea_ice/diff/b40.1850.cam5-lens.'+str(
    f"{m:02d}")+'.dSIA_perc.nc'
# sia_diff_perc.to_netcdf(pathout)

#### SAT

# read in ctrl ICEFRAC
sat_ctrl = xr.open_mfdataset('/dx01/janoski/cesm/output/b40.1850.cam5-lens.ctrl.'+
                             str(f"{m:02d}")+'.*.h1_TREFHT.nc',
                              preprocess=preprocess,use_cftime=True,combine='nested',concat_dim='ens').TREFHT
# 4xCO2
sat_exp = xr.open_mfdataset('/dx01/janoski/cesm/output/b40.1850.cam5-lens.4xCO2.'+
                             str(f"{m:02d}")+'.*.h1_TREFHT.nc',
                              preprocess=preprocess,use_cftime=True,combine='nested',concat_dim='ens').TREFHT
# take difference
sat_diff = sat_exp - sat_ctrl

pathout = '/dx02/janoski/cesm/sea_ice/diff/b40.1850.cam5-lens.'+str(
    f"{m:02d}")+'.dSAT.nc'
sat_diff.to_netcdf(pathout)