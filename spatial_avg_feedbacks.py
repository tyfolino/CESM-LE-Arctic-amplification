# import statements
import xarray as xr
import numpy as np
import dask as da
import matplotlib.pyplot as plt
import xesmf as xe
from dask.diagnostics import ProgressBar
import sys

# read in cell area file for taking spatial averages
areacella = xr.open_dataarray('/dx01/janoski/cesm/output/gridarea.nc')
area=areacella
# create function for taking spatial averages, while weighting for latitude
def spav(ds_in, areacella=areacella, lat_bound_s = -91, lat_bound_n = 91):
    """
    Use xarray/numpy to calculate spatial average while weighting for latitude.
    
    Keyword arguments:
    ds_in -- Dataset or DataArray to take the average of with coords lat and lon
    areacella -- Dataset or DataArray containing the grid cell area with coords lat and lon
    lat_bound_s -- float, Southern boundary of area to average
    lat_bound_n -- float, Northern boundary of area to average
    """
    ds_in = ds_in.sel(lat=slice(lat_bound_s,lat_bound_n))
    areacella = areacella.sel(lat=slice(lat_bound_s,lat_bound_n))
    out = (ds_in*(areacella/areacella.sum(dim=['lat','lon']))).sum(dim=['lat','lon'])
    return(out)

def preprocess(ds):
    dsnew = ds.copy()
    dsnew['time'] = np.arange(1,731,1)
    return dsnew
# ====================================================================================================

# read in all 100 ensemble members lazily, overwriting the time axis in order to concatenate
# along new ensemble dimension
# get starting month
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()

lead = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
    f"{m:02d}")+'.*.h1_'

AHT = xr.open_mfdataset(lead+'AHT_resi.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').AHT
dEdt = xr.open_mfdataset(lead+'AHT_resi.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').dEdt
alb_toa = xr.open_mfdataset(lead+'albedo.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FSNT
alb_sfc = xr.open_mfdataset(lead+'albedo.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FSNS
clw_toa = xr.open_mfdataset(lead+'cloud_LW_TOA.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').cloud_LW
clw_sfc = xr.open_mfdataset(lead+'cloud_LW_sfc.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').cloud_LW
csw_toa = xr.open_mfdataset(lead+'cloud_SW_TOA.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').cloud_SW
csw_sfc = xr.open_mfdataset(lead+'cloud_SW_sfc.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').cloud_SW
dTS = xr.open_mfdataset(lead+'dTS.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').TS
dSAT = xr.open_mfdataset(lead+'dSAT.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').dSAT
lapse = xr.open_mfdataset(lead+'lapse_tropo_TOA.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FLNT
planck = xr.open_mfdataset(lead+'Planck_tropo_TOA.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FLNT
qlw_toa = xr.open_mfdataset(lead+'q_tropo.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FLNT
qsw_toa = xr.open_mfdataset(lead+'q_tropo.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FSNT
qlw_sfc = xr.open_mfdataset(lead+'q_tropo.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FLNS
qsw_sfc = xr.open_mfdataset(lead+'q_tropo.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FSNS
Sconv = xr.open_mfdataset(lead+'Sconv.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').Sconv
Fsfc = xr.open_mfdataset(lead+'Fsfc.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').dFsfc
dLHFLX = xr.open_mfdataset(lead+'Fsfc.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').dLHFLX
dSHFLX = xr.open_mfdataset(lead+'Fsfc.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').dSHFLX
Ta = xr.open_mfdataset(lead+'Ta_tropo_sfc.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FLNS
Ts = xr.open_mfdataset(lead+'Ts_tropo_sfc.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').FLNS
T_strato = xr.open_mfdataset(lead+'T_strato_TOA.nc',preprocess=preprocess,combine='nested',
                            concat_dim='ens').FLNT
qlw_strato = xr.open_mfdataset(lead+'q_strato.nc',preprocess=preprocess,combine='nested',
                            concat_dim='ens').FLNT
qsw_strato = xr.open_mfdataset(lead+'q_strato.nc',preprocess=preprocess,combine='nested',
                            concat_dim='ens').FSNT
Wconv = xr.open_mfdataset(lead+'Wconv.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens').Wconv
dSIC = xr.open_mfdataset(lead+'dSIC.nc',preprocess=preprocess,combine='nested',
                         concat_dim='ens').dSIC

# read in CO2 RF
fin = xr.open_dataset('/dx02/janoski/cesm/inst_rf/CAM4_inst_4xCO2_RF_2yr_regridded.nc',
                      use_cftime=True)

if(m==1):
    CO2SW_TOA = fin.FSNT
    CO2LW_TOA = -1 * fin.FLNT
    CO2SW_sfc = fin.FSNS
    CO2LW_sfc = -1 * fin.FLNS
elif(m==7):
    CO2SW_TOA = fin.FSNT.roll(time=-181,roll_coords=True)
    CO2LW_TOA = -1 * fin.FLNT.roll(time=-181,roll_coords=True)
    CO2SW_sfc = fin.FSNS.roll(time=-181,roll_coords=True)
    CO2LW_sfc = fin.FLNS.roll(time=-181,roll_coords=True)
CO2_TOA = CO2SW_TOA + CO2LW_TOA
CO2_sfc = CO2SW_SFC + CO2LW_SFC
CO2_TOA['time'] = Wconv.time.data
CO2_sfc['time'] = Wconv.time.data

# we also need the land fraction
lnd = xr.open_dataset('/dx01/janoski/cesm/output/CESM-LE_landfrac.nc').landfrac
lnd = lnd.fillna(0)
lnd['lat'] = CO2_TOA.lat.data
lnd['lon'] = CO2_sfc.lon.data

# ==========================================================================================================

# begin spatial averaging

if(m==1):
    pathout = '/dx02/janoski/cesm/spat_avg_feedbacks/all/b40.1850.cam5-lens.01.glb_'
elif(m==7):
    pathout = '/dx02/janoski/cesm/spat_avg_feedbacks/all/b40.1850.cam5-lens.07.glb_'
    
print(pathout)

to_avg = [AHT,dEdt,alb_toa,alb_sfc,clw_toa,clw_sfc,csw_toa,csw_sfc,dTS,dSAT,
         lapse,planck,qlw_toa,qsw_toa,qlw_sfc,qsw_sfc,Sconv,Fsfc,dLHFLX,dSHFLX,
         Ta,Ts,T_strato,qlw_strato,qsw_strato,Wconv,CO2_TOA]
out_names = ['AHT.nc','dEdt.nc','alb_TOA.nc','alb_sfc.nc','cloud_LW_TOA.nc','cloud_LW_sfc.nc',
             'cloud_SW_TOA.nc','cloud_SW_sfc.nc','dTS.nc','dSAT.nc','lapse_tropo_TOA.nc','planck_tropo_TOA.nc','q_LW_tropo_TOA.nc',
             'q_SW_tropo_TOA.nc','q_LW_tropo_sfc.nc','q_SW_tropo_sfc.nc','Sconv.nc','Fsfc.nc','dLHFLX.nc','dSHFLX.nc',
             'Ta_tropo_sfc.nc','Ts_tropo_sfc.nc','T_strato_TOA.nc','q_LW_strato_TOA.nc','q_SW_strato_TOA.nc',
             'Wconv.nc','CO2_TOA.nc']

for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area).compute().to_netcdf(pathout+out_names[i])
        
pathout = pathout.replace('glb','arc')
for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area,lat_bound_s=70).compute().to_netcdf(pathout+out_names[i])
        
pathout = pathout.replace('arc','ant')
for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area,lat_bound_n=-70).compute().to_netcdf(pathout+out_names[i])
        
print(pathout)
        
pathout = pathout.replace('ant','glb').replace('all','land_only')
for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area.where(lnd>=0.5)).compute().to_netcdf(pathout+out_names[i])
        
print(pathout)
        
pathout = pathout.replace('glb','arc')
for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area.where(lnd>=0.5),lat_bound_s=70).compute().to_netcdf(pathout+out_names[i])
        
pathout = pathout.replace('arc','ant')
for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area.where(lnd>=0.5),lat_bound_n=-70).compute().to_netcdf(pathout+out_names[i])
        
print(pathout)
        
pathout = pathout.replace('ant','glb').replace('land_only','ocean_only')
for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area.where(lnd<0.5)).compute().to_netcdf(pathout+out_names[i])
        
print(pathout)
        
pathout = pathout.replace('glb','arc')
for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area.where(lnd<0.5),lat_bound_s=70).compute().to_netcdf(pathout+out_names[i])
        
print(pathout)

pathout = pathout.replace('arc','ant')
for i in range(len(to_avg)):
    with ProgressBar():
        spav(to_avg[i],area.where(lnd<0.5),lat_bound_n=-70).compute().to_netcdf(pathout+out_names[i])
        
print(pathout)

#==========================explicit AHT======================

# this portion is a little different because the Arctic AHT is calculated specifically at 70N
# so there is no spatial averaging needed.

Fwall = xr.open_mfdataset(lead+'AHT_70N_expl.nc',preprocess=preprocess,combine='nested',
                     concat_dim='ens')
AHT_expl = Fwall.MSE
Wconv_expl = Fwall.VQ
Sconv_expl = Fwall.VZ+Fwall.VT

if(m==1):
    pathout = '/dx02/janoski/cesm/spat_avg_feedbacks/all/b40.1850.cam5-lens.01.glb_'
elif(m==7):
    pathout = '/dx02/janoski/cesm/spat_avg_feedbacks/all/b40.1850.cam5-lens.07.glb_'
    
to_avg = [AHT_expl,Wconv_expl,Sconv_expl]
out_names = ['AHT_expl.nc','Wconv_expl.nc','Sconv_expl.nc']

for i in range(len(to_avg)):
    (to_avg[i]*0).to_netcdf(pathout+out_names[i])
    
pathout = pathout.replace('glb','arc')
for i in range(len(to_avg)):
    to_avg[i].to_netcdf(pathout+out_names[i])