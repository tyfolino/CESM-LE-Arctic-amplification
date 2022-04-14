# We will be using the adjusted cloud radiative effect method for determining cloud feedbacks
# All-sky and clear-sky feedbacks are made in Albedo/Q/T_feedbacks.py

# By: Ty Janoski
# Edited: 12.28.21

# import statements

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

#####################################################################################
##                                  CAM4                                           ##
#####################################################################################

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

if(m==7):
    inst_RF = xr.open_dataset('/dx02/janoski/cesm/inst_rf/CAM4_inst_4xCO2_RF_2yr_regridded.nc').roll(time=-181,roll_coords=True)
    end=51
else:
    inst_RF = xr.open_dataset('/dx02/janoski/cesm/inst_rf/CAM4_inst_4xCO2_RF_2yr_regridded.nc')
    end=101

dRFLW_TOA = np.array(-1 * (inst_RF.FLNTC - inst_RF.FLNT))
dRFSW_TOA = np.array(inst_RF.FSNTC - inst_RF.FSNT)
dRFLW_SFC = np.array(-1 * (inst_RF.FLNSC - inst_RF.FLNS))
dRFSW_SFC = np.array(inst_RF.FSNSC - inst_RF.FSNS)

#### note: convention at TOA is that positive = downward radiative flux (so negative values = cooling)
# print('Doing CAM4 TOA cloud feedbacks...')

# for e in range(1,end,1):
#     if(e%5==0):
#         print(e)
#     # read in T feedbacks
#     lapse = xr.open_dataset('/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_lapse_tropo_TOA.nc',
#                             use_cftime=True)
#     planck = xr.open_dataset('/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_Planck_tropo_TOA.nc',
#                              use_cftime=True)
#     dT = (lapse.FLNTC - lapse.FLNT) + (planck.FLNTC - planck.FLNT)

#     # Q feedbacks
#     Q = xr.open_dataset('/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(f"{m:02d}")+'.'
#                         +str(f"{e:02d}")+'.h1_q_tropo.nc',
#                         use_cftime=True)
#     dQ_LW = Q.FLNTC - Q.FLNT
#     dQ_SW = Q.FSNTC - Q.FSNT

#     # albedo
#     α = xr.open_dataset('/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(f"{m:02d}")+'.'+
#                         str(f"{e:02d}")+'.h1_albedo.nc',
#                         use_cftime=True)
#     dα = (α.FSNTC - α.FSNT).fillna(0)
    
#     # calc change in CRE
#     FLNT = read_in('ctrl',m,e,'FLNT').FLNT
#     FLNTC = read_in('ctrl',m,e,'FLNTC').FLNTC
#     FSNT = read_in('ctrl',m,e,'FSNT').FSNT
#     FSNTC = read_in('ctrl',m,e,'FSNTC').FSNTC

#     CRE_ctrl_SW = (FSNT - FSNTC) # it's SW - LW because LW is positive upwards (we want it positive downwards)
#     CRE_ctrl_LW = -1*(FLNT - FLNTC)

#     FLNT = read_in('4xCO2',m,e,'FLNT').FLNT
#     FLNTC = read_in('4xCO2',m,e,'FLNTC').FLNTC
#     FSNT = read_in('4xCO2',m,e,'FSNT').FSNT
#     FSNTC = read_in('4xCO2',m,e,'FSNTC').FSNTC

#     CRE_exp_SW = (FSNT - FSNTC)
#     CRE_exp_LW = -1*(FLNT - FLNTC)

#     dCRE_SW = CRE_exp_SW - CRE_ctrl_SW
#     dCRE_LW = CRE_exp_LW - CRE_ctrl_LW

#     # basis for adjustment method is to adjust the change in CRE for environmental masking of other radiative perturbations
#     # albedo terms has some nans, so just turn them into 0s
#     cloud_SW = dCRE_SW + dQ_SW + dα + dRFSW_TOA
#     cloud_LW = dCRE_LW + dQ_LW + dT + dRFLW_TOA
    
#     out_SW = xr.merge([cloud_SW.rename('cloud_SW'),dCRE_SW.rename('dCRE_SW'),dQ_SW.rename('dQ_SW'),dα.rename('dα')])
#     out_LW = xr.merge([cloud_LW.rename('cloud_LW'),dCRE_LW.rename('dCRE_LW'),dQ_LW.rename('dQ_LW'),dT.rename('dT')])
    
#     pathout_SW = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
#         f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_cloud_SW_TOA.nc'
#     pathout_LW = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
#         f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_cloud_LW_TOA.nc'
    
#     out_SW.to_netcdf(pathout_SW)
#     out_LW.to_netcdf(pathout_LW)
    


## note: convention at sfc is that positive = upwards radiative flux (so negative values = cooling)
print('Doing CAM4 sfc cloud feedbacks...')
for e in range(1,end,1):
    if(e%5==0):
        print(e)
    # read in T feedbacks
    sfc_T = xr.open_dataset('/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_Ts_tropo_sfc.nc',use_cftime=True)
    atm_T = xr.open_dataset('/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_Ta_tropo_sfc.nc',use_cftime=True)

    dT = (sfc_T.FLNSC - sfc_T.FLNS) + (atm_T.FLNSC - atm_T.FLNS)
    
#   Q feedbacks
    Q = xr.open_dataset('/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_q_tropo.nc',use_cftime=True)
    dQ_LW = Q.FLNSC - Q.FLNS
    dQ_SW = Q.FSNSC - Q.FSNS

    # albedo
    α = xr.open_dataset('/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
            f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_albedo.nc',use_cftime=True)
    dα = (α.FSNSC - α.FSNS).fillna(0)
    
    # calc change in CRE
    FLNS = read_in('ctrl',m,e,'FLNS').FLNS
    FLNSC = read_in('ctrl',m,e,'FLNSC').FLNSC
    FSNS = read_in('ctrl',m,e,'FSNS').FSNS
    FSNSC = read_in('ctrl',m,e,'FSNSC').FSNSC

    CRE_ctrl_SW = (FSNS - FSNSC) # it's SW - LW because LW is positive upwards (we want it positive downwards)
    CRE_ctrl_LW = -1*(FLNS - FLNSC)

    FLNS = read_in('4xCO2',m,e,'FLNS').FLNS
    FLNSC = read_in('4xCO2',m,e,'FLNSC').FLNSC
    FSNS = read_in('4xCO2',m,e,'FSNS').FSNS
    FSNSC = read_in('4xCO2',m,e,'FSNSC').FSNSC

    CRE_exp_SW = (FSNS - FSNSC)
    CRE_exp_LW = -1*(FLNS - FLNSC)

    dCRE_SW = CRE_exp_SW - CRE_ctrl_SW
    dCRE_LW = CRE_exp_LW - CRE_ctrl_LW

    # basis for adjustment method is to adjust the change in CRE for environmental masking of other radiative perturbations
    # albedo terms has some nans, so just turn them into 0s
    cloud_SW = dCRE_SW + dQ_SW + dα +  dRFSW_SFC
    cloud_LW = dCRE_LW + dQ_LW + dT + dRFLW_SFC
    
    out_SW = xr.merge([cloud_SW.rename('cloud_SW'),dCRE_SW.rename('dCRE_SW'),dQ_SW.rename('dQ_SW'),dα.rename('dα')])
    out_LW = xr.merge([cloud_LW.rename('cloud_LW'),dCRE_LW.rename('dCRE_LW'),dQ_LW.rename('dQ_LW'),dT.rename('dT')])
    
    pathout_SW = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_cloud_SW_sfc.nc'
    pathout_LW = '/dx02/janoski/cesm/vert_int_feedbacks/b40.1850.cam5-lens.'+str(
        f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_cloud_LW_sfc.nc'
    
    out_SW.to_netcdf(pathout_SW)
    out_LW.to_netcdf(pathout_LW)