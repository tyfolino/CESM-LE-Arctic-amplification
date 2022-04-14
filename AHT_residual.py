# Purpose: Calculate the moist static energy (MSE) transport for output from the CESM model - both to verify control climate
# and to find difference after 4xCO2

# By: Ty Janoski
# Updated: 01.18.22

# import statments
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import scipy as sp
from cftime import DatetimeNoLeap
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
    return(xr.open_dataset(filein,chunks=None))

# because CESM output uses hybrid sigma vertical coordinates, the upper levels are isobaric surfaces
upper = np.array([  3.64346569,   7.59481965,  14.35663225,  24.61222   ,
        38.26829977,  54.59547974,  72.01245055,  87.82123029,
       103.31712663, 121.54724076, 142.99403876, 168.22507977])
# for lower level, we will go by 25 hPa
plevs = np.concatenate([upper,np.arange(200,1001,25)])

# get midpoints (= interface levels) for layer thickness later
mids = (plevs[1:] + plevs[:-1])/2

# get starting month
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()

# read in topography
topo = xr.open_dataset('/dx01/janoski/cesm/output/orog_fx_CCSM4_decadal1961_r0i0p0.nc').orog

if(m==1):
    end=101
elif(m==7):
    end=51

# do this for ctrl annd 4xCO2 simulations
for s in ['ctrl','4xCO2']:
    # iterate through each ensemble member
    for e in range(1,end,1):
        print(e)

        # read in control simulations
        T = read_in(s,m,e,'T')
        Q = read_in(s,m,e,'Q')
        U = read_in(s,m,e,'U')
        V = read_in(s,m,e,'V')
        time = T.time
        lat = T.lat
        lon = T.lon

        # even though topo and the 3D variables are on the same horizontal grid, it seems that there are some extremely small
        # differences induced by truncation. We will just overwrite topo's lat with the 3D variables' lat.
        topo['lat'] = T.lat.data

        # calculate total energy
        TE = (T.T * 1004) + (Q.Q*2.5e6) + (U.U**2 + V.V**2)**0.5 + topo*9.81

        # save some space
        del Q, U, V

        # regrid to plevs
        TE_regrid = ngl.vinth2p(TE,T.hyam,T.hybm,plevs,T.PS,1,T.P0/100,1,kxtrp=False)
        TE_regrid[(TE_regrid==1e30)] = np.nan

        dp = np.diff(np.concatenate([np.zeros((730,1,192,288)),np.broadcast_to(np.expand_dims(mids,axis=[0,2,3]),(730,44,192,288)),
                np.expand_dims(np.array(T.PS),axis=1)/100],axis=1),axis=1)
        dp[(dp<0)] = np.nan

        del T, TE
        # vertically integrate
        # 100 is to convert hPa to Pa
        TE_vi = np.nansum(TE_regrid * dp * 100,axis=1)/9.81
        
        del TE_regrid, dp

        # turn into dataarray
        da = xr.DataArray(data=TE_vi,dims=['time','lat','lon'],coords=dict(
        lon=(["lon"], lon.data),
        lat=(["lat"], lat.data),
        time=time,
            ),
                attrs=dict(
                units="W/m^2",
            ),
        ).rename('TE_vi')

        # take time derivative
        dEdt = da.differentiate(coord='time',edge_order=1,datetime_unit='s').rename('dEdt')
        del TE_vi

        # remove TOA and surface fluxes
        SHFLX = read_in(s,m,e,'SHFLX').SHFLX
        LHFLX = read_in(s,m,e,'LHFLX').LHFLX
        FSNS = -1 *read_in(s,m,e,'FSNS').FSNS # this is so that up = positive
        FLNS = read_in(s,m,e,'FLNS').FLNS

        FLNT = -1 * read_in(s,m,e,'FLNT').FLNT
        FSNT = read_in(s,m,e,'FSNT').FSNT

        AHT = dEdt - (SHFLX + LHFLX + FSNS + FLNS) - (FLNT + FSNT)
        AHT = AHT.rename('AHT')
        
        del SHFLX, LHFLX, FSNS, FLNS, FLNT, FSNT

        # save out
        out = xr.merge([da,AHT,dEdt])
        pathout = '/dx02/janoski/cesm/ctrl_4xCO2_transports/b40.1850.cam5-lens.'+s+'.'+str(f"{m:02d}")+'.'+str(f"{e:02d}")+'.h1_Fconv_resi.nc'
        out.to_netcdf(pathout)

        del AHT, dEdt, da, out
