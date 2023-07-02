# Purpose: Read in feedbacks calculated in spatial_avg_feedbacks.py and
# normalize by the globally averaged Planck feedback.

# By: Ty Janoski
# Updated: 05.19.23

# import statements
import xarray as xr
import numpy as np

# request user input to determine which month to use
m = int(input('Enter 1 or 7 for January or July initializations, respectively: '))
if(m!=1 and m!=7):
    print("That's not 1 or 7.")
    sys.exit()
    
r = str(input('Enter arc, glb, or ant60 for the Arctic, global, or Antarcic average, respectively: '))
if(r!='arc' and r != 'glb' and r != 'ant60'):
    print('Invalid input.')
    sys.exit()
    
# read in global average Planck feedback and surface air temperature
planck = xr.open_dataarray('/dx02/janoski/cesm/spat_avg_feedbacks/all/b40.1850.cam5-lens.'+
                         str(f"{m:02d}"+'.glb_planck_tropo_TOA.nc'))
dSAT = xr.open_dataarray('/dx02/janoski/cesm/spat_avg_feedbacks/all/b40.1850.cam5-lens.'+
                         str(f"{m:02d}"+'.glb_dSAT.nc'))
Ts = xr.open_dataarray('/dx02/janoski/cesm/spat_avg_feedbacks/all/b40.1850.cam5-lens.'+
                         str(f"{m:02d}"+'.glb_Ts_tropo_sfc.nc'))

# normalization factor is the planck feedback divided by the
# change in SAT (in our framework, at least)
norm = -1 * planck.mean(dim='ens')/dSAT.mean(dim='ens')
glb_Ts = Ts.mean(dim='ens')/dSAT.mean(dim='ens')

# create list to iterate through all variable files
domains = ['all','land_only','ocean_only']
for d in domains:
    lead = '/dx02/janoski/cesm/spat_avg_feedbacks/'+d+'/b40.1850.cam5-lens.'+str(
        f"{m:02d}"+'.'+r+'_')

    file_endings = ['AHT.nc','dEdt.nc','alb_TOA.nc','alb_sfc.nc','cloud_LW_TOA.nc',
                    'cloud_LW_sfc.nc','cloud_SW_TOA.nc','cloud_SW_sfc.nc',
                    'lapse_tropo_TOA.nc','planck_tropo_TOA.nc',
                    'q_LW_tropo_TOA.nc','q_SW_tropo_TOA.nc','q_LW_tropo_sfc.nc',
                    'q_SW_tropo_sfc.nc','Sconv.nc','Fsfc.nc','dLHFLX.nc',
                    'dSHFLX.nc','Ta_tropo_sfc.nc','Ts_tropo_sfc.nc',
                    'T_strato_TOA.nc','q_LW_strato_TOA.nc','q_SW_strato_TOA.nc',
                    'Wconv.nc','CO2_TOA.nc','CO2_sfc.nc']

    # iterate through files

    for f in file_endings:
        print(f[:-3])
        # open file, rename variable
        fin = xr.open_dataarray(lead+f)
        fin = fin.rename(f[:-3]).copy()

        # planck is defined as the local deviation
        if(f=='planck_tropo_TOA.nc'):
            dSAT = xr.open_dataarray(lead+'dSAT.nc')
            local = fin/dSAT.mean(dim='ens')
            normalized = (local + norm)*dSAT.mean(dim='ens')/norm
        # elif(f=='Ts_tropo_sfc.nc'):
        #     dSAT = xr.open_dataarray(lead+'dSAT.nc')
        #     local = fin/dSAT.mean(dim='ens')
        #     normalized = (local - glb_Ts) * dSAT.mean(dim='ens')/norm
        else:
            # normalize the variable
            normalized = fin/norm

        # save it out
        pathout = lead.replace('feedbacks/','feedbacks/norm_feedbacks/') + f
        print(pathout)
        normalized.to_netcdf(pathout)
        
# create list to iterate through all variable files
#domains = ['all']
#for d in domains:
#    lead = '/dx02/janoski/cesm/spat_avg_feedbacks/'+d+'/b40.1850.cam5-lens.'+str(
#        f"{m:02d}"+'.'+r+'_')

#    file_endings = ['CO2_TOA.nc','CO2_sfc.nc']

    # iterate through files

#    for f in file_endings:
#        print(f[:-3])
#        # open file, rename variable
#        fin = xr.open_dataarray(lead+f)
#        fin = fin.rename(f[:-3]).copy()

        # planck is defined as the local deviation
#        if(f=='planck_tropo_TOA.nc'):
#            dSAT = xr.open_dataarray(lead+'dSAT.nc')
#            local = fin/dSAT.mean(dim='ens')
#            normalized = (local + norm)*dSAT.mean(dim='ens')/norm
#        else:
            # normalize the variable
#            normalized = fin/norm

        # save it out
#        pathout = lead.replace('feedbacks/','feedbacks/norm_feedbacks/') + f
#        print(pathout)
#        normalized.to_netcdf(pathout)
