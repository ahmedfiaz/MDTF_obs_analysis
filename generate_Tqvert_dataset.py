'''
NAME: generate_bint_dataset.py

PURPOSE: To generate save co-incident values of TRMM precip. cwv and T_hat
         Only use trop. ocean data. Create single-D arrays for each month.
AUTHOR: Fiaz Ahmed

DATE: 11/24/19

'''
### !!! Run this to produce annual average if the files are already written !!!!

import numpy as np
import glob
from netCDF4 import Dataset
from numpy import dtype
import datetime as dt
import itertools
from dateutil.relativedelta import relativedelta
from mpi4py import MPI
from sys import exit
from glob import glob
from bin_parameters import * 

diri_pcp='/glade/p/univ/p35681102/fiaz/T3B42/'
diri_Tq='/glade/p/univ/p35681102/fiaz/erai_data/regridded/era_Tqml/'
diri_surfp='/glade/p/univ/p35681102/fiaz/erai_data/regridded/era_surfp/'
# diri_cwv='/glade/p/univ/p35681102/fiaz/erai_data/layer_moist_static_energy/'

months_jja=[6,7,8]
months_djf=[12,1,2]

list_trmm=[]

strt_date=dt.datetime(2002,1,1)
end_date=dt.datetime(2002,1,31)

mask_land=np.copy(lsm)
mask_ocean=np.copy(lsm)

mask_land[mask_land!=0]=np.nan
mask_land[mask_land==0]=1.
mask_ocean[mask_ocean!=1]=np.nan


while strt_date<=end_date:
    d1=strt_date.strftime("%Y%m")
    fname=diri_pcp+'TRMM.3B42.'+str(d1)+'*.nc'
    list_temp=(glob(fname))
    list_temp.sort()
    list_trmm.append(list_temp)
    strt_date+=relativedelta(months=1)

    
chain1=itertools.chain.from_iterable(list_trmm)
list_trmm= (list(chain1))


f=Dataset(list_trmm[0],'r')
lat=f.variables['latitude'][:]
lon=f.variables['longitude'][:]
f.close()


### Scatter Jobs #####
comm=MPI.COMM_WORLD
print(comm.rank)

def split(container, count):
    """
    Simple function splitting a container into equal length chunks.
    Order is not preserved but this is potentially an advantage depending on
    the use case.
    """
    return [container[_i::count] for _i in range(count)]

if comm.rank == 0:
        jobs=list_trmm
        jobs=split(jobs,comm.size)

else:
        jobs=None

### Scatter jobs across cores
jobs=comm.scatter(jobs,root=0)


# cape_dry=np.zeros((lat.size,lon.size))
# subsat_dry=np.zeros((lat.size,lon.size))
# bint_dry=np.zeros((lat.size,lon.size))
# 
# cape_clim=np.zeros((lat.size,lon.size))
# subsat_clim=np.zeros((lat.size,lon.size))
# bint_clim=np.zeros((lat.size,lon.size))

# for i in list_trmm:
for i in jobs:

    print('Opening TRMM')

    d1=i[-12:-6]
    dts=dt.datetime.strptime(str(d1), "%Y%m")

    f=Dataset(i,'r')
    prc=f.variables['precip_trmm'][:]
    f.close()
    
    prc[prc<0]=np.nan
    
    d2=dts.strftime("%Y-%m")
    fname2=diri_Tq+'era_vertTq_regridded_'+d2+'*'
    list2=(glob(fname2))
    list2.sort()
    
    fname3=diri_surfp+'era_surfp_regridded_'+d2+'*'
    list3=glob(fname3)
    list3.sort()
    
    f2=Dataset(list2[0],'r')
    lev=f2.variables['level'][:]
    f2.close()


    print('Opening ERA-I')

    for l,(k,m) in enumerate(zip(list2,list3)):
#     for l,k in enumerate(list2):

        i1,i2=l*4,(l+1)*4
        f=Dataset(k,'r')
        lv_HYBL1_a, lv_HYBL1_b=f.variables['lv_HYBL1_a'][:],f.variables['lv_HYBL1_b'][:]
        Ttemp,qtemp=f.variables['T'][:],f.variables['q'][:]
        Ttemp=np.swapaxes(Ttemp,0,1)#*mask_land 
        qtemp=np.swapaxes(qtemp,0,1)#*mask_land 
        f.close()

        f=Dataset(m,'r')
        surfp=f.variables['pres'][:]
        f.close()

        ### Slice oceans and preserve vertical structure information ###
        
        prc_temp=prc[i1:i2,...]*mask_land
        
        mask_ind=np.where(np.isfinite(prc_temp))
        prc_slice=prc[mask_ind[0],mask_ind[1],mask_ind[2]]
        T_slice=Ttemp[:,mask_ind[0],mask_ind[1],mask_ind[2]]
        q_slice=qtemp[:,mask_ind[0],mask_ind[1],mask_ind[2]]
        surfp_slice=surfp[mask_ind]
        
        print 'SAVING FILE'
        fout='/glade/p/univ/p35681102/fiaz/erai_data/cwv_that_erai_trmm_ocn/daily_files/'
        filo=fout+'Tqvert_ocns_'+str(k[-13:-3])+'.nc'

        try:ncfile.close()
        except:pass
    
        ncfile = Dataset(filo, mode='w', format='NETCDF4')
    
        ncfile.createDimension('lev',lv_HYBL1_a.size)
        ncfile.createDimension('hor',prc_slice.size)
    
#         lv = ncfile.createVariable('lev',dtype('float32').char,('lev'))
#         hor = ncfile.createVariable('hor',dtype('float32').char,('hor'))
    
#         pred_pcp_std = ncfile.createVariable('predicted_precip_90p',dtype('float32').char,('lat','lon'),zlib=True)
        lv_HYBL1_a_var=ncfile.createVariable('lv_HYBL1_a',dtype('float32').char,('lev'),zlib=True)
        lv_HYBL1_b_var=ncfile.createVariable('lv_HYBL1_b',dtype('float32').char,('lev'),zlib=True)
        prc_var = ncfile.createVariable('precip_ocn',dtype('float32').char,('hor'),zlib=True)
        surfp_var=ncfile.createVariable('surfp_ocn',dtype('float32').char,('hor'),zlib=True)
        temp_var= ncfile.createVariable('temp_ocn',dtype('float32').char,('lev','hor'),zlib=True)
        sphum_var= ncfile.createVariable('sphum_ocn',dtype('float32').char,('lev','hor'),zlib=True)
    
#         tm[:]=np.arange(4)
        lv_HYBL1_a_var[:]=lv_HYBL1_a   
        lv_HYBL1_b_var[:]=lv_HYBL1_b   
#         hor[:]=np.arange(prc_var.size)
        prc_var[:]=prc_slice
        surfp_var[:]=surfp_slice
        temp_var[:]=T_slice
        sphum_var[:]=q_slice
        
        ncfile.close()
        
        print('DONE SAVING')

    


        # 
# 
# 
#         exit()
# 
#         if l==0:
#             Tocn=np.copy(T_slice)
#             qocn=np.copy(q_slice)
#             prc_ocn=np.copy(prc_slice)
#             surfp_ocn=np.copy(surfp_slice)
#             
#         elif l>0:
# #             print(Ttemp.shape,T_slice.shape)
#             Tocn=np.hstack((Tocn,T_slice))
#             qocn=np.hstack((qocn,q_slice))
#             prc_ocn=np.hstack((prc_ocn,prc_slice))
#             surfp_ocn=np.hstack((surfp_ocn,surfp_slice))
#                     
#     print(Tocn.shape,qocn.shape,prc_ocn.shape,surfp_ocn.shape)
#     
#     dts=dt.datetime.strptime(str(d1), "%Y%m")



#     fout='/glade/p/univ/p35681102/fiaz/predicted_protected_precip/'+'global_protected_dilute_counts_'+strt_end_yr+'.nc'
#     fout='/glade/p/univ/p35681102/fiaz/erai_data/cwv_that_erai_trmm_ocn/'
    ##### SAVE FILE ######
# 
#     try:ncfile.close()
#     except:pass
# 
#     ncfile = Dataset(fout, mode='w', format='NETCDF4')
# 
#     ncfile.createDimension('lat',lat.size)
#     ncfile.createDimension('lon',lon.size)
#     ncfile.createDimension('time',None)
# 
#     lt = ncfile.createVariable('lat',dtype('float32').char,('lat'))
#     ln = ncfile.createVariable('lon',dtype('float32').char,('lon'))
#     # tm = ncfile.createVariable('time',dtype('float32').char,('time'))
# 
#     # pred_pcp_std = ncfile.createVariable('predicted_precip_90p',dtype('float32').char,('lat','lon'),zlib=True)
#     cnts_dil = ncfile.createVariable('counts_dil_bc_exceedance',dtype('float32').char,('lat','lon'),zlib=True)
#     cnts_dil_jja = ncfile.createVariable('counts_dil_bc_exceedance_jja',dtype('float32').char,('lat','lon'),zlib=True)
#     cnts_dil_djf = ncfile.createVariable('counts_dil_bc_exceedance_djf',dtype('float32').char,('lat','lon'),zlib=True)
# 
#     cnts_prot = ncfile.createVariable('counts_prot_bc_exceedance',dtype('float32').char,('lat','lon'),zlib=True)
#     cnts_prot_jja = ncfile.createVariable('counts_prot_bc_exceedance_jja',dtype('float32').char,('lat','lon'),zlib=True)
#     cnts_prot_djf = ncfile.createVariable('counts_prot_bc_exceedance_djf',dtype('float32').char,('lat','lon'),zlib=True)
# 
# 
#     # tm[:]=np.arange(4)
#     lt[:]=lat
#     ln[:]=lon
# 
#     cnts_dil[:]=counts_dil
#     cnts_dil_jja[:]=counts_dil_jja
#     cnts_dil_djf[:]=counts_dil_djf
# 
#     cnts_prot[:]=counts_prot
#     cnts_prot_jja[:]=counts_prot_jja
#     cnts_prot_djf[:]=counts_prot_djf
# 
#     ncfile.close()
    
   #  print('SAVING FILE')
# 
#     fout_T='/glade/p/univ/p35681102/fiaz/erai_data/cwv_that_erai_trmm_ocn/'+'Tvert_ocns.'+d1
#     fout_q='/glade/p/univ/p35681102/fiaz/erai_data/cwv_that_erai_trmm_ocn/'+'qvert_ocns.'+d1
#     fout_Psurfp='/glade/p/univ/p35681102/fiaz/erai_data/cwv_that_erai_trmm_ocn/'+'prec_surfp_ocns.'+d1
# 
# #     np.savez_compressed(fout,Tvert_ocean=Tocn,qvert_ocean=qocn,
# #     prc_trmm_ocean=prc_ocn,surfp_ocean=surfp_ocn)
# 
# #     print(fout)
# #     print(type(Tocn.filled()),type(qocn.filled()))
# #     print(Tocn.shape,Tocn.filled().shape)
# #     print(np.max(qocn.filled()),np.min(qocn.filled()))
#     np.savez_compressed(fout_T,
#     Tvert_ocean=Tocn.filled())
#     
#     np.savez_compressed(fout_q,
#     qvert_ocean=qocn.filled())
#     
#     np.savez_compressed(fout_Psurfp,
#     prc_trmm_ocean=prc_ocn.filled(),
#     surfp_ocean=surfp_ocn.filled())





    
