'''
PURPOSE: To use the ERA-I grib files of temperature, surface pressure and specific humidity          to create: 
         and generate two layer theta_e variables that can combine to produce a buoyancy-like
         variable

AUTHOR: Fiaz Ahmed

DATE: 03/7/17
'''

import numpy as np
from netCDF4 import Dataset
from glob import glob
import datetime as dt
from dateutil.relativedelta import relativedelta
import time
import itertools
from mpi4py import MPI
from sys import exit
from numpy import dtype
from parameters import *
from vert_cython import vert_integ_variable_bl,vert_integ_exneri_variable_bl,vert_integ_lt_variable_bl
from vert_cython import find_closest_index

####### MASK ########
f=Dataset('/glade/u/home/fiaz/era_process/binning/land_sea_mask_era.nc','r')
lsm=np.asarray(f.variables['LSMASK'][:],dtype='float')
latm,lonm=f.variables['lat'][:],f.variables['lon'][:]
f.close()

mask_ocean=np.copy(lsm)
mask_ocean[mask_ocean!=1]=np.nan

mask_land=np.copy(lsm)
mask_land[mask_land!=0]=np.nan
mask_land[mask_land==0]=1.

### Create list of file to process

strt_date=dt.datetime(2013,6,4)
# end_date=dt.datetime(2001,9,2)
end_date=dt.datetime(2014,12,31)

dts=[]

dirc='/glade/p/univ/p35681102/fiaz/erai_data/regridded/'

dirc1=dirc+'era_T/'
dirc2=dirc+'era_q/'
dirc3=dirc+'era_surfTq/'
dirc4=dirc+'era_surfp/'
dirc5=dirc+'era_geop/'
    
list1=[]
list2=[]
list3=[]
list4=[]
    
while strt_date<=end_date:

    d1=strt_date.strftime("%Y-%m-%d")
    yr=strt_date.strftime("%Y")
    dts.append(strt_date)
    fname=dirc1+'era_vertT_regridded_'+d1+'*'
    list_temp=(glob(fname))
    list_temp.sort()
    list1.append(list_temp)
    strt_date+=relativedelta(days=1)
    
chain1=itertools.chain.from_iterable(list1)
list1= (list(chain1))

f=Dataset(list1[0],'r')
lev=f.variables['level'][:]
lat=f.variables['lat'][:]
lon=f.variables['lon'][:]
f.close()

# print(lev,type(lev))
lev=np.asarray(lev) 

i200=np.argmin(abs(lev-200))
i500=np.argmin(abs(lev-500))
i1000=np.argmin(abs(lev-1000))

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
        jobs=list1
        jobs=split(jobs,comm.size)

else:
        jobs=None

# Scatter jobs across cores
jobs=comm.scatter(jobs,root=0)
# 
s=time.time()

# Use this if MPI is off:
# jobs=list1

for j in jobs:
    d1=j[-13:-3]
    print(j)

    fname2=dirc2+'era_vertq_regridded_'+d1+'*'
    fname3=dirc3+'era_surfTq_regridded_'+d1+'*'
    fname4=dirc4+'era_surfp_regridded_'+d1+'*'
    fname5=dirc5+'era_vertZ_regridded_'+d1+'*'

    list2=(glob(fname2))
    list3=(glob(fname3))
    list4=(glob(fname4))
    list5=(glob(fname5))

    ## LOAD TEMPERATURE ##

    print('LOADING TEMP.')
    f=Dataset(j,'r')
    t=f.variables['T'][:]
    f.close()
    
    ## LOAD SPECIFIC HUMIDITY (UNITS OF K)##

    print('LOADING SP.HUM.')
    f=Dataset(list2[0],'r')
    q=f.variables['Q'][:]
    f.close()

    ## LOAD 2 meter TEMP. & HUMIDITY ##

    print('LOADING SURF. T & Q')
    f=Dataset(list3[0],'r')
    dt2m=f.variables['Q'][:]
    t2m=f.variables['T'][:]
    f.close()
    
    ## LOAD surface geopotential ##

    print('LOADING SURF. GEOPOTENTIAL')
    f=Dataset(dirc+'/era_surf_geo_regridded.nc','r')
    Zsurf=f.variables['geop_surf'][:]
    f.close()
    Zsurf=Zsurf[None,:,:]
    
    ## LOAD SURFACE PRESSURE ######

    print('LOADING SURF. PRESSURE')
    f=Dataset(list4[0],'r')
    pres=f.variables['pres'][:]*1e-2 # Converting to hPa            
    f.close()    
        
    pbl_top=pres-100 ## The sub-cloud layer is 150 mb thick ##
    low_top=np.zeros_like(pres)
    low_top[:]=500
#     low_top=pbl_top-400 ## The next layer is 350 mb thick ##

    pbl_top=np.float_(pbl_top.flatten())
    low_top=np.float_(low_top.flatten())
    lev=np.float_(lev)

    print('PREPPING LAYER AVERAGING-2')
    
    pbl_ind=np.zeros(pbl_top.size,dtype=np.int64)
    low_ind=np.zeros(low_top.size,dtype=np.int64)

    find_closest_index(pbl_top,lev,pbl_ind)
    find_closest_index(low_top,lev,low_ind)

    pres_3d=np.zeros_like(t)
    pres_3d[:]=pres[:,None,:,:]

    levels=np.zeros_like(t)        
    levels[:]=lev[None,:,None,None]
    
     ##### Calculate qsat #####

    Tk0 = 273.15 # Reference temperature.
    Es0 = 610.7 # Vapor pressure at Tk0.
    Lv0 = 2500800 # Latent heat of evaporation at Tk0.
    cpv = 1869.4 # Isobaric specific heat capacity of water vapor at tk0.
    cl = 4218.0 # Specific heat capacity of liquid water at tk0.
    R = 8.3144 # Universal gas constant.
    Mw = 0.018015 # Molecular weight of water.
    Rv = R/Mw # Gas constant for water vapor.
    Ma = 0.028964 # Molecular weight of dry air.
    Rd = R/Ma # Gas constant for dry air.
    epsilon = Mw/Ma 
    g = 9.80665 
    cpd=1004.

    ### TROPOSPHERIC MOISTURE AND TEMPERATURE CALCULATIONS ###


    ### Saturation specific humidity and mixing ratio ###
    print('ESTIMATING qSAT.')

#     Es=(Es0*(t/Tk0)**((cpv-cl)/Rv))*np.exp((Lv0+(cl-cpv)*Tk0)/Rv*(1/Tk0-1./t))
    Es=es_calc(t)
    ws=(epsilon)*(Es/levels)    
#     qsat=(Lv/Cp)*ws/(1+ws)
    qsat=ws/(1+ws)
    
    w=q/(1-q)
    e=w*levels/(epsilon+w) # vapor pressure
    
    ### END ###

    ### SURFACE MOISTURE AND TEMPERATURE CALCULATIONS ###

    ### Get 2 metre specific humidity ###

    ### 2m specific humidity and mixing ratio ###
    print('ESTIMATING 2m SP.HUM.')

    eps=(Es0*(dt2m/Tk0)**((cpv-cl)/Rv))*np.exp((Lv0+(cl-cpv)*Tk0)/Rv*(1/Tk0-1./dt2m))
    eps=eps*1e-2
    wps=(epsilon)*(eps/pres)
#     qps=(Lv/Cp)*wps/(1+wps) 
    qps=wps/(1+wps) 

    ### 2m qsat and saturation specific humidity ###

    print('ESTIMATING 2m qSAT')

    Esat_ps=(Es0*(t2m/Tk0)**((cpv-cl)/Rv))*np.exp((Lv0+(cl-cpv)*Tk0)/Rv*(1/Tk0-1./t2m))
    Esat_ps=Esat_ps*1e-2
    wsat_ps=(epsilon)*(Esat_ps/pres)
#     qsat_ps=(Lv/Cp)*wsat_ps/(1+wsat_ps) 
    qsat_ps=wsat_ps/(1+wsat_ps) 
    
    #### END ########
    
    print('FILLING SUB-SURFACE PRESSURE LEVELS')

    ##### Fill all pressure level below surface with surface values ###

    var_4d=np.zeros_like(t)

    ## Surface temperature
    var_4d[:]=t2m[:,None,:,:]
    t[levels>=pres_3d]=var_4d[levels>=pres_3d]

    ## Surface specific humidity 
    var_4d[:]=qps[:,None,:,:]
    q[levels>=pres_3d]=var_4d[levels>=pres_3d]
    
    ### Surface pressure
    dps=np.zeros_like(pres) ## Thickness of surface layer if surface_pres>1000 mb
    dps[pres>levels[:,-1,:,:]]=pres[pres>levels[:,-1,:,:]]-levels[:,-1,:,:][pres>levels[:,-1,:,:]]
    levels[levels>=pres_3d]=pres_3d[levels>=pres_3d]
    dp=np.diff(levels,axis=1)
    
    print('COMPUTING THETA_E')

    ### Surface theta_e values ###
    
    psd=pres-eps # partial pressure of dry air
    rh_ps=qps/qsat_ps
    theta_e_ps=(t2m*(Po/psd)**(Rd/cpd))*((rh_ps)**(-wps*Rv/cpd))*np.exp(Lv0*wps/(cpd*t2m))
    theta_e_sat_ps=(t2m*(Po/psd)**(Rd/cpd))*np.exp(Lv0*wsat_ps/(cpd*t2m))
    theta_e_sub_sat_ps=theta_e_sat_ps-theta_e_ps

    ## RH ###
    rh=q/qsat
    pd=levels-e # partial pressure of dry air

   #Calculate theta_e
    theta_e=(t*(Po/pd)**(Rd/cpd))*((rh)**(-w*Rv/cpd))*np.exp(Lv0*w/(cpd*t))
   #Saturated theta_e
    theta_e_sat=(t*(Po/pd)**(Rd/cpd))*np.exp(Lv0*ws/(cpd*t))

    ###### INTEGRATED QUANTITIES ###### 
    mse=theta_e
    mse_sat=theta_e_sat
    mse_surf=theta_e_ps
    mse_sat_surf=theta_e_sat_ps
    
    print('PREPPING LAYER AVERAGING-3')

    pbl_top=np.asarray([lev[i] for i in pbl_ind])
    low_top=np.asarray([lev[i] for i in low_ind])

    pbl_top=pbl_top.flatten()
    pbl_ind=pbl_ind.flatten()
    low_ind=low_ind.flatten()


    mse=np.swapaxes(mse,0,1)
    mse=mse.reshape(lev.size,-1)

    ### Get the lower layer temp. and qsat ###
    t=np.swapaxes(t,0,1)
    t=t.reshape(lev.size,-1)

    qsat=np.swapaxes(qsat,0,1)
    qsat=qsat.reshape(lev.size,-1)
    ##########################################

    mse_sat=np.swapaxes(mse_sat,0,1)
    mse_sat=mse_sat.reshape(lev.size,-1)
    mse_surf=mse_surf.flatten()

    pres=pres.flatten()

    dps=np.float64(dps.flatten())

    dp=np.swapaxes(dp,0,1)
    dp=dp.reshape(lev.size-1,-1)


    ### MSE VERT. ###

    mse_vert=(mse[1:,:]+mse[:-1,:])*0.5
    mse_lt_vert=(mse[1:,:]+mse[:-1,:])*0.5
    mse_sat_vert=(mse_sat[1:,:]+mse_sat[:-1,:])*0.5
    mse_last=(mse[-1,:]+mse_surf)*0.5*dps

    qsat_vert=(qsat[1:,:]+qsat[:-1,:])*0.5
    t_vert=(t[1:,:]+t[:-1,:])*0.5

    mse_bl=np.zeros((pres.size))
    mse_lt=np.zeros((pres.size))

    qsat_lt=np.zeros((pres.size))
    t_lt=np.zeros((pres.size))

    mse_sat_lt=np.zeros((pres.size))
    mse_mt=np.zeros(pres.size)
    mse_sat_mt=np.zeros(pres.size)

    t1=time.time()
    mse_vert=(np.float_(mse_vert))
    qsat_vert=(np.float_(qsat_vert))
    t_vert=(np.float_(t_vert))

    mse_surf=(np.float_(mse_surf))
    mse_last=(np.float_(mse_last))
    mse_sat=(np.float_(mse_sat))
    pbl_ind=np.int_(pbl_ind)
    lev=np.float_(lev)
    dp=np.float_(dp)
    low_ind=np.int_(low_ind)

    print('VERTICAL INTEGRATION')

#     vert_integ_variable_bl(mse_vert,mse_last,mse_sat,
#     pbl_ind,lev,dp,mse_bl,mse_lt,mse_sat_lt,low_ind)

#     vert_integ_variable_bl(mse_vert,mse_last,mse_sat,
#     pbl_ind,lev,dp,mse_bl,mse_lt,mse_sat_lt,mse_mt,mse_sat_mt,
#     i500,low_ind)

    vert_integ_variable_bl(mse_vert,mse_last,mse_sat,
    pbl_ind,lev,dp,mse_bl,mse_lt,mse_sat_lt,low_ind)
    
    mse_bl=mse_bl.reshape(4,lat.size,lon.size)
    mse_lt=mse_lt.reshape(4,lat.size,lon.size)
    mse_sat_lt=mse_sat_lt.reshape(4,lat.size,lon.size)
    mse_mt=mse_mt.reshape(4,lat.size,lon.size)
    mse_sat_mt=mse_sat_mt.reshape(4,lat.size,lon.size)

    pbl_top=pbl_top.reshape(4,lat.size,lon.size)
    low_top=low_top.reshape(4,lat.size,lon.size)
    pres=pres.reshape(4,lat.size,lon.size)

    mse_bl/=(pres-pbl_top)
    mse_lt/=(pbl_top-low_top)
    mse_sat_lt/=(pbl_top-low_top)

    print('Maximum:')
    print('thetae BL:',np.nanmax(mse_bl),np.nanmin(mse_bl))
    print('thetae LT:',np.nanmax(mse_lt),np.nanmin(mse_lt))
    print('thetae SAT LT:',np.nanmax(mse_sat_lt),np.nanmin(mse_sat_lt))

#     print('EXNERI BL:',np.nanmax(exneri_bl),np.nanmin(exneri_bl))
#     print('EXNERI LT:',np.nanmax(exneri_lt),np.nanmin(exneri_lt))
    
    print('---------------------')

    print('Mean:')

    print('thetae BL:',np.nanmean(mse_bl))
    print('thetae LT:',np.nanmean(mse_lt))
    print('thetae SAT LT:',np.nanmean(mse_sat_lt))

#     exit()


#     print('MSE BL:',mse_bl.max(),mse_bl.min())
#     print('MSE LT:',mse_lt.max(),mse_lt.min())
#     print('MSE SAT LT:',mse_sat_lt.max(),mse_sat_lt.min())

    print('SAVING FILE')

#     fout='/glade/p/univ/p35681102/fiaz/erai_data/layer_moist_static_energy/'+'era_layer_Lqcptgz_'+d1+'.nc'
    fout='/glade/p/univ/p35681102/fiaz/erai_data/layer_thetae/'+'era_2layers_thetae_'+d1+'.nc'

    # fout='/glade/scratch/fiaz/era_processed/'+'era_trop_0.25_'+d1+'.nc'

    ##### SAVE FILE ######

    try:ncfile.close()
    except:pass

    ncfile = Dataset(fout, mode='w', format='NETCDF4')

    ncfile.createDimension('time',q.shape[0])
    ncfile.createDimension('lat',lat.size)
    ncfile.createDimension('lon',lon.size)

    dy = ncfile.createVariable('time','i4',('time'))
    lt = ncfile.createVariable('lat',dtype('float32').char,('lat'))
    ln = ncfile.createVariable('lon',dtype('float32').char,('lon'))

    dy.units = "days since 1800-01-01 00:00" ;
    dy.long_name = "initial time" ;
    dy.calendar = "gregorian"

    mbl = ncfile.createVariable('thetae_bl',dtype('float32').char,('time','lat','lon'),zlib=True)
    mlt = ncfile.createVariable('thetae_lt',dtype('float32').char,('time','lat','lon'),zlib=True)
    mslt = ncfile.createVariable('thetae_sat_lt',dtype('float32').char,('time','lat','lon'),zlib=True)
    mmt = ncfile.createVariable('thetae_mt',dtype('float32').char,('time','lat','lon'),zlib=True)
    msmt = ncfile.createVariable('thetae_sat_mt',dtype('float32').char,('time','lat','lon'),zlib=True)

#     ebl = ncfile.createVariable('exneri_bl',dtype('float32').char,('time','lat','lon'),zlib=True)
#     elt = ncfile.createVariable('exneri_lt',dtype('float32').char,('time','lat','lon'),zlib=True)

    mbl.description="thetae averaged from surface to 100 mb up" 
    mlt.description="thetae averaged from 100 mb above surface to 500 mb "
    mslt.description="SATURATED thetae averaged from 100 mb above surface to 500 mb "

    dy[:]=np.asarray(t.shape[0])
    lt[:]=lat
    ln[:]=lon

    mbl[:]=mse_bl
    mlt[:]=mse_lt
    mslt[:]=mse_sat_lt
    mmt[:]=mse_mt
    msmt[:]=mse_sat_mt

    ncfile.close()

    print('FILE WRITTEN')
    et=time.time()
    #print et-s
