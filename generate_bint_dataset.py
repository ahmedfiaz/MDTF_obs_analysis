'''
NAME: generate_bint_dataset.py

PURPOSE: To generate B_int data from ERA-I, and also include corresponding TRMM precip.
         To save:
            a) B_int
            b) CAPE component
            c) Subsat component
            d) weights w_b and w_L
            e) precip.

AUTHOR: Fiaz Ahmed

DATE: 08/23/19

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

ref_thetae=340 ## reference theta_e in K to convert buoy. to temp units
gravity=9.8 ##
thresh_pres=700 ## Filter all point below this surface pressure in hPa

diri_pcp='/glade/p/univ/p35681102/fiaz/T3B42/'
diri_buoy='/glade/p/univ/p35681102/fiaz/erai_data/layer_thetae/'
diri_surfp='/glade/p/univ/p35681102/fiaz/erai_data/regridded/era_surfp/'
# diri_buoy='/glade/p/univ/p35681102/fiaz/erai_data/layer_moist_static_energy/'

months_jja=[6,7,8]
months_djf=[12,1,2]

list_trmm=[]

strt_date=dt.datetime(2013,7,1)
end_date=dt.datetime(2014,12,31)


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
print comm.rank

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
    
    d2=dts.strftime("%Y-%m")
    fname2=diri_buoy+'era_2layers_thetae_'+d2+'*'
    list2=(glob(fname2))
    list2.sort()
    
    fname3=diri_surfp+'era_surfp_regridded_'+d2+'*'
    list3=(glob(fname3))
    list3.sort()

    thetae_bl=np.zeros_like(prc)
    thetae_lt=np.zeros_like(prc)
    thetae_sat_lt=np.zeros_like(prc)
    surfp=np.zeros_like(prc)

        
    print('Opening ERA-I')

    for l,(k,m) in enumerate(zip(list2,list3)):
    
        i1,i2=l*4,(l+1)*4
        f=Dataset(k,'r')
        thetae_bl[i1:i2,...]=f.variables['thetae_bl'][:]
        thetae_lt[i1:i2,...]=f.variables['thetae_lt'][:]
        thetae_sat_lt[i1:i2,...]=f.variables['thetae_sat_lt'][:]
        f.close()
        
        f=Dataset(m,'r')
        surfp[i1:i2,...]=f.variables['pres'][:]/100
        f.close()

    print(np.nanmax(thetae_bl),np.nanmax(thetae_lt),np.nanmax(thetae_sat_lt))
    print(np.nanmin(thetae_bl),np.nanmin(thetae_lt),np.nanmin(thetae_sat_lt))

    delta_pl=surfp-100-500
    delta_pb=100
    
    wb=(delta_pb/delta_pl)*np.log((delta_pl+delta_pb)/delta_pb)
    wl=1-wb
    
    ### Filter all points whose surface pressure is < 750 mb ###
    ### these include high-altitude points. Orographic rain is dominated by convergence
    ### and perhaps less likely by the thermodynamics (?). A future framework to capture
    ### the dynamics will be useful.
    
    wb[surfp<thresh_pres]=np.nan
    wl[surfp<thresh_pres]=np.nan
    
    cape=ref_thetae*(thetae_bl-thetae_sat_lt)/thetae_sat_lt
    subsat=ref_thetae*(thetae_sat_lt-thetae_lt)/thetae_sat_lt    
    bint=gravity*(wb*(thetae_bl-thetae_sat_lt)/thetae_sat_lt-wl*(thetae_sat_lt-thetae_lt)/thetae_sat_lt)
    
    cape[surfp<thresh_pres]=np.nan
    subsat[surfp<thresh_pres]=np.nan
    bint[surfp<thresh_pres]=np.nan
    
    
    print(np.nanmax(bint),np.nanmin(bint))
    print(np.nanmax(cape),np.nanmin(cape))
    print(np.nanmax(subsat),np.nanmin(subsat))
    print(np.nanmax(wb),np.nanmin(wb))
    print(np.nanmax(wl),np.nanmin(wl))
    
    
    dts=dt.datetime.strptime(str(d1), "%Y%m")
    
    print('SAVING FILE')

    fout='/glade/p/univ/p35681102/fiaz/erai_data/buoy_erai_trmm/'+'bint_precip.'+d1+'.nc'

    ##### SAVE FILE ######

    try:ncfile.close()
    except:pass

    ncfile = Dataset(fout, mode='w', format='NETCDF4')

    ncfile.createDimension('time',prc.shape[0])
    ncfile.createDimension('lat',lat.size)
    ncfile.createDimension('lon',lon.size)

    ncfile.description="All points with surface pressure < 600 hPa are Nans"


    dy = ncfile.createVariable('time','i4',('time'))
    lt = ncfile.createVariable('lat',dtype('float32').char,('lat'))
    ln = ncfile.createVariable('lon',dtype('float32').char,('lon'))

    dy.units = "hours since"+str(dts)+";"
    dy.long_name = "initial time" ;
    dy.calendar = "gregorian"

    by = ncfile.createVariable('B_int',dtype('float32').char,('time','lat','lon'),zlib=True)
    cpe = ncfile.createVariable('cape',dtype('float32').char,('time','lat','lon'),zlib=True)
    sst = ncfile.createVariable('subsat',dtype('float32').char,('time','lat','lon'),zlib=True)
    wtb = ncfile.createVariable('wb',dtype('float32').char,('time','lat','lon'),zlib=True)
    prcp = ncfile.createVariable('precip',dtype('float32').char,('time','lat','lon'),zlib=True)

    by.description="B_int (in m/s^2)" 
    cpe.description="Undilute buoyancy--CAPE-like (in K)"
    sst.description="Subsaturation--effects of plume dilution (in K)"
    wtb.description="Fractional weighting of the CAPE-like component"
    prcp.description="TRMM 3B42 precip (in mm/hr)"

    dy[:]=np.asarray(prc.shape[0])*6.0
    lt[:]=lat
    ln[:]=lon

    by[:]=bint
    cpe[:]=cape
    sst[:]=subsat
    wtb[:]=wb
    prcp[:]=prc

    ncfile.close()
    print('FILE WRITTEN')
    
   
    



    
