import numpy as np
import scipy.io
from numba import jit,autojit
from scipy.interpolate import NearestNDInterpolator


def generate_region_mask(region_mask_filename, lat, lon):
    
    print("   Generating region mask..."),

    # Load & Pre-process Region Mask
    matfile=scipy.io.loadmat(region_mask_filename)
    lat_m=matfile["lat"]
    lon_m=matfile["lon"] # 0.125~359.875 deg
    region=matfile["region"]
    lon_m=np.append(lon_m,np.reshape(lon_m[0,:],(-1,1))+360,0)
    lon_m=np.append(np.reshape(lon_m[-2,:],(-1,1))-360,lon_m,0)
    region=np.append(region,np.reshape(region[0,:],(-1,lat_m.size)),0)
    region=np.append(np.reshape(region[-2,:],(-1,lat_m.size)),region,0)

    LAT,LON=np.meshgrid(lat_m,lon_m,sparse=False,indexing="xy")
    LAT=np.reshape(LAT,(-1,1))
    LON=np.reshape(LON,(-1,1))
    REGION=np.reshape(region,(-1,1))

    LATLON=np.squeeze(np.array((LAT,LON)))
    LATLON=LATLON.transpose()

    regMaskInterpolator=NearestNDInterpolator(LATLON,REGION)

    # Interpolate Region Mask onto Model Grid using Nearest Grid Value
#     pr_netcdf=Dataset(model_netcdf_filename,"r")
#     lon=np.asarray(pr_netcdf.variables[lon_var][:],dtype="float")
#     lat=np.asarray(pr_netcdf.variables[lat_var][:],dtype="float")
#     pr_netcdf.close()
    if lon[lon<0.0].size>0:
        lon[lon[lon<0.0]]+=360.0
    lat=lat[np.logical_and(lat>=-20.0,lat<=20.0)]

    LAT,LON=np.meshgrid(lat,lon,sparse=False,indexing="xy")
    LAT=np.reshape(LAT,(-1,1))
    LON=np.reshape(LON,(-1,1))
    LATLON=np.squeeze(np.array((LAT,LON)))
    LATLON=LATLON.transpose()
    REGION=np.zeros(LAT.size)
    for latlon_idx in np.arange(REGION.shape[0]):
        REGION[latlon_idx]=regMaskInterpolator(LATLON[latlon_idx,:])
    REGION=np.reshape(REGION.astype(int),(-1,lat.size))
    
    print("...Generated!")

    return REGION
  
@jit(nopython=True)   
def convecTransLev2_binThetae(lon_idx, REGION, NUMBER_CAPE_BIN, NUMBER_SUBSAT_BIN, 
NUMBER_BINT_BIN, CAPE, SUBSAT, BINT, RAIN, p0, p1, p2, q0, q1, q2):
 
    for lat_idx in np.arange(SUBSAT.shape[1]):
        subsat_idx=SUBSAT[:,lat_idx,lon_idx]
        cape_idx=CAPE[:,lat_idx,lon_idx]
        bint_idx=BINT[:,lat_idx,lon_idx]
        rain=RAIN[:,lat_idx,lon_idx]
        reg=REGION[lon_idx,lat_idx]
        if reg>0:
            for time_idx in np.arange(SUBSAT.shape[0]):
                if (cape_idx[time_idx]<NUMBER_CAPE_BIN and cape_idx[time_idx]>=0 
                and subsat_idx[time_idx]<NUMBER_SUBSAT_BIN and subsat_idx[time_idx]>=0):
                    p0[subsat_idx[time_idx],cape_idx[time_idx]]+=1
                    p1[subsat_idx[time_idx],cape_idx[time_idx]]+=rain[time_idx]
                    p2[subsat_idx[time_idx],cape_idx[time_idx]]+=rain[time_idx]**2

                if (bint_idx[time_idx]<NUMBER_BINT_BIN and bint_idx[time_idx]>=0):
                    q0[bint_idx[time_idx]]+=1
                    q1[bint_idx[time_idx]]+=rain[time_idx]
                    q2[bint_idx[time_idx]]+=rain[time_idx]**2
                
