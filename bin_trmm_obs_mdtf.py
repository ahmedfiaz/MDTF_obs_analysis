'''
PURPOSE: To bin the precipitation vs. Bint and a two-dimensional binnning
         using observed precipitation data and the same code as used for
         MDTF.

AUTHOR: Fiaz Ahmed  

DATE: 01/13/19 

'''

import datetime as dt
import numpy as np
from netCDF4 import Dataset
from dateutil.relativedelta import relativedelta
import glob,itertools
from sys import exit
#import h5py
# from bin_parameters import * ## Import all the region indices and relevant variables
# import bin_cython1
from dateutil.relativedelta import relativedelta
import time
# from mpi4py import MPI
import pickle
from mdtf_binning_funcs import generate_region_mask,convecTransLev2_binThetae
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

print 'STARTING SCRIPT'

#from ncdump import ncdump

### 0=ocean, 1=land, 2=lake, 3=small island, 4=ice shelf

##### Set Bin Parameters here ######
BINT_BIN_WIDTH=0.01 # default=0.3 (following satellite retrieval product)
BINT_RANGE_MAX=1.51 # default=90 (75 for satellite retrieval product)
BINT_RANGE_MIN=-1.5 # default=90 (75 for satellite retrieval product)

# Bin width and intervals for CAPE and SUBSAT.
# In units of K
CAPE_RANGE_MIN=-40.0
CAPE_RANGE_MAX=17.0
CAPE_BIN_WIDTH=1.0

SUBSAT_RANGE_MIN=0.0
SUBSAT_RANGE_MAX=42.0
SUBSAT_BIN_WIDTH=1.0            
#############################
cape_bin_center=np.arange(CAPE_RANGE_MIN,CAPE_RANGE_MAX+CAPE_BIN_WIDTH,CAPE_BIN_WIDTH)
subsat_bin_center=np.arange(SUBSAT_RANGE_MIN,SUBSAT_RANGE_MAX+SUBSAT_BIN_WIDTH,SUBSAT_BIN_WIDTH)
bint_bin_center=np.arange(BINT_RANGE_MIN,BINT_RANGE_MAX+BINT_BIN_WIDTH,BINT_BIN_WIDTH)

NUMBER_CAPE_BIN=cape_bin_center.size
NUMBER_SUBSAT_BIN=subsat_bin_center.size
NUMBER_BINT_BIN=bint_bin_center.size

#############################

# Allocate Memory for Arrays
P0=np.zeros((NUMBER_SUBSAT_BIN,NUMBER_CAPE_BIN))
P1=np.zeros((NUMBER_SUBSAT_BIN,NUMBER_CAPE_BIN))
P2=np.zeros((NUMBER_SUBSAT_BIN,NUMBER_CAPE_BIN))
PE=np.zeros((NUMBER_SUBSAT_BIN,NUMBER_CAPE_BIN))

Q0=np.zeros((NUMBER_BINT_BIN))
Q1=np.zeros((NUMBER_BINT_BIN))
Q2=np.zeros((NUMBER_BINT_BIN))


dts=[]

list1=[]
list2=[]
list_trmm=[]

## Load paths and file name information ##
# dirc='/glade/scratch/fiaz/regridded/buoy_pbl/'

dirc='/glade/p/univ/p35681102/fiaz/erai_data/buoy_erai_trmm/'
fil='bint_precip.'

 ##############################
 
dts=[]

strt_date=dt.datetime(2002,1,1)
end_date=dt.datetime(2014,12,31)

date_str=str(strt_date.timetuple().tm_year)+'_'+str(end_date.timetuple().tm_year)

BIN_OUTPUT_DIR='/glade/u/home/fiaz/buoy_maps/files/'
BIN_OUTPUT_FILENAME='trmm3B42_erai_'+date_str+".convecTransLev2"

while strt_date<=end_date:
    d1=strt_date.strftime("%Y%m")
    fname=dirc+fil+str(d1)+'*.nc'
    list_temp=(glob.glob(fname))
    list_temp.sort()
    list_trmm.append(list_temp)
    strt_date+=relativedelta(months=1)

chain1=itertools.chain.from_iterable(list_trmm)
list1= (list(chain1))

jobs=list1

f=Dataset(list1[0],'r')
lat,lon=f.variables['lat'][:],f.variables['lon'][:]
f.close()

ind_lat=np.where(abs(lat)<20)[0]

region_mask_dir='/glade/u/home/fiaz/buoy_maps/files/'
region_mask_file='region_0.25x0.25_costal2.5degExcluded.mat'

REGION=generate_region_mask(region_mask_dir+region_mask_file, lat, lon)


for j in jobs:

    print j[-9:-3]

    date=dt.datetime.strptime(j[-9:-3],"%Y%m")
    d1=date.strftime("%Y-%m")
    d2=date.strftime("%Y%m")
    month=date.timetuple().tm_mon


    ##### LOAD BINNING VARIABLES #####

    f=Dataset(j,"r")    
    Bint_era=np.asarray(f.variables['B_int'][:,ind_lat,:])
    cape_era=np.asarray(f.variables['cape'][:,ind_lat,:])
    subsat_era=np.asarray(f.variables['subsat'][:,ind_lat,:])
    prc_trmm=np.asarray(f.variables['precip'][:,ind_lat,:])
    f.close()
                  
    ### Compute a value of integrated buoyancy

    ### Start binning
    SUBSAT=(subsat_era-SUBSAT_RANGE_MIN)/SUBSAT_BIN_WIDTH-0.5
    SUBSAT=SUBSAT.astype(int)
    
    CAPE=(cape_era-CAPE_RANGE_MIN)/CAPE_BIN_WIDTH-0.5
    CAPE=CAPE.astype(int)
    
    BINT=(Bint_era-BINT_RANGE_MIN)/(BINT_BIN_WIDTH)+0.5
    BINT=BINT.astype(int)

    
#     print(CAPE[0,0,0],cape_era[0,0,0])

    RAIN=prc_trmm       
    RAIN[RAIN<0]=0 # Sometimes models produce negative rain rates
        
        # Binning is structured in the following way to avoid potential round-off issue
    #  (an issue arise when the total number of events reaches about 1e+8)
    for lon_idx in np.arange(SUBSAT.shape[2]):
        p0=np.zeros((NUMBER_SUBSAT_BIN,NUMBER_CAPE_BIN))
        p1=np.zeros((NUMBER_SUBSAT_BIN,NUMBER_CAPE_BIN))
        p2=np.zeros((NUMBER_SUBSAT_BIN,NUMBER_CAPE_BIN))
        pe=np.zeros((NUMBER_SUBSAT_BIN,NUMBER_CAPE_BIN))

        q0=np.zeros((NUMBER_BINT_BIN))
        q1=np.zeros((NUMBER_BINT_BIN))
        q2=np.zeros((NUMBER_BINT_BIN))
    
        convecTransLev2_binThetae(lon_idx, REGION, NUMBER_CAPE_BIN, NUMBER_SUBSAT_BIN, 
        NUMBER_BINT_BIN, CAPE, SUBSAT, BINT, RAIN, p0, p1, p2, q0, q1, q2)

        P0+=p0
        P1+=p1
        P2+=p2
    
        Q0+=q0
        Q1+=q1
        Q2+=q2


    print("...Complete for current files!")

print("   Total binning complete!")
    
    
# Save Binning Results
bin_output_netcdf=Dataset(BIN_OUTPUT_DIR+BIN_OUTPUT_FILENAME+".nc","w",format="NETCDF4")
        
bin_output_netcdf.description="Convective Onset Buoyancy Statistics for "+"TRMM 3B42-ERAI"
bin_output_netcdf.source="Convective Onset Buoyancy Statistics Diagnostic Package \
- as part of the NOAA Model Diagnostic Task Force (MDTF) effort"

subsat_dim=bin_output_netcdf.createDimension("subsat",len(subsat_bin_center))
subsat_var=bin_output_netcdf.createVariable("subsat",np.float64,("subsat",))
subsat_var.units="K"
subsat_var[:]=subsat_bin_center

cape_dim=bin_output_netcdf.createDimension("cape",len(cape_bin_center))
cape=bin_output_netcdf.createVariable("cape",np.float64,("cape"))
cape.units="K"
cape[:]=cape_bin_center

bint_dim=bin_output_netcdf.createDimension("bint",len(bint_bin_center))
bint_var=bin_output_netcdf.createVariable("bint",np.float64,("bint",))
bint_var.units="K"
bint_var[:]=bint_bin_center


p0=bin_output_netcdf.createVariable("P0",np.float64,("subsat","cape"))
print(p0.shape,P0.shape)
p0[:,:]=P0

p1=bin_output_netcdf.createVariable("P1",np.float64,("subsat","cape"))
p1.units="mm/hr"
p1[:,:]=P1

p2=bin_output_netcdf.createVariable("P2",np.float64,("subsat","cape"))
p2.units="mm^2/hr^2"
p2[:,:]=P2


q0=bin_output_netcdf.createVariable("Q0",np.float64,("bint"))
q0[:]=Q0

q1=bin_output_netcdf.createVariable("Q1",np.float64,("bint"))
q1.units="mm/hr"
q1[:]=Q1

q2=bin_output_netcdf.createVariable("Q2",np.float64,("bint"))
q2.units="mm^2/hr^2"
q2[:]=Q2


bin_output_netcdf.close()

print("   Binned results saved as "+BIN_OUTPUT_DIR+BIN_OUTPUT_FILENAME+".nc!")





    
