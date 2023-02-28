from dask_jobqueue import PBSCluster
from dask.distributed import Client
import time

cluster = PBSCluster(
    queue="casper",
    walltime="03:00:00",
    project="P48500028",
    memory="30GB",
    cores=1,
    processes=1,
)

cluster.scale(12)

client=Client(cluster)
time.sleep(30) # wait 30 seconds to give all dask workers time to populate
print('Cluster created and assigned to dask client')

import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import timedelta,datetime
from dateutil.relativedelta import relativedelta


#--------------------------------------------------------------------------
# Function for writing ERA5 data to yearly input files for GARD


era5pl_template = "/glade/collections/rda/data/ds633.0/e5.oper.an.pl/{yr}{mth}/e5.oper.an.pl.{var_id}.{dt}00_{dt}23.grb"
era5sfc_template = "/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/{yr}{mth}/e5.oper.an.sfc.{var_id}.{dt}00_{dt2}23.grb"

var_formats = {'v':"128_132_v.ll025uv",
               'u':"128_131_u.ll025uv",
               'w':"128_135_w.ll025sc",
               'q':"128_133_q.ll025sc",
               't':"128_130_t.ll025sc",
               'tcrw':"228_089_tcrw.ll025sc",
                }

# new variable names to match those from CESM LENS 2
newvarname = {'v':'V',
               'u':'U',
               'w':'W',
               'q':'Q',
               't':'T',
               'tcrw':'PRECT',
                }

def createERA5Dataset(yr_st,yr_end,varlist):
    
    m=0
    for var in varlist:
        var_files = []
        dt = datetime(yr_st,1,1)
        while dt.year <= yr_end:
            yr = dt.year
            mth = "%.02d"%dt.month
            if var in ('u','v','w','q','t'):
                var_files.append(era5pl_template.format(yr=yr,mth=mth,var_id=var_formats[var],
                                                        dt=dt.strftime('%Y%m%d')))
                dt = dt + timedelta(days=1)
                
            elif var == 'tcrw':
                dt2 = dt + relativedelta(months=+1) - timedelta(days=1)
                var_files.append(era5sfc_template.format(yr=yr,mth=mth,var_id=var_formats[var],
                                                         dt=dt.strftime('%Y%m%d'),dt2=dt2.strftime('%Y%m%d')))
                
                dt = dt + relativedelta(months=+1)
        
        vardata = xr.open_mfdataset(var_files,concat_dim='time',combine='nested',
                                 backend_kwargs={"indexpath":""},parallel=True).sel(latitude=slice(50,20),
                                               longitude=slice(360-120,360-60))[var]
        
        if var in ('u','v','w','q','t'):
            vardata = vardata.isel(isobaricInhPa=20).drop_vars('isobaricInhPa') #~450 mb
            vardata = vardata.resample(time='1D').mean()

        else:
            vardata = vardata.resample(time='1D').sum()
        
        
        if m==0:
            era5_ds = vardata.to_dataset()
            era5_ds = era5_ds.rename({var:newvarname[var]})
        else:
            era5_ds = era5_ds.assign(var=vardata)
            era5_ds = era5_ds.rename({'var':newvarname[var]})
        
        m+=1
    
    era5_ds = era5_ds.drop_vars(('number','step','surface'))
    
    if yr_st!=yr_end:
        outfile = '/glade/scratch/shartke/gard/era5/era5_daily_%d_%d.nc'%(yr_st,yr_end)
    else:
        outfile = '/glade/scratch/shartke/gard/era5/era5_daily_%d.nc'%yr_st
    era5_ds.to_netcdf(outfile)    
      


#--------------------------------------------------------------------------
# Function for writing CESM LENS2 data to decadal input files for GARD

        
cesm_template = "/glade/campaign/cgd/cesm/CESM2-LE/atm/proc/tseries/day_1/{var}/b.e21.B{scen}{forcing}.f09_g17.LE2-{styr}.0{ens}.cam.h{i}.{var}.{yr1}0101-{yr2}1231.nc"
scen="HIST"
f="cmip6"

def createCESM2Dataset(yr,styr,varlist,enslist):
    
    yr1 = yr-yr%10
    
    for e in enslist:
        m=0
        for var in varlist:
            if var in ('U','V','T','Q'):
                ds = xr.open_dataset(cesm_template.format(var=var,scen=scen,forcing=f,styr=styr,
                                                          ens="%.02d"%e,yr1=yr1,yr2=yr1+9,i=6))[var]
                # select data over CONUS at ~450 mb level
                vardata = ds.sel(lev=ds.lev[19],lat=slice(20,50),lon=slice(360.-120.,360.-60.)) # ,time=slice(str(yr),str(yr))
                vardata = vardata.drop_vars('lev')
            elif var in ('PSL','PRECT'):
                ds = xr.open_dataset(cesm_template.format(var=var,scen=scen,forcing=f,styr=styr,
                                                          ens="%.02d"%e,yr1=yr1,yr2=yr1+9,i=1))[var]
                vardata = ds.sel(lat=slice(20,50),lon=slice(360.-120.,360.-60.)) # ,time=slice(str(yr),str(yr))
                # convert m/s to mm/d
                vardata = vardata*3600*24*1000
            
            if m==0:
                cesm_ds = vardata.to_dataset()
            else:
                cesm_ds = cesm_ds.assign(var=vardata)
                cesm_ds = cesm_ds.rename({'var':var})
            
            m+=1
        
        outfile = '/glade/scratch/shartke/gard/cesmlens2/cesm_daily_%d_%d_%d_%.02d.nc'%(yr,yr+9,styr,e)
        cesm_ds.to_netcdf(outfile)    

        
#--------------------------------------------------------------------------
    
# Note: Generating the ERA5 datasets will take the bulk of the time for this program

print(datetime.now())
styr = 1301 # 1231, 1251, 1281, or 1301
createCESM2Dataset(1960,styr,['U','V','W','Q','T','PRECT'],np.arange(1,3))
createCESM2Dataset(1970,styr,['U','V','W','Q','T','PRECT'],np.arange(1,3))
print('CESM LENS2 datasets complete at: ',datetime.now())

for yr in (1980,1981,1982,1983):
    createERA5Dataset(yr,yr,['u','v','w','q','t','tcrw'])
    print('ERA5 %s dataset complete at: '%yr,datetime.now())




# now you should be able to train GARD using 1980-1999 ERA5 data
# and predict downscaled 1960-1979 precip or temp using CESM LENS2 data
