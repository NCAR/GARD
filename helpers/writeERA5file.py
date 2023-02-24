import os
import numpy as np
import pandas as pd
import xarray as xr
#import qfgrib
from datetime import timedelta,datetime
from dateutil.relativedelta import relativedelta


era5pl_template = "/glade/collections/rda/data/ds633.0/e5.oper.an.pl/{yr}{mth}/e5.oper.an.pl.{var_id}.{dt}00_{dt}23.grb"
era5sfc_template = "/glade/collections/rda/data/ds633.0/e5.oper.an.sfc/{yr}{mth}/e5.oper.an.sfc.{var_id}.{dt}00_{dt2}23.grb"
w_format = "128_135_w.ll025sc"
v_format = "128_132_v.ll025uv"
u_format = "128_131_u.ll025uv"
q_format = "128_133_q.ll025sc"
t_format = "128_130_t.ll025sc"
tcrw_format = "228_089_tcrw.ll025sc"



var_formats = {'v':v_format,
               'u':u_format,
               'w':w_format,
               'q':q_format,
               't':t_format,
               'tcrw':tcrw_format,
                }

# new variable names to match those from CESM LENS 2
newvarname = {'v':'V',
               'u':'U',
               'w':'W',
               'q':'Q',
               't':'T',
               'tcrw':'PRECT',
                }

def createInputDataset(yr_st,yr_end,varlist):
    
    m=0
    for var in varlist:
        var_files = []
        dt = datetime(yr_st,1,1)
        while dt.year <= yr_end:
            yr = dt.year
            mth = "%.02d"%dt.month
            if var in ('u','v','w','q','t'):
                var_files.append(era5pl_template.format(yr=yr,mth=mth,var_id=var_formats[var],dt=dt.strftime('%Y%m%d')))
                dt = dt + timedelta(days=1)
                
            elif var == 'tcrw':
                dt2 = dt + relativedelta(months=+1) - timedelta(days=1)
                var_files.append(era5sfc_template.format(yr=yr,mth=mth,var_id=var_formats[var],
                                                         dt=dt.strftime('%Y%m%d'),dt2=dt2.strftime('%Y%m%d')))
                
                dt = dt + relativedelta(months=+1)
        
        vardata = xr.open_mfdataset(var_files,concat_dim='time',combine='nested',
                                 backend_kwargs={"indexpath":""}).sel(latitude=slice(50,20),
                                               longitude=slice(360-120,360-60))[var]
        
        if var in ('u','v','w','q','t'):
            vardata = vardata.isel(isobaricInhPa=20).drop_vars('isobaricInhPa') #~450 mb
            vardata = vardata.resample(time='1D').mean()

        else:
            vardata = vardata.resample(time='1D').sum()
        
        
        if m==0:
            era5_ds = vardata.to_dataset()
        else:
            era5_ds = era5_ds.assign(var=vardata)
        
        era5_ds = era5_ds.rename({'var':newvarname[var]})
        
        m+=1
    
    if yr_st!=yr_end:
        outfile = '/glade/scratch/shartke/gard/era5/era5_daily_%d_%d.nc'%(yr_st,yr_end)
    else:
        outfile = '/glade/scratch/shartke/gard/era5/era5_daily_%d.nc'%yr_st
    era5_ds.to_netcdf(outfile)    
      

    
print(datetime.now())
createInputDataset(1982,1982,['u','v','w','q','t','tcrw'])
print('1982 dataset complete at: ',datetime.now())
createInputDataset(1983,1983,['u','v','w','q','t','tcrw'])
print('1983 dataset complete at: ',datetime.now())