import datetime

import numpy as np

import mygis
from bunch import Bunch

def write_file(date,info,erai):
    """writes ERAi input data to a netcdf file"""

    filename=str(date).replace(" ","_")
    dims    = ("time", "level","lat","lon")
    dims2dt = ("time", "lat","lon")

    extra_vars=[]
    # 3D variables
    # cloud,ice,qv,u,v,t,p
    # 2D variables
    # hgt,latent_heat,PBL_height,sensible_heat,sfc_hgt (sfc_hgt not used currently should be ~the same as hgt)
    # atts=Bunch(long_name="Cloud liquid water content",units="kg kg**-1", coordinates='latitude longitude')
    # extra_vars.append(Bunch(name="cloud",data=erai["cloud"],dims=dims,dtype="f",attributes=atts))
    #
    # atts=Bunch(long_name="Cloud ice water content",units="kg kg**-1", coordinates='latitude longitude')
    # extra_vars.append(Bunch(name="ice",data=erai["ice"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Relative Humidity",units="[]", coordinates='latitude longitude z time')
    extra_vars.append(Bunch(name="rh",data=erai["rh"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="U (E/W) wind speed",units="m s**-1", coordinates='latitude longitude z time')
    extra_vars.append(Bunch(name="u",data=erai["u"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="V (N/S) wind speed",units="m s**-1", coordinates='latitude longitude z time')
    extra_vars.append(Bunch(name="v",data=erai["v"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Potential Temperature",units="kg kg**-1", coordinates='latitude longitude z time')
    extra_vars.append(Bunch(name="theta",data=erai["t"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Pressure",units="Pa", coordinates='latitude longitude z time')
    extra_vars.append(Bunch(name="p",data=erai["p"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Atmospheric Elevation",units="m", coordinates='latitude longitude', axis="Z")
    extra_vars.append(Bunch(name="z",data=erai["z"],dims=dims,dtype="f",attributes=atts))

    atts=Bunch(long_name="Topographic Height",units="m", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="hgt",data=erai["sfc_hgt"],dims=dims[2:],dtype="f",attributes=atts))

    atts=Bunch(long_name="Total Precipitation",units="kg m**-2", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="pcp",data=erai["precip_total"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Convective Precipitation",units="kg m**-2", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="pcp_conv",data=erai["precip_conv"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface solar radiation (downwards)",units="W m**-2", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="swdown",data=erai["sw"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface longwave radiation (downwards)",units="W m**-2", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="lwdown",data=erai["lw"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Latent Heat flux (positive up)",units="W m**-2", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="latent_heat",data=erai["latent_heat"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Surface Sensible Heat flux (positive up)",units="W m**-2", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="sensible_heat",data=erai["sensible_heat"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Planetary Boundary Layer Height",units="m", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="PBL_height",data=erai["PBL_height"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Skin Temperature",units="K", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="tskin",data=erai["tskin"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Minimum 2m Temperature",units="K", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="t2min",data=erai["tmin"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="Maximum 2m Temperature",units="K", coordinates='latitude longitude time')
    extra_vars.append(Bunch(name="t2max",data=erai["tmax"],dims=dims2dt,dtype="f",attributes=atts))

    atts=Bunch(long_name="latitude",units="degrees", axis='Y')
    extra_vars.append(Bunch(name="latitude",data=info.lat_data,dims=dims[2:],dtype="f",attributes=atts))

    atts=Bunch(long_name="longitude",units="degrees", axis='X')
    extra_vars.append(Bunch(name="longitude",data=info.lon_data,dims=dims[2:],dtype="f",attributes=atts))

    time_since_1900 = date - datetime.datetime(1900,1,1,0,0,0)
    time = time_since_1900.days + np.float64(time_since_1900.seconds/86400.0)
    atts=Bunch(long_name="time",units="days since 1900-01-01",calendar='gregorian', axis="T")
    extra_vars.append(Bunch(name="time",data=time,dims=(dims[0],),dtype="d",attributes=atts))

    # for e in extra_vars:
    #     print(e.name, e.data.shape, e.dims)

    qvatts=Bunch(long_name="Specific Humidity",units="kg kg**-1", coordinates='latitude longitude z time')

    # write to output file
    mygis.write(filename=filename,varname="qv",data=erai.qv,attributes=qvatts,dtype="f",dims=dims,
                  extravars=extra_vars,history=" Produced by erai2gard "+info.version)
