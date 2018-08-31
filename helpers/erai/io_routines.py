import subprocess,os
import numpy as np
import mygis
from bunch import Bunch

global last_precip
last_precip = -1

sfcvarlist = ["TP_GDS4_SFC", "CP_GDS4_SFC", "2T_GDS4_SFC","SSHF_GDS4_SFC","SLHF_GDS4_SFC","Z_GDS4_SFC","BLH_GDS4_SFC","SSRD_GDS4_SFC","STRD_GDS4_SFC", "SKT_GDS4_SFC"]
gard_sfc_var = ["precip_total","precip_conv","t2m","sensible_heat","latent_heat","hgt_98","PBL_height","sw","lw", "tskin"]

atmvarlist = ["Z_GDS4_ISBL","T_GDS4_ISBL","Q_GDS4_ISBL","R_GDS4_ISBL","lv_ISBL0"]
gard_atm_var = ["gph","t","qv","rh","Plev"]

atmuvlist = ["U_GDS4_ISBL","V_GDS4_ISBL"]
gard_uv_var = ["u","v"]

converted_sfc_files = []
sfc_ncfiles = dict()

def grib2nc(erai_file,varlist,output_dir):
    """convert a grib file to a netcdf file"""
    outputfile=output_dir+erai_file.split("/")[-1]+".nc"
    if not os.path.isfile(outputfile):
        try:
            print("Converting: "+erai_file.split("/")[-1])
            print("ncl_convert2nc "+erai_file+" -e grb -L -v "+",".join(varlist)+" -o "+output_dir)
            os.system("ncl_convert2nc "+erai_file+" -e grb -L -v "+",".join(varlist)+" -o "+output_dir +"&> /dev/null")
        except:
            print("ERROR: ncl_convert2nc is not available in your $PATH (?)")
            print("     On many HPC systems you can run `module load ncl`.")
            print("     On other systems you may need to install NCL.")

    return outputfile

def find_sfc_file(time,info):
    file_base= info.sfcdir+info.sfcfile
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_M_","{0:02}".format(time.month))
    file_base= file_base.replace("_D_","{0:02}".format(time.day))
    if (time.hour>=12):
        hour="12"
        offset=round(time.hour/3)-4
    else:
        hour="00"
        offset=round(time.hour/3)
    file_base= file_base.replace("_h_",hour)
    return file_base,offset

def find_atm_file(time,info):
    file_base= info.atmdir+info.atmfile
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_M_","{0:02}".format(time.month))
    file_base= file_base.replace("_D_","{0:02}".format(time.day))
    atm_file = file_base.replace("_h_","{0:02}".format(time.hour))

    file_base= info.atmdir+info.uvfile
    file_base= file_base.replace("_Y_",str(time.year))
    file_base= file_base.replace("_M_","{0:02}".format(time.month))
    file_base= file_base.replace("_D_","{0:02}".format(time.day))
    uv_file  = file_base.replace("_h_","{0:02}".format(time.hour))

    return uv_file,atm_file

def load_sfc(time,info):
    """load surface forcing from a grib file (or netcdf file if it has been converted previously)"""
    global last_precip

    inputfile, offset = find_sfc_file(time,info)
    if last_precip > offset: last_precip = -1

    try:
        nc_file = sfc_ncfiles[inputfile]
    except KeyError:
        nc_file = grib2nc(inputfile, sfcvarlist, info.nc_file_dir)
        sfc_ncfiles[inputfile] = nc_file

    outputdata = Bunch()
    for s,v in zip(gard_sfc_var, sfcvarlist):
        # print(nc_file)
        nc_data = mygis.read_nc(nc_file,v,returnNCvar=True)
        input_data = nc_data.data[:,info.ymin:info.ymax,info.xmin:info.xmax]

        if ((s=="latent_heat") or (s=="sensible_heat") or (s=="sw") or (s=="lw")):
            # if offset>=3:
            #     input_data[offset,...]-=input_data[offset-3,...]
            #     input_data[offset,...]/3.0
            # if offset>=2:
            #     input_data[offset,...] -= input_data[offset-2,...]
            #     input_data[offset,...] /= 2
            # elif offset==1:
            if offset>0:
                input_data[offset,...] -= input_data[offset-1,...]

        if (v=="2T_GDS4_SFC"):
            outputdata["tmax"] = input_data[int(offset):int(offset+2),:,:].max(axis=0)
            outputdata["tmin"] = input_data[int(offset):int(offset+2),:,:].min(axis=0)

        if ((s=="precip_conv") or (s=="precip_total")):

            offset += 1 # needs to accumulate over the entire time period, assumes two time steps... :(
            if last_precip!=-1:
                input_data[offset,...] -= input_data[last_precip,...]

        outputdata[s] = input_data[int(offset),:,:]

        if ((s=="precip_conv") or (s=="precip_total")):
            offset -= 1 # needs to be reset for other variables

        nc_data.ncfile.close()

    last_precip = offset+1
    return outputdata

def load_atm(time,info):
    """Load atmospheric variable from a GRIB file"""
    uvfile, scfile = find_atm_file(time, info)
    uvnc_file = grib2nc(uvfile, atmuvlist, info.nc_file_dir)
    scnc_file = grib2nc(scfile, atmvarlist, info.nc_file_dir)

    outputdata = Bunch()
    for s,v in zip(gard_uv_var, atmuvlist):
        nc_data = mygis.read_nc(uvnc_file, v, returnNCvar=True)
        if len(nc_data.data.shape)==3:
            outputdata[s] = nc_data.data[:,info.ymin:info.ymax,info.xmin:info.xmax]
        else:
            outputdata[s] = nc_data.data[info.ymin:info.ymax,info.xmin:info.xmax]
        nc_data.ncfile.close()

    for s,v in zip(gard_atm_var, atmvarlist):
        nc_data = mygis.read_nc(scnc_file, v, returnNCvar=True)

        if len(nc_data.data.shape)==3:
            outputdata[s] = nc_data.data[:,info.ymin:info.ymax,info.xmin:info.xmax]
        elif len(nc_data.data.shape)==2:
            outputdata[s] = nc_data.data[info.ymin:info.ymax,info.xmin:info.xmax]
        elif len(nc_data.data.shape)==1:
            outputdata[s] = nc_data.data[:]
        else:
            try:
                outputdata[s] = nc_data.data[:]
            except:
                outputdata[s] = nc_data.data.get_value()

        nc_data.ncfile.close()

    return outputdata


def load_data(time,info):
    """docstring for load_data"""
    sfc = load_sfc(time, info)
    atm = load_atm(time, info)
    return Bunch(sfc=sfc, atm=atm)
