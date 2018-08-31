import numpy as np
from bunch import Bunch

R=8.3144621 # J/mol/K
cp=29.19 # J/mol/K   =1.012 J/g/K
g=9.81 # m/s^2

def subset_pressure_levels(output_data):
    levels = [0,6,11,15,20,26,28,31]

    output_data.u   = output_data.u[:,levels,:,:]
    output_data.v   = output_data.v[:,levels,:,:]
    output_data.p   = output_data.p[:,levels,:,:]
    output_data.z   = output_data.z[:,levels,:,:]
    output_data.t   = output_data.t[:,levels,:,:]
    output_data.qv  = output_data.qv[:,levels,:,:]
    output_data.rh  = output_data.rh[:,levels,:,:]



# gard_atm_var=["u","v","gph","t","qv","ln_p_sfc","cloud","ice","sigma"]
def convert_atm(data):
    output_data=Bunch()
    output_data.u   = data.u[np.newaxis, ::-1,::-1,:]           # m/s
    output_data.v   = data.v[np.newaxis, ::-1,::-1,:]           # m/s
    output_data.hgt = data.gph[-1,::-1,:]/g                     # (m^2/s^2) / (m/s^2) = m

    # calculate pressure in Pa from ln(sfc_press) and hybrid sigma coordinates
    output_data.p = np.zeros((output_data.u.shape))
    nlevels = len(data.Plev)
    for i in range(nlevels):
        output_data.p[0,i,:,:] = data.Plev[nlevels - i - 1]

    output_data.z     = data.gph[np.newaxis,::-1,::-1,:]/9.81           # m
    output_data.t     = data.t[np.newaxis,::-1,::-1,:]                  # K
    output_data.qv    = data.qv[np.newaxis,::-1,::-1,:]                 # kg/kg
    output_data.rh    = data.rh[np.newaxis,::-1,::-1,:]              # kg/kg

    subset_pressure_levels(output_data)

    return output_data

# gard_sfc_var=["sensible_heat","latent_heat","hgt_98","PBL_height"]
def convert_sfc(data):

    dt = 3.0 * 60.0 * 60.0
    output_data=Bunch()

    output_data.precip_total    = data.precip_total[np.newaxis,::-1,:]      # kg/m^2
    output_data.precip_conv     = data.precip_conv[np.newaxis,::-1,:]       # kg/m^2
    output_data.sensible_heat   = data.sensible_heat[np.newaxis,::-1,:]/dt  # W/m^2
    output_data.latent_heat     = data.latent_heat[np.newaxis,::-1,:]/dt    # W/m^2
    output_data.sfc_hgt         = data.hgt_98[::-1,:]/g                     # (m^2/s^2) / (m/s^2) = m
    output_data.PBL_height      = data.PBL_height[np.newaxis,::-1,:]        # m
    output_data.tskin           = data.tskin[np.newaxis,::-1,:]             # K
    output_data.tmax            = data.tmax[np.newaxis,::-1,:]              # K
    output_data.tmin            = data.tmin[np.newaxis,::-1,:]              # K
    output_data.sw              = data.sw[np.newaxis,::-1,:] / dt   # convert from Joules to W /m^2
    output_data.lw              = data.lw[np.newaxis,::-1,:] / dt   # convert from Joules to W /m^2

    return output_data

def era2gard(data):
    output_data=Bunch()
    atm=convert_atm(data.atm)
    sfc=convert_sfc(data.sfc)

    for k in atm.keys():
        output_data[k]=atm[k]
    for k in sfc.keys():
        output_data[k]=sfc[k]

    return output_data
