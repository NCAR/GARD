import xarray as xr


min_var_list = ["t2min"]
max_var_list = ["t2max"]
sum_var_list = ["pcp", "pcp_conv"]

def make_daily(info):

    ds = xr.open_mfdataset(str(info.start_date.year)+"*")

    data_vars = dict()

    for v in ds.variables:
        if "time" in ds[v].dims:
            if v=="time":
                pass
            else:
                grouped = ds[v].resample(time="1D")
                if v in min_var_list:
                    data_vars[v] = grouped.min(dim="time")
                elif v in max_var_list:
                    data_vars[v] = grouped.max(dim="time")
                elif v in sum_var_list:
                    data_vars[v] = grouped.sum(dim="time")
                else:
                    data_vars[v] = grouped.mean(dim="time")
        else:
            data_vars[v] = ds[v]

    newds = xr.Dataset(data_vars = data_vars)

    newds.to_netcdf(str(info.start_date.year)+"_daily.nc")
