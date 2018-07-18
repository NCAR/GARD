#!/usr/bin/env python

import xarray as xr
import numpy as np
from scipy import stats

varname = "PRECT"
output_base = "output/gard_out_"

has_pop = True

mean_file = output_base+varname+".nc"
err_file = output_base+varname+"_errors.nc"
pop_file = output_base+varname+"_logistic.nc"

# scrfs_file = "conus_scrfs_150yrs.nc"
scrfs_file = "AR_scrfs.nc"
output_file = "processed_{}.nc".format(varname)


def main():

    # open the primary data files
    mean = xr.open_dataset(mean_file)
    err = xr.open_dataset(err_file)
    scrfs = xr.open_dataset(scrfs_file)

    # apply the transform that was used within GARD
    d1 = np.cbrt(mean[varname])

    # get the random variables into a format that is appropriate in case the shapes don't match
    rand_tmp = scrfs.p_rand_uniform[:d1.shape[0], :d1.shape[1], :d1.shape[2]]

    # if there is a probability of precipitation (PoP) or other threshold value specified
    # then we need to do some extra work to use that prediction
    if has_pop:
        # load the PoP dataset
        pop = xr.open_dataset(pop_file)
        popval = pop[varname+"_exceedence_probability"]

        # figure out where the probability was not exceeded
        noprecip = rand_tmp.values < (1-popval.values)
        # rescale the random numbers to account for the limited portion of the cdf they should occupy
        rand_tmp.values -= (1-popval.values)
        rand_tmp /= popval.values

    # convert uniform random numbers to normally distributed values
    rand_tmp.values[rand_tmp.values<=1e-20] = 1e-20
    errors = stats.norm.ppf(rand_tmp.values)

    # if you get a lot of errors from inifite numbers, uncomment these lines
    # errors[errors<-10] = -10
    # errors[errors>10] = 10

    # compute the resulting value by added the error term and retransforming (**3)
    result = (d1 + (err[varname+"_error"].values * errors))**3  # Changed 6/6/18 from pr_error to PRECT_error

    if has_pop:
        # the error term could set values that would otherwise have had precip to <0
        result.values[result.values<1e-10] = 1e-10
        # also mask any values that fell did not meet the threshold criteria
        result.values[noprecip] = 0

    # finally write the results to an output file
    result.to_netcdf(output_file)

if __name__ == '__main__':
    main()
