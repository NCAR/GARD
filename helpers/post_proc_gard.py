#!/usr/bin/env python
"""
SYNOPSIS

    post_proc_gard.py [-h] [--verbose] [-v, --version]
                    [-b base_filename]
                    [-o output_filename]
                    [-v variable_name]
                    [-s SCRF_filename]
                    [-t transform]
                    [--has_pop]
                    [--offset offset]

DESCRIPTION

    TODO This describes how to use this script.
    This docstring will be printed by the script if there is an error or
    if the user requests help (-h or --help).

EXAMPLES

    post_proc_gard.py -b output/gard_out_ -o processed_gard.nc -v pcp -t cuberoot --has_pop
    post_proc_gard.py -b output/gard_out_ -o processed_gard.nc -v tmin -t None

AUTHOR

    Ethan Gutmann - gutmann@ucar.edu

LICENSE

    This script is in the public domain.

VERSION


"""
from __future__ import absolute_import, print_function, division

import sys
import os
import traceback
import argparse

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import xarray as xr
import numpy as np
from scipy import stats

global verbose
verbose=False


def process_file(varname = "pcp",
                 output_base = "output/gard_out_",
                 has_pop = False,
                 scrfs_file = "conus_scrfs_150yrs.nc",
                 output_file = None,
                 scrf_offset = 0,
                 transform="cuberoot"):

    mean_file = output_base+varname+".nc"
    err_file = output_base+varname+"_errors.nc"
    pop_file = output_base+varname+"_logistic.nc"

    if output_file is None:
        output_file = "processed_{}.nc".format(varname)

    if verbose:
        if has_pop:
            print("Reading files:\n Mean:{}\n err:{}\n scrf:{}\n pop:{}".format(mean_file, err_file, scrfs_file, pop_file))
        else:
            print("Reading files:\n Mean:{}\n err:{}\n scrf:{}".format(mean_file, err_file, scrfs_file))
        print("Using transform: "+transform)
        print("Outputfile = "+output_file)
        print("varname = "+varname)


    if verbose: print("open the primary data files")
    mean = xr.open_dataset(mean_file)
    err = xr.open_dataset(err_file)
    scrfs = xr.open_dataset(scrfs_file)

    if verbose: print("apply the transform that was used within GARD")
    if transform=="cuberoot":
        d1 = np.cbrt(mean[varname])
    elif transform=="squareroot":
        d1 = np.sqrt(mean[varname])
    elif transform=="log":
        d1 = np.log(mean[varname])
    elif transform=="None":
        d1 = mean[varname]
    else:
        raise KeyError("Unknown transform:"+transform)



    if verbose: print("Get the random variables into a format that is appropriate in case the shapes don't match")
    rand_tmp = scrfs.p_rand_uniform[scrf_offset:scrf_offset+d1.shape[0], :d1.shape[1], :d1.shape[2]]

    # if there is a probability of precipitation (PoP) or other threshold value specified
    # then we need to do some extra work to use that prediction
    if has_pop:
        if verbose: print("load the PoP open_dataset")
        pop = xr.open_dataset(pop_file)
        popval = pop[varname+"_exceedence_probability"]

        if verbose: print("figure out where the probability was not exceeded")
        noprecip = rand_tmp.values < (1-popval.values)
        if verbose: print("rescale the random numbers to account for the limited portion of the cdf they should occupy")
        rand_tmp.values -= (1-popval.values)
        rand_tmp /= popval.values

        if verbose: print("convert uniform random numbers to normally distributed values")
        rand_tmp.values[rand_tmp.values<=1e-20] = 1e-20
        errors = stats.norm.ppf(rand_tmp.values)

        # if you get a lot of errors from inifite numbers, uncomment these lines
        # errors[errors<-10] = -10
        # errors[errors>10] = 10
    else:
        if verbose: print("Getting error term")
        errors = scrfs.t_rand[scrf_offset:scrf_offset+d1.shape[0], :d1.shape[1], :d1.shape[2]]

    if verbose: print("compute the resulting value by added the error term and retransforming")
    result = (d1 + (err[varname+"_error"].values * errors))

    if transform=="cuberoot":
        result = result**3
    elif transform=="squareroot":
        result = result**2
    elif transform=="log":
        result = np.exp(result)
    elif transform=="None":
        pass
    else:
        print("Unknown transform:"+transform)


    if has_pop:
        # the error term could set values that would otherwise have had precip to <0
        result.values[result.values<1e-10] = 1e-10
        if verbose: print("Mask any values that fell did not meet the threshold criteria")
        result.values[noprecip] = 0

    if verbose: print("Write the results to an output file")
    result.to_netcdf(output_file)

if __name__ == '__main__':
    try:
        parser= argparse.ArgumentParser(description='Post process GARD output',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument('-b', nargs="?", dest="base_filename", action='store',
                default= "output/gard_out_", help="root of GARD output filenames to read")
        parser.add_argument('-o', nargs="?", dest="output_filename", action='store',
                default= None, help="Name of output file to create")
        parser.add_argument('-v', nargs="?", dest="variable_name", action='store',
                default= "pcp", help="Name of variable name to process")
        parser.add_argument('-s', nargs="?", dest="SCRF_filename", action='store',
                default= "conus_scrfs_150yrs.nc", help="Name of SCRF file to read")
        parser.add_argument('-t', nargs="?", dest="transform", action='store',
                default="cuberoot", help="Transformation to apply [cuberoot, squareroot, log, None]")
        parser.add_argument('--offset', nargs="?", dest="offset", action='store',
                default="0", help="offset into SCRF fields to readily create ensembles")
        parser.add_argument('--has_pop', action='store_true',
                default=False, help='verbose output', dest='has_pop')

        parser.add_argument('--version',action='version',
                version='post_proc_gard 1.0')
        parser.add_argument ('--verbose', action='store_true',
                default=False, help='verbose output', dest='verbose')
        args = parser.parse_args()

        verbose = args.verbose

        exit_code = process_file(varname = args.variable_name,
                                 output_base = args.base_filename,
                                 has_pop = args.has_pop,
                                 scrfs_file = args.SCRF_filename,
                                 output_file = args.output_filename,
                                 scrf_offset = int(args.offset),
                                 transform = args.transform)
        if exit_code is None:
            exit_code = 0
        sys.exit(exit_code)
    except KeyboardInterrupt as e: # Ctrl-C
        raise e
    except SystemExit as e: # sys.exit()
        raise e
    except Exception as e:
        print('ERROR, UNEXPECTED EXCEPTION')
        print(str(e))
        traceback.print_exc()
        os._exit(1)
