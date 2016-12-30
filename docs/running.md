
## Requirements

GARD has the following dependencies:

1. A Fortran compiler. GARD has been used with the following compilers:
  - gfortran (version 4.9.2)
  - ifort (version 15.0.2)

  GARD has also been compiled using these compilers:
  - PGI (pgf)

1. LAPACK — Linear Algebra PACKage.
1. netCDF4 - Network Common Data Form.

*Note: GARD allocates memory to the stack. Users should set the "The maximum stack size" to "unlimited" prior to building/running GARD. `ulimit -s unlimited`*

## Building GARD

GARD is built using a standard `makefile`.  From the command line, simply run the following command:

```
cd GARD/src/
make
```

Depending on which computer you are running GARD on, you may need to edit one or more of the following variables:

- `FC`: Fortran compiler
- `NCDF_PATH`: path to netCDF installation (see `nc-config --prefix` )
- `LAPACK_PATH`: path to lapack installation
- `INSTALL_DIR`: path to install GARD executable
- `RM`: path to unix `rm` command
- `CP`: path to unix `cp` command

## Running GARD

After building GARD, it is run on the command line following this syntax:

```
./gard downscale_options.txt
```

## Useful commands
Use the following to generate a list of e.g. GEFS precipitation files for input.

    ls -1 gefs/2010/*/apcp_sfc_*_mean.nc | sed 's/*//g;s/$/"/g;s/^/"/g'>gefs_pr_file.txt
