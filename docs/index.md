# Generalized Analog Regression Downscaling (GARD) code

This code is designed to provide a simple statistical downscaling method relying on
regressions and statistical transformations from various inputs (e.g. precipitation,
humidity, wind, PCA, etc.) to various outputs (e.g. precipitation, temperature, etc.)

## Useful commands
Use the following to generate a list of e.g. GEFS precipitation files for input.

    ls -1 gefs/2010/*/apcp_sfc_*_mean.nc | sed 's/*//g;s/$/"/g;s/^/"/g'>gefs_pr_file.txt


## Requirements

GARD has the following dependencies:

1. A Fortran compiler. GARD has been used with the following compilers:
  - gfortran (version 4.9.2)
  - ifort (version 15.0.2)

  GARD has also been compiled using these compilers:
  - PGI (pgf)

1. LAPACK — Linear Algebra PACKage.
1. netCDF4 - Network Common Data Form.

## Developing
T.B.D

## Reference
T.B.D.
