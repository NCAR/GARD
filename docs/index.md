# Generalized Analog Regression Downscaling (GARD) code

This code is designed to provide a simple statistical downscaling method relying on
regressions and statistical transformations from various inputs (e.g. precipitation,
humidity, wind, PCA, etc.) to various outputs (e.g. precipitation, temperature, etc.).  

Components of the core algorithms in GARD are documented in Clark and Hay (2004) and Clark and Slater (2006); the application to climate downscaling with additional enhancements is currently being documented in a new paper in prep.

## Documentation

A very basic set of documentation on how to run and configure GARD is slowly being developed.  Please see specific documentation below and in the sidebar on the left.

- [Running GARD](running.md)
- [Configuration Options](config.md)

## Developing
- [Development notes](Development/code.md)

## Input Data
- [Data description](Datasets/data.md)

## Reference
Clark, M. P., & Hay, L. E. (2004). Use of medium-range numerical weather prediction model output to produce forecasts of streamflow. _Journal of Hydrometeorology_, **5**(1), 15–32. [doi:10.1175/1525-7541(2004)005<0015:UOMNWP>2.0.CO;2](http://doi.org/10.1175/1525-7541(2004)005<0015:UOMNWP>2.0.CO;2)

Clark, M. P., & Slater, A. G. (2006). Probabilistic quantitative precipitation estimation in complex terrain. _Journal of Hydrometeorology_, **7**(1), 3–22. [doi:10.1175/JHM474.1](http://doi.org/10.1175/JHM474.1)
