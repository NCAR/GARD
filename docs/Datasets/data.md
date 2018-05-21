# Data

To run the En-GARD downscaling package, three primary sets of data required, observed "truth", training data, and predictor data. The relationship between the training data and the observations will be used to downscale the predictor data. One of the primary advantages of En-GARD over more traditional "Model Output Statistics" (MOS) downscaling and bias correction methods (such as LOCA or BCSD), is the ability to incorporate atmospheric circulation properties into the downscaling method, for example upper level winds. As such it is recommended that the training and predictor variables include variables other than just precipitation and temperature.  

## Observations

The observations you supply to GARD represent the "truth" that GARD is trying to match.  These observations can be from a few stations or from a gridded dataset such as the [Newman](http://dx.doi.org/doi:10.5065/D6TH8JR2) or [Maurer](http://www.engr.scu.edu/~emaurer/gridded_obs/index_gridded_obs.html) datasets.  Regardless of the source of the data, the files must be in netcdf format with time, latitude, and longitude, dimensions (in that order). Latitude and Longitude can be somewhat irregular (e.g. in a lambert conformal conic projection), but if latitude and longitude are 1D variables, GARD will assume them to define a regular rectangular grid. Input data can be broken into as many files as desired, GARD will take a list of files as input, and will assume that the files are specified in chronological order.

The output data will be provided on the same grid and units as the input observations.

### Example
    netcdf nldas_met_update.obs.daily.pr.1979 {
    dimensions:
    	lon = 462 ;
    	lat = 222 ;
    	time = UNLIMITED ; // (365 currently)
    variables:
    	float lon(lon) ;
    		lon:standard_name = "longitude" ;
    		lon:units = "degrees_east" ;
    	float lat(lat) ;
    		lat:standard_name = "latitude" ;
    		lat:units = "degrees_north" ;
    	double time(time) ;
    		time:standard_name = "time" ;
    		time:units = "days since 1940-01-01 00:00:00" ;
    		time:calendar = "standard" ;
    	float pr(time, lat, lon) ;
    		pr:units = "mm/d" ;
    		pr:FillValue = 1.e+20f ;

## Training Data

The training data should be similar to the predictor data, but should be synchronous in time with the observations.  Because predictor data often from from free running climate models, there is no expectation that the predictors are synchronous with the observations.  For this reason, a reanalysis dataset such as [ERA-interim](https://www.ecmwf.int/en/forecasts/datasets/archive-datasets/reanalysis-datasets/era-interim) or [CFSR](http://cfs.ncep.noaa.gov/cfsr/) is preferred. Such datasets can have a wide variety of variables available, and it is left to the user to decide which variables are most appropriate for their application.  While GARD can (optionally) match the statistical properties of the training data by quantile mapping the predictor data to the training data, it is best if the two datasets (training and predictor) are as similar as possible.  For this reason, it is not recommended that one use, e.g. the highest resolution reanalysis dataset, such as the new ~30 km [ERA5](https://www.ecmwf.int/en/forecasts/datasets/archive-datasets/reanalysis-datasets/era5) data, when downscaling a 1 to 2.5 degree climate model. A better result is likely to be achieved by first running a low resolution regional climate model with both the training and predictor datasets as input as is often done in the [CORDEX](http://www.cordex.org) project.  The outputs of this regional climate model will then be more consistent in both the training and the predictor datasets, making the results from GARD likely to be more accurate. GARD can also be used for weather forecasting using, e.g. the [Global Ensemble Forecast System (GEFS)](https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/global-ensemble-forecast-system-gefs) for both training and predictor data.

Regardless of the training data used, the data files should be netcdf with at least time, latitude, and longitude dimensions.  Because atmospheric data often come in three spatial dimensions, GARD can also select a specific level from a four dimensional data set (e.g. time, level, latitude, longitude).  Input data can be broken into as many files as desired, GARD will take a list of files as input, and will assume that the files are specified in chronological order.

### Example
    netcdf wrf_daily_1990_data {
    dimensions:
            lat = 71 ;
            lon = 107 ;
            time = 365 ;
    variables:
        float XLAT(lat, lon) ;
                XLAT:description = "LATITUDE, SOUTH IS NEGATIVE" ;
                XLAT:units = "degree_north" ;
        float XLONG(lat, lon) ;
                XLONG:description = "LONGITUDE, WEST IS NEGATIVE" ;
                XLONG:units = "degree_east" ;
        double XTIME(time) ;
                XTIME:units = "minutes since 1979-01-01 00:00:00" ;
                XTIME:description = "minutes since 1979-01-01 00:00:00" ;
                XTIME:calendar = "standard" ;
        float PREC_TOT(time, lat, lon) ;
                PREC_TOT:description = "ACCUMULATED CUMULUS PRECIPITATION OVER prec_acc_dt PERIODS OF TIME" ;
                PREC_TOT:units = "mm" ;
        float U(time, lat, lon) ;
                U:description = "eastward-wind component" ;
                U:units = "m s-1" ;
        float V(time, lat, lon) ;
                V:description = "northward-wind component" ;
                V:units = "m s-1" ;
        float T(time, lat, lon) ;
                T:description = "potential temperature" ;
                T:units = "K" ;
        float QVAPOR(time, lat, lon) ;
                QVAPOR:description = "Water vapor mixing ratio" ;
                QVAPOR:units = "kg kg-1" ;


## Predictor Data

The predictor data to be used should be as similar as possible to the training data (see training data description), and should provide the same variables (e.g. zonal winds at 500 hPa). These data typically come from global climate models from the Coupled Models Intercomparison Project (CMIP). CMIP5 data can be accessed through the [Earth System Grid Federation](https://esgf-node.llnl.gov/projects/cmip5/). Another approach is to used regional climate model output, either run on your own or from a large project such as the [Coordinated Regional Downscaling Experiment (CORDEX)](http://www.cordex.org).  For example in North America, the [North America CORDEX](https://www.earthsystemgrid.org/project/NA-CORDEX.html) project has data from a variety of simulations with a nearly 50 and 25 km grid spacing.

### Example

    GARD predictor data should look similar to the training data.

    The variable names and grid spacing do not need to be the same, but they should overlap.

    Even though variables names do not need to be the same, they should represent the same quantity.
