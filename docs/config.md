# Configuring GARD

After building GARD, it is run on the command line following this syntax:

```
./gard downscale_options.txt
```

where `downscale_options.txt` is a Fortran style namelist (see example in GARD/run).

The namelist options are described in the tables below:

## Parameters

| Name                            | Type    | Required? | Default              | Description                                                                      |
|---------------------------------|---------|-----------|----------------------|----------------------------------------------------------------------------------|
| training_file                   | string  | No        | Options File Name    | Used to specify a separate namelist file for the training section shown below    |
| observation_file                | string  | No        | Options File Name    | Used to specify a separate namelist file for the observation section shown below |
| prediction_file                 | string  | No        | Options File Name    | Used to specify a separate namelist file for the prediction section shown below  |
| output_file                     | string  | No        | downscaled_output.nc | Downscaled output filename.                                                      |
| pass_through_var                | integer | No        | 1                    | Select the variable to be passed through as the output if pass_through=True      |
| n_analogs                       | integer | No (1)    | -1                   | The number of analogs to include                                                 |
| n_log_analogs                   | integer | No (1)    | -1                   | The number of analogs to include when computing threshold exceedence from analogs|
| analog_threshold                | real    | No (1)    | -1                   | If specified, GARD will compute the probability of exceeding this threshold      |
| start_date                      | string  | Yes       | n/a                  | Start date for prediction period                                                 |
| end_date                        | string  | Yes       | n/a                  | End date for prediction period                                                   |
| start_train                     | string  | Yes       | n/a                  | Start date for training period                                                   |
| end_train                       | string  | Yes       | n/a                  | End date for training period                                                     |
| start_transform                 | string  | Yes       | n/a                  | Start date for transform period                                                  |
| end_transform                   | string  | Yes       | n/a                  | End date for transform period                                                    |
| pure_regression                 | logical | No        | FALSE                | Use pure regression downscaling approach                                         |
| pure_analog                     | logical | No        | FALSE                | Use pure analog downscaling approach                                             |
| analog_regression               | logical | No        | TRUE                 | Use analog regression downscaling approach                                       |
| pass_through                    | logical | No        | False                | Pass a given predictor variable through as the output without downscaling        |
| sample_analog                   | logical | No        | FALSE                |                                                                                  |
| logistic_from_analog_exceedance | logical | No        | FALSE                |                                                                                  |
| logistic_threshold              | logical | No        | -9999                |                                                                                  |
| weight_analogs                  | logical | No        | TRUE                 |                                                                                  |
| debug                           | logical | No        | TRUE                 | prints more output at runtime and outputs files including the coefficients used in each analog regression (or analog values) as well as the predictor data |
| interactive                     | logical | No        | TRUE                 | Print downscaling status as a percentage on the command line                     |

Notes:

1.  Must specify one of n_analogs, n_log_analogs, or analog_threshold

## Training Parameters

See description of [training data](Datasets/data.md) suggestions.

| Name                  | Type    | Required? | Default | Description                                                                              |
|-----------------------|---------|-----------|---------|------------------------------------------------------------------------------------------|
| name                  | string  | Yes       | n/a     | Name of the training parameters dataset                                                  |
| preloaded             | string  | No        | n/a     | filepath of the preloaded training parameters dataset                                    |
| interpolation_method  | integer | No        | 1       | nearest neighbor= 1, bilenear =2                                                         |
| normalization_method  | integer | No        | 0       | no normalization = 0, mean/stddev from: training data = 1                                |
| time_indices          | integer | Yes       | -1      | list of timesteps in file to aggregate over (GEFS only)                                  |
| time_weights          | real    | No        | 1       | list of averaging weights for individual time indices
| agg_method            | integer | No        | 0       | per variable aggregation method when aggregating over time_indices: mean = 0, minimum = 1, maximum = 2, sum = 3 (GEFS only) |
| nvars                 | integer | Yes       | -1      | number of variables to be used in training                                               |
| data_type             | string  | Yes       | n/a     | dataset type: GEFS or GCM                                                                |
| lat_name              | string  | Yes       | n/a     | netCDF variable name for latitude                                                        |
| lon_name              | string  | Yes       | n/a     | netCDF variable name for longitude                                                       |
| time_name             | string  | Yes       | n/a     | netCDF variable name for time                                                            |
| nfiles                | integer | Yes       | -1      | number of files in each file list                                                        |
| input_transformations | integer | No        | 0       | no transform = 0, quantile mapping = 1, log transform = 2, cube root = 3, fifth root = 4 |
| var_names             | string  | Yes       | n/a     | variables names to use in training (one for each variable)                               |
| file_list             | string  | Yes       | n/a     | path to file containing a list of training files (one for each variable)                 |
| selected_time         | integer | No        | -1      | if set, only this time step will be read from each input file (GEFS only)                |
| calendar              | string  | Yes       | standard|                                                                                          |
| calendar_start_year   | integer | No        | 1900    |                                                                                          |
| timezone_offset       | real    | No        | 0       | hours to offset time variable (e.g. `time_name`) to account for timezone.                |


## Prediction Parameters

See description of [predictor data](Datasets/data.md) suggestions.

| Name                  | Type    | Required? | Default | Description                                                                              |
|-----------------------|---------|-----------|---------|------------------------------------------------------------------------------------------|
| name                  | string  | Yes       | n/a     | Name of the prediction parameters dataset                                                |
| preloaded             | string  | No        | n/a     | filepath of the preloaded prediction parameters dataset                                  |
| interpolation_method  | integer | No        | 1       | nearest neighbor= 1, bilenear =2                                                         |
| normalization_method  | integer | No        | 0       | no normalization = 0, mean/stddev from: prediction data = 1, training data = 2           |
| time_indices          | integer | Yes       | -1      | list of timesteps in file to aggregate over (GEFS only)                                  |
| time_weights          | real    | No        | 1       | list of averaging weights for individual time indices                                    |
| agg_method            | integer | No        | 0       | per variable aggregation method when aggregating over time_indices: mean = 0, minimum = 1, maximum = 2, sum = 3 (GEFS only) |
| nvars                 | integer | Yes       | -1      | number of prediction parameters to use in downscaling                                    |
| data_type             | string  | Yes       | n/a     | dataset type: GEFS or GCM                                                                |
| lat_name              | string  | Yes       | n/a     | netCDF variable name for latitude                                                        |
| lon_name              | string  | Yes       | n/a     | netCDF variable name for longitude                                                       |
| time_name             | string  | Yes       | n/a     | netCDF variable name for time                                                            |
| nfiles                | integer | Yes       | -1      | number of files in each file list                                                        |
| transformations       | integer | No        | 0       | no transform = 0, quantile mapping = 1, log transform = 2, cube root = 3, fifth root = 4 |
| input_transformations | integer | No        | 0       | no transform = 0, quantile mapping = 1, log transform = 2, cube root = 3, fifth root = 4 |
| var_names             | string  | Yes       | n/a     | variables names to use in prediction (one for each variable)                             |
| file_list             | string  | Yes       | n/a     | path to file containing a list of prediction filepaths (one for each variable)           |
| selected_time         | integer | No        | -1      | if set, only this time step will be read from each input file (GEFS only)                |
| calendar              | string  | Yes       | standard|                                                                                          |
| calendar_start_year   | integer | No        | 1900    |                                                                                          |
| timezone_offset       | real    | No        | 0       | hours to offset time variable (e.g. `time_name`) to account for timezone.                |

## Observation Parameters

See description of [observational data](Datasets/data.md) suggestions.

| Name                  | Type    | Required? | Default       | Description                                                                              |
|-----------------------|---------|-----------|---------------|------------------------------------------------------------------------------------------|
| name                  | string  | Yes       | n/a           | Name of the observation parameters dataset                                               |
| preloaded             | string  | No        | n/a           | filepath of the preloaded observation parameters dataset                                 |
| nvars                 | integer | Yes       | -1            | number of observation variables to downscale (currently must be 1)                       |
| nfiles                | integer | Yes       | -1            | number of files in each file list                                                        |
| data_type             | string  | Yes       | n/a           | dataset type, typically "obs"                                                            |
| lat_name              | string  | Yes       | n/a           | netCDF variable name for latitude                                                        |
| lon_name              | string  | Yes       | n/a           | netCDF variable name for longitude                                                       |
| time_name             | string  | Yes       | n/a           | netCDF variable name for time                                                            |
| input_transformations | integer | No        | 0             | no transform = 0, quantile mapping = 1, log transform = 2, cube root = 3, fifth root = 4 |
| var_names             | string  | Yes       | n/a           | variables names to use in obs dataset                                                    |
| file_list             | string  | Yes       | n/a           | path to file containing a list of obs filepaths                                          |
| calendar              | string  | Yes       | standard      |                                                                                          |
| calendar_start_year   | integer | No        | 1900          |                                                                                          |
